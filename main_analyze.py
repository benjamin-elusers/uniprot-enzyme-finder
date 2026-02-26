# main_analyze.py
import argparse
import json
import time
import gzip
from collections import defaultdict, Counter
from tqdm import tqdm

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.csv as pv
from config.settings import FULL_VERSION
from config.settings import FINAL_OUTPUT_FILE, STATS_FILE

def get_headers(filepath):
    """Reads just the first line to get the exact column names."""
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        header_line = f.readline().strip('\n')
        return header_line.split('\t')

def main():
    parser = argparse.ArgumentParser(description="Stage 3: Generate Coverage and Stats Reports")
    parser.add_argument("--input", "-i", required=True, help="Input final TSV.gz file")
    parser.add_argument("--out_stats", default=f"uniprot_stats_summary_{FULL_VERSION}.json", help="Output JSON stats file")
    parser.add_argument("--out_cov", default=f"uniprot_field_coverage_{FULL_VERSION}.tsv", help="Output TSV coverage file")
    parser.add_argument("-t", "--total", type=int, default=None, help="Total rows (for exact ETA progress bar)")
    args = parser.parse_args()

    # =========================================================================
    #   DATA STRUCTURES
    # =========================================================================
    
    coverage_counts = Counter()
    
    stats = {
        "Total_Entries": 0,
        "Reviewed_Count": 0,
        "Unreviewed_Count": 0,
        "Classes": defaultdict(lambda: {
            "Total": 0,
            "Reviewed": 0,
            "Confidence": Counter(),
            "Enzyme_Type": Counter(),
            "Functional_Category": Counter(),
            "EC_Completeness": Counter(),
            "Evidence_Flags": Counter()
        })
    }

    print(f"[Start] Analyzing {args.input} (C++ Streaming Mode)...")
    start_time = time.time()

    # =========================================================================
    #   PYARROW CHUNKED STREAMING (Low Memory, High Speed)
    # =========================================================================
    headers = get_headers(args.input)
    
    # Force all columns to strings to prevent inference crashes on sparse data
    column_types = {col: pa.string() for col in headers}
    
    parse_options = pv.ParseOptions(delimiter='\t')
    # Stream in 50MB chunks to keep RAM usage extremely low
    read_options = pv.ReadOptions(block_size=50 * 1024 * 1024) 
    convert_options = pv.ConvertOptions(column_types=column_types)

    with pv.open_csv(args.input, parse_options=parse_options, read_options=read_options, convert_options=convert_options) as reader:
        
        with tqdm(total=args.total, desc="Analyzing", unit=" rows") as pbar:
            for batch in reader:
                # Convert the 50MB C++ block into Python dictionaries
                records = batch.to_pylist()
                
                for row in records:
                    stats["Total_Entries"] += 1
                    # PyArrow parses booleans as strings "true" / "false"
                    is_reviewed = (row.get("reviewed") == "true")
                    
                    if is_reviewed:
                        stats["Reviewed_Count"] += 1
                    else:
                        stats["Unreviewed_Count"] += 1

                    # 1. Coverage Tracking
                    for field, value in row.items():
                        if value and str(value).strip():
                            coverage_counts[field] += 1

                    # 2. Statistics Tracking
                    p_class = row.get("primary_class")
                    if not p_class: p_class = "Non_Enzyme"
                    
                    class_stats = stats["Classes"][p_class]
                    
                    class_stats["Total"] += 1
                    if is_reviewed:
                        class_stats["Reviewed"] += 1
                        
                    class_stats["Confidence"][row.get("confidence", "None")] += 1
                    class_stats["Enzyme_Type"][row.get("enzyme_type", "None")] += 1
                    class_stats["Functional_Category"][row.get("functional_category", "Unclassified")] += 1
                    class_stats["EC_Completeness"][row.get("ec_completeness", "None")] += 1
                    
                    # Tally Boolean Flags
                    if row.get("flag_has_ec") == "TRUE": class_stats["Evidence_Flags"]["Has_EC"] += 1
                    if row.get("flag_has_active_site") == "TRUE": class_stats["Evidence_Flags"]["Has_Active_Site"] += 1
                    if row.get("flag_has_kinetics") == "TRUE": class_stats["Evidence_Flags"]["Has_Kinetics"] += 1
                    if row.get("flag_has_reaction") == "TRUE": class_stats["Evidence_Flags"]["Has_Reaction"] += 1
                    if row.get("flag_has_experimental") == "TRUE": class_stats["Evidence_Flags"]["Has_Experimental"] += 1
                    if row.get("flag_is_pseudo") == "TRUE": class_stats["Evidence_Flags"]["Is_Pseudo"] += 1
                
                # Update progress bar
                pbar.update(batch.num_rows)
                
                # Free memory instantly
                del records
                del batch

    # =========================================================================
    #   FORMAT AND WRITE OUTPUTS
    # =========================================================================
    
    # Write JSON Statistics
    print(f"\n[Write] Saving statistics to {args.out_stats}...")
    with open(args.out_stats, 'w') as f_json:
        clean_stats = {
            "Total_Entries": stats["Total_Entries"],
            "Reviewed_Count": stats["Reviewed_Count"],
            "Unreviewed_Count": stats["Unreviewed_Count"],
            "Classes": {}
        }
        for cls_name, cls_data in stats["Classes"].items():
            clean_stats["Classes"][cls_name] = {
                "Total": cls_data["Total"],
                "Reviewed": cls_data["Reviewed"],
                "Confidence": dict(cls_data["Confidence"]),
                "Enzyme_Type": dict(cls_data["Enzyme_Type"]),
                "Functional_Category": dict(cls_data["Functional_Category"]),
                "EC_Completeness": dict(cls_data["EC_Completeness"]),
                "Evidence_Flags": dict(cls_data["Evidence_Flags"])
            }
        json.dump(clean_stats, f_json, indent=4)

    # Write TSV Coverage
    print(f"[Write] Saving coverage report to {args.out_cov}...")
    with open(args.out_cov, 'w') as f_cov:
        total_rows = stats["Total_Entries"]
        f_cov.write("Field\tCount_Found\tTotal_Entries\tPercentage_Found\tMissing_Count\n")
        
        for field in headers:
            count = coverage_counts.get(field, 0)
            pct = (count / total_rows * 100) if total_rows > 0 else 0
            missing = total_rows - count
            f_cov.write(f"{field}\t{count}\t{total_rows}\t{pct:.2f}%\t{missing}\n")

    print(f"[Done] Analysis completed in {(time.time() - start_time) / 60:.2f} minutes.")

if __name__ == "__main__":
    main()
