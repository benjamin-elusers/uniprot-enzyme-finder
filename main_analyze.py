# main_analyze.py

import argparse
import csv
import json
import time
from collections import defaultdict, Counter
from tqdm import tqdm
from utils.io import smart_open

def main():
    parser = argparse.ArgumentParser(description="Stage 3: Generate Coverage and Stats Reports")
    parser.add_argument("--input", "-i", required=True, help="Input final TSV.gz file")
    parser.add_argument("--out_stats", default="uniprot_stats_summary_v14.json", help="Output JSON stats file")
    parser.add_argument("--out_cov", default="uniprot_field_coverage_v14.tsv", help="Output TSV coverage file")
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

    print(f"[Start] Analyzing {args.input}...")
    start_time = time.time()

    # =========================================================================
    #   SINGLE PASS STREAMING
    # =========================================================================
    
    with smart_open(args.input, 'rt') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames

        for row in tqdm(reader, desc="Generating Reports", unit=" proteins"):
            stats["Total_Entries"] += 1
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
            p_class = row.get("primary_class", "Non_Enzyme")
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

    # =========================================================================
    #   FORMAT AND WRITE OUTPUTS
    # =========================================================================
    
    # Write JSON Statistics
    print(f"[Write] Saving statistics to {args.out_stats}...")
    with open(args.out_stats, 'w') as f_json:
        # Convert defaultdicts and Counters to standard dicts for clean JSON serialization
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
        
        for field in fieldnames:
            count = coverage_counts.get(field, 0)
            pct = (count / total_rows * 100) if total_rows > 0 else 0
            missing = total_rows - count
            f_cov.write(f"{field}\t{count}\t{total_rows}\t{pct:.2f}%\t{missing}\n")

    print(f"[Done] Analysis completed in {(time.time() - start_time) / 60:.2f} minutes.")

if __name__ == "__main__":
    main()