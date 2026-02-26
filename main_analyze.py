# main_analyze.py
import argparse
import json
import time
import os
from collections import defaultdict, Counter
from tqdm import tqdm

import pyarrow as pa
import pyarrow.dataset as ds

from utils.logger import setup_logger
from config.settings import VERSION_DIR, FINAL_OUTPUT_FILE, STATS_FILE

def main():
    # Set default dataset directory based on your unified settings
    ENZ_DATASET_DIR = FINAL_OUTPUT_FILE.replace('.parquet', '_dataset')
    
    parser = argparse.ArgumentParser(description="Stage 3: Generate Coverage and Stats Reports")
    parser.add_argument("--input", "-i", default=ENZ_DATASET_DIR, help="Input Parquet dataset directory")
    parser.add_argument("--out_stats", default=STATS_FILE, help="Output JSON stats file")
    parser.add_argument("--out_cov", default=os.path.join(VERSION_DIR, "uniprot_field_coverage.tsv"), help="Output TSV coverage file")
    args = parser.parse_args()

    # Setup Logging
    log_file = os.path.join(VERSION_DIR, "analysis.log")
    logger = setup_logger("pipeline_logger", log_file)
    logger.info(f"Starting analysis on Parquet Dataset: {args.input}")

    if not os.path.exists(args.input):
        logger.error(f"Dataset directory not found: {args.input}")
        return
    
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

    start_time = time.time()

    # =========================================================================
    #   PYARROW DATASET STREAMING (Ultra-low Memory)
    # =========================================================================
    # Load the entire directory as a single logical dataset
    dataset = ds.dataset(args.input, format="parquet")
    headers = dataset.schema.names
    
    # Parquet stores metadata, so counting total rows across all files is instant
    total_rows = dataset.count_rows()
    logger.info(f"Dataset loaded. Total entries to process: {total_rows:,}")

    try:
        # Stream the dataset in manageable batches
        with tqdm(total=total_rows, desc="Analyzing", unit=" rows") as pbar:
            for batch in dataset.to_batches():
                records = batch.to_pylist()
                
                for row in records:
                    stats["Total_Entries"] += 1
                    
                    # PyArrow parses our string booleans
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
                    if not p_class: 
                        p_class = "Non_Enzyme"
                    
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
                
                # Push the progress bar forward by exactly the chunk size
                pbar.update(batch.num_rows)
                
                # Force memory cleanup for the batch
                del records
                del batch

    except Exception as e:
        logger.error(f"Error during dataset streaming: {e}")
        raise

    # =========================================================================
    #   FORMAT AND WRITE OUTPUTS
    # =========================================================================
    
    # Write JSON Statistics
    logger.info(f"Saving statistics to {args.out_stats}...")
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
    logger.info(f"Saving coverage report to {args.out_cov}...")
    with open(args.out_cov, 'w') as f_cov:
        total_processed = stats["Total_Entries"]
        f_cov.write("Field\tCount_Found\tTotal_Entries\tPercentage_Found\tMissing_Count\n")
        
        for field in headers:
            count = coverage_counts.get(field, 0)
            pct = (count / total_processed * 100) if total_processed > 0 else 0
            missing = total_processed - count
            f_cov.write(f"{field}\t{count}\t{total_processed}\t{pct:.2f}%\t{missing}\n")

    elapsed = (time.time() - start_time) / 60
    logger.info(f"Analysis completed successfully in {elapsed:.2f} minutes.")

if __name__ == "__main__":
    main()
