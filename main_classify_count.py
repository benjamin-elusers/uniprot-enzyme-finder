# main_classify_count.py

import os
import argparse
import json
import multiprocessing
import time
import itertools
import logging
from collections import Counter
from tqdm import tqdm

import pyarrow as pa
import pyarrow.parquet as pq

from core.classifier import EnzymeClassifier
from core.models import COLUMNS, build_parquet_record
from utils.io import load_checkpoint save_checkpoint chunk_reader, get_estimated_total
from utils.stats import generate_stats_report
from utils.obo_loader import load_catalytic_go_ids
from utils.logger import setup_logger
from config.settings import CHECKPOINT_FILE, FULL_VERSION, FINAL_OUTPUT_FILE, VERSION_DIR, OBO_FILE
from config.keywords import CATALYTIC_KW_MAP

worker_context = {}

def worker_init(valid_go_ids, save_non_enzymes):
    """Initializes the classifier ONCE per worker process."""
    worker_context['classifier'] = EnzymeClassifier(valid_catalytic_go_ids=valid_go_ids)
    worker_context['save_non_enzymes'] = save_non_enzymes

def classify_chunk(lines):
    """Worker: Decodes JSON, classifies, appends to lists of dictionaries."""
    classifier = worker_context['classifier']
    save_non_enzymes = worker_context['save_non_enzymes']
    logger = logging.getLogger("pipeline_logger")
    
    enz_records, non_enz_records = [], []
    stats = Counter()
    enz_count, non_enz_count = 0, 0

    for line in lines:
        try:
            entry = json.loads(line)
        except json.JSONDecodeError as e: 
            logger.debug(f"JSON decode error: {e}")
            continue
            
        feats = entry.get('features', {})
        ec_list = entry.get('lists', {}).get('ec_numbers', [])
        feats['has_catalytic_kw'] = any(k in CATALYTIC_KW_MAP for k in entry.get('lists', {}).get('kw_ids', []))
        
        primary_class, stats_flags = classifier.classify_entry(entry, feats)
        fam, cls, comp = classifier.get_ec_hierarchy(ec_list)
        
        has_rxn = 1 if stats_flags.get("flag_has_reaction") else 0
        has_act = 1 if stats_flags.get("flag_has_active_site") else 0
        has_kin = 1 if stats_flags.get("flag_has_kinetics") else 0
        has_kw = 1 if feats.get("has_catalytic_kw") else 0
        has_exp = 1 if stats_flags.get("flag_has_experimental") else 0
        
        combo_sig = f"{primary_class} | EC:{comp}, Rxn:{has_rxn}, ActSite:{has_act}, Kin:{has_kin}, KW:{has_kw}, Exp:{has_exp}"
        stats[combo_sig] += 1
        stats[f"TOTAL_{primary_class}"] += 1

        if primary_class == "Non_Enzyme" and not save_non_enzymes:
            non_enz_count += 1
            continue

        conf = "None"
        if primary_class == "Enzyme": conf = "High"
        elif primary_class == "Suspected_Enzyme": conf = "Medium"
        elif primary_class == "Pseudo_Enzyme": conf = "Bronze"

        cat = classifier.assign_category(entry.get('lists', {}).get('kw_ids', []), ec_list, primary_class, feats.get('has_text_mining', False))

        row = build_parquet_record(entry, feats, primary_class, conf, cat, fam, cls, comp, stats_flags)
        
        if primary_class == "Non_Enzyme":
            non_enz_records.append(row)
            non_enz_count += 1
        else:
            enz_records.append(row)
            enz_count += 1
            
    return (enz_count, enz_records, non_enz_count, non_enz_records, stats)


def main():
    ENZ_DATASET_DIR = FINAL_OUTPUT_FILE.replace('.parquet', '_dataset')
    NON_ENZ_DATASET_DIR = os.path.join(VERSION_DIR, "uniprot_non_enzymes_dataset")

    parser = argparse.ArgumentParser(description="Stage 2: Classify Enzymes")
    parser.add_argument("--input", "-i", required=True, help="Input raw JSONL.gz")
    parser.add_argument("--out_enzymes_dir", "-o", default=ENZ_DATASET_DIR, help="Output enzymes Parquet directory")
    parser.add_argument("--out_non_enzymes_dir", "-n", default=NON_ENZ_DATASET_DIR, help="Output non-enzymes Parquet directory")
    parser.add_argument("--out_stats", "-s", default=os.path.join(VERSION_DIR, "uniprot_combinatorial_stats.json"), help="Output Stats JSON")
    parser.add_argument("--save_non_enzymes", action="store_true", help="Pass this flag to save the non-enzymes to disk.")
    args = parser.parse_args()

    # Setup Logging
    log_file = os.path.join(VERSION_DIR, "classification.log")
    logger = setup_logger("pipeline_logger", log_file)
    logger.info("Starting UniProt classification pipeline...")

    os.makedirs(args.out_enzymes_dir, exist_ok=True)
    if args.save_non_enzymes:
        os.makedirs(args.out_non_enzymes_dir, exist_ok=True)

    # --- Load Checkpoint Dictionary ---
    checkpoint = load_checkpoint(CHECKPOINT_FILE)
    fname = os.path.basename(args.input)
    
    # Initialize Tracking for this JSONL file if it doesn't exist
    if fname not in checkpoint:
        estimated_total = get_estimated_total() # Combined total from relnotes
        checkpoint[fname] = {
            "entries_done": 0,
            "total_entries": estimated_total,
            "remaining_entries": estimated_total,
            "completed": False
        }
    
    file_info = checkpoint[fname]
    
    if file_info.get("completed", False):
        logger.info(f"[Skip] {fname} is already fully classified.")
        return

    start_line = file_info["entries_done"]

    # --- Worker Initialization ---
    valid_go_ids = load_catalytic_go_ids(OBO_FILE)
    num_workers = max(1, multiprocessing.cpu_count() - 1)
    schema = pa.schema([(col, pa.string()) for col in COLUMNS])
    pool = multiprocessing.Pool(num_workers, initializer=worker_init, initargs=(valid_go_ids, args.save_non_enzymes), maxtasksperchild=100)
    
    logger.info(f"Classifying from {args.input} with {num_workers} workers...")
    
    start_time = time.time()
    total_enz_written, total_non_enz_written = 0, 0
    global_stats = Counter()
    
    stream = chunk_reader(args.input, chunk_size=2000, skip_to=start_line)
    
    # Generate unique parquet file parts so resumes don't overwrite previous data
    enz_file = os.path.join(args.out_enzymes_dir, f"part_{start_line}.parquet")
    non_enz_file = os.path.join(args.out_non_enzymes_dir, f"part_{start_line}.parquet")

    writer_enz = pq.ParquetWriter(enz_file, schema, compression='snappy')
    writer_non = pq.ParquetWriter(non_enz_file, schema, compression='snappy') if args.save_non_enzymes else None

    try:
        with tqdm(total=file_info["total_entries"], desc="Classifying", unit=" proteins") as pbar:
            pbar.update(start_line) 
            
            while True:
                super_batch = list(itertools.islice(stream, num_workers * 4))
                if not super_batch: break
                    
                entries_in_batch = sum(len(chunk) for chunk in super_batch)
                results = pool.map(classify_chunk, super_batch)
                
                for e_cnt, e_records, ne_cnt, ne_records, chunk_stats in results:
                    global_stats.update(chunk_stats)

                    if e_cnt > 0:
                        total_enz_written += e_cnt
                        writer_enz.write_table(pa.Table.from_pylist(e_records, schema=schema))
                    
                    if ne_cnt > 0 and args.save_non_enzymes:
                        total_non_enz_written += ne_cnt
                        writer_non.write_table(pa.Table.from_pylist(ne_records, schema=schema))

                # Update State
                file_info["entries_done"] += entries_in_batch
                file_info["remaining_entries"] = max(0, file_info["total_entries"] - file_info["entries_done"])
                
                # Push the progress bar forward
                pbar.update(entries_in_batch)

                # Periodically save state to disk
                if pbar.n % (num_workers * 10 * 2000) == 0: 
                    save_checkpoint(CHECKPOINT_FILE, checkpoint)

        # Mark file as complete once the loop finishes successfully
        logger.info(f"[Done] Finished classifying {fname}.")
        file_info["completed"] = True
        file_info["remaining_entries"] = 0
        file_info["total_entries"] = file_info["entries_done"] 
        save_checkpoint(CHECKPOINT_FILE, checkpoint)
    except KeyboardInterrupt:
        logger.warning("Interrupted by user. Checkpoint saved.")
    finally:
        writer_enz.close()
        if writer_non: writer_non.close()
        pool.close()
        pool.join()
        save_checkpoint(CHECKPOINT_FILE, checkpoint)

    generate_stats_report(global_stats, args.out_stats)

    elapsed = (time.time() - start_time) / 60
    logger.info(f"Done! Time elapsed: {elapsed:.2f} minutes")
    logger.info(f"Wrote {total_enz_written} Enzymes to {args.out_enzymes_dir}")

if __name__ == "__main__":
    main()