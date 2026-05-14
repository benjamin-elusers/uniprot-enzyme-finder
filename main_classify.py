# main_classify.py

import os
import argparse
import json
import multiprocessing
import time
import itertools
import logging
from tqdm import tqdm

import pyarrow as pa
import pyarrow.parquet as pq

from core.classifier import EnzymeClassifier
from core.models import COLUMNS, BOOL_COLS, INT_COLS, build_parquet_record, 
from utils.io import load_checkpoint, save_checkpoint, chunk_reader, get_estimated_total
from utils.obo_loader import load_catalytic_go_ids, load_eco_definitions
from utils.logger import setup_logger
from config.settings import CHECKPOINT_FILE, FINAL_OUTPUT_FILE, VERSION_DIR, GO_OBO_FILE, ECO_OBO_FILE, RAW_DATA_DIR, NUM_WORKERS, BATCH_SIZE

# =============================================================================
# WORKER CONTEXT & EXECUTION
# =============================================================================

worker_context = {}

def worker_init(valid_go_ids, eco_ids, save_non_enzymes):
    """Initializes the classifier ONCE per worker process."""
    worker_context['classifier'] = EnzymeClassifier(valid_catalytic_go_ids=valid_go_ids, eco_defs=eco_ids)
    worker_context['save_non_enzymes'] = save_non_enzymes

def classify_chunk(lines):
    """Worker: Decodes JSON, classifies, appends to lists of dictionaries."""
    classifier = worker_context['classifier']
    save_non_enzymes = worker_context['save_non_enzymes']
    logger = logging.getLogger("pipeline_logger")
    
    enz_records, non_enz_records = [], []

    for line in lines:
        try:
            entry = json.loads(line)
        except json.JSONDecodeError as e: 
            logger.debug(f"JSON decode error: {e}")
            continue

        try:
            features = entry.get('features', {})
            primary_class, row_metadata = classifier.classify_entry(entry, features)
                
            # Determine Context
            functional_cat = classifier.assign_category(
                entry['lists'].get('kw_ids', []),
                entry['lists'].get('ec_numbers', []),
                primary_class,
                features.get('has_sec_met_text', False)
            )
            
            ec_fam, ec_cls, ec_comp = classifier.get_ec_hierarchy(entry['lists'].get('ec_numbers', []))

            # Build row using the renamed row_metadata dict
            record = build_parquet_record(
                entry, features, primary_class, 
                row_metadata.get('confidence_tier', 'Unrated'), 
                functional_cat, ec_fam, ec_cls, ec_comp, row_metadata
            )

            # Route to appropriate list
            if primary_class in ["Enzyme", "Inferred_Enzyme", "Suspected_Enzyme", "Pseudo_Enzyme"]:
                enz_records.append(record)
            elif save_non_enzymes:
                non_enz_records.append(record)

        except Exception as e:
            logger.error(f"Failed to process entry: {e}")
            
    return enz_records, non_enz_records

# =============================================================================
# MODULAR PIPELINE HELPER FUNCTIONS
# =============================================================================

def parse_arguments():
    """Parses command line arguments."""
    ENZ_DATASET_DIR = FINAL_OUTPUT_FILE.replace('.parquet', '_dataset')
    NON_ENZ_DATASET_DIR = os.path.join(VERSION_DIR, "uniprot_non_enzymes_dataset")

    parser = argparse.ArgumentParser(description="Stage 2: Classify Enzymes (Parquet Export)")
    parser.add_argument("--input", "-i", default=RAW_DATA_DIR, help="Input raw JSONL dataset directory")
    parser.add_argument("--out_enzymes_dir", "-o", default=ENZ_DATASET_DIR, help="Output Parquet directory for enzymes")
    parser.add_argument("--out_non_enzymes_dir", "-n", default=NON_ENZ_DATASET_DIR, help="Output Parquet directory for non-enzymes")
    parser.add_argument("--save_non_enzymes", action="store_true", help="Save non-enzymes to disk")
    
    return parser.parse_args()

def setup_directories(args):
    """Ensures output directories exist."""
    os.makedirs(args.out_enzymes_dir, exist_ok=True)
    if args.save_non_enzymes:
        os.makedirs(args.out_non_enzymes_dir, exist_ok=True)

def get_or_create_checkpoint(input_dir):
    """Loads checkpoint and initializes tracking for the input directory if missing."""
    checkpoint = load_checkpoint(CHECKPOINT_FILE)
    fname = os.path.basename(os.path.normpath(input_dir))
    
    if fname not in checkpoint:
        estimated_total = get_estimated_total() # Global total ~200M
        checkpoint[fname] = {
            "entries_done": 0,
            "total_entries": estimated_total,
            "remaining_entries": estimated_total,
            "completed": False
        }
    return checkpoint, fname

def manage_parquet_writers(action, schema, start_line, args, old_writer_enz=None, old_writer_non=None):
    """Helper to cleanly open, close, and rollover Parquet writers."""
    if action == "close":
        if old_writer_enz: old_writer_enz.close()
        if old_writer_non: old_writer_non.close()
        return None, None
        
    elif action == "open" or action == "rollover":
        if action == "rollover":
            if old_writer_enz: old_writer_enz.close()
            if old_writer_non: old_writer_non.close()

        enz_file = os.path.join(args.out_enzymes_dir, f"part_{start_line}.parquet")
        writer_enz = pq.ParquetWriter(enz_file, schema, compression='snappy')
        
        writer_non = None
        if args.save_non_enzymes:
            non_enz_file = os.path.join(args.out_non_enzymes_dir, f"part_{start_line}.parquet")
            writer_non = pq.ParquetWriter(non_enz_file, schema, compression='snappy')
            
        return writer_enz, writer_non

def process_classification_stream(args, pool, schema, checkpoint, fname):
    """The core execution loop for chunking, mapping, and writing."""
    file_info = checkpoint[fname]
    start_line = file_info["entries_done"]
    
    stream = chunk_reader(args.input, chunk_size=BATCH_SIZE, skip_to=start_line)
    writer_enz, writer_non = manage_parquet_writers("open", schema, start_line, args)
    
    MAX_ENTRIES_PER_FILE = 20_000_000  # Split files every 20M entries
    current_file_start = start_line
    total_enz_written = 0

    try:
        with tqdm(total=file_info["total_entries"], desc="Classifying", unit=" proteins") as pbar:
            pbar.update(start_line)             
            
            while True:
                # 50 chunks * BATCH_SIZE (e.g. 2000) = 100,000 entries per loop
                super_batch = list(itertools.islice(stream, 50))
                if not super_batch: break

                entries_in_batch = sum(len(chunk) for chunk in super_batch)
                results = pool.map(classify_chunk, super_batch)
                
                # Write Results
                for e_records, ne_records in results:
                    if e_records:
                        total_enz_written += len(e_records)
                        writer_enz.write_table(pa.Table.from_pylist(e_records, schema=schema))
                    
                    if args.save_non_enzymes and ne_records:
                        writer_non.write_table(pa.Table.from_pylist(ne_records, schema=schema))

                # Update State
                file_info["entries_done"] += entries_in_batch
                file_info["remaining_entries"] = max(0, file_info["total_entries"] - file_info["entries_done"])
                pbar.update(entries_in_batch)

                # Check Rollover Limit
                if (file_info["entries_done"] - current_file_start) >= MAX_ENTRIES_PER_FILE:
                    writer_enz, writer_non = manage_parquet_writers(
                        "rollover", schema, file_info['entries_done'], args, writer_enz, writer_non
                    )
                    current_file_start = file_info["entries_done"]

                # Checkpoint Save
                if pbar.n % (NUM_WORKERS * 10 * BATCH_SIZE) == 0: 
                    save_checkpoint(CHECKPOINT_FILE, checkpoint)

        # Finalize success
        logging.getLogger("pipeline_logger").info(f"[Done] Finished classifying {fname}.")
        file_info["completed"] = True
        file_info["remaining_entries"] = 0
        file_info["total_entries"] = file_info["entries_done"] 

    except KeyboardInterrupt:
        logging.getLogger("pipeline_logger").warning("Interrupted by user. Saving Checkpoint.")
        
    finally:
        manage_parquet_writers("close", schema, None, args, writer_enz, writer_non)
        save_checkpoint(CHECKPOINT_FILE, checkpoint)
        
    return total_enz_written

# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main():
    start_time = time.time()
    args = parse_arguments()
    setup_directories(args)

    logger = setup_logger("pipeline_logger", os.path.join(VERSION_DIR, "classification.log"))
    logger.info("Starting UniProt Classification Pipeline...")

    # 1. State Check
    checkpoint, fname = get_or_create_checkpoint(args.input)
    if checkpoint[fname].get("completed", False):
        logger.info(f"[Skip] {fname} is already fully classified.")
        return

    # 2. Resource Initialization
    valid_go_ids = load_catalytic_go_ids(GO_OBO_FILE)


    # 2. Build the correct schema dynamically
    schema_fields = []
    for col in COLUMNS:
        if col in BOOL_COLS:
            schema_fields.append((col, pa.bool_()))
        elif col in INT_COLS:
            schema_fields.append((col, pa.int32()))
        else:
            schema_fields.append((col, pa.string()))
    schema = pa.schema(schema_fields)
    #schema = pa.schema([(col, pa.string()) for col in COLUMNS])
    eco_ids = load_eco_definitions(ECO_OBO_FILE)

    pool = multiprocessing.Pool(
        processes=NUM_WORKERS,
        initializer=worker_init, 
        initargs=(valid_go_ids, eco_ids, args.save_non_enzymes),
        maxtasksperchild=100
    )
    
    logger.info(f"Classifying from {args.input} with {NUM_WORKERS} workers...")

    # 3. Execution
    total_enz_written = process_classification_stream(args, pool, schema, checkpoint, fname)

    # 4. Teardown
    pool.close()
    pool.join()

    elapsed = (time.time() - start_time) / 60
    logger.info(f"Done! Time elapsed: {elapsed:.2f} minutes")
    logger.info(f"Wrote {total_enz_written} Enzymes to {args.out_enzymes_dir}")

if __name__ == "__main__":
    main()