# main_extract.py

import argparse
import multiprocessing
import os
import time
from tqdm import tqdm

from core.parser import parse_worker
from utils.io import fast_xml_stream, writer_listener, load_extraction_checkpoint, save_extraction_checkpoint, get_estimated_total
from utils.logger import setup_logger
from config.settings import NUM_WORKERS, BATCH_SIZE, RAW_DATA_FILE, CHECKPOINT_FILE, VERSION_DIR, INPUT_FILES

def main():
    parser = argparse.ArgumentParser(description="Stage 1: Extract Raw Features")
    parser.add_argument("--input", "-i", nargs='+', default=INPUT_FILES, help="Input XML.gz files")
    args = parser.parse_args()

    # Setup Logging
    log_file = os.path.join(VERSION_DIR, "extraction.log")
    logger = setup_logger("pipeline_logger", log_file)
    logger.info("Starting UniProt XML extraction pipeline...")

    # Load Checkpoint
    checkpoint = load_extraction_checkpoint(CHECKPOINT_FILE)


    # Setup Writer Queue
    manager = multiprocessing.Manager()
    queue = manager.Queue(maxsize=50)
    pool = multiprocessing.Pool(NUM_WORKERS)
    
    # Start Writer Process
    write_mode = 'at' if checkpoint['entries_done'] > 0 else 'wt'
    writer = multiprocessing.Process(target=writer_listener, args=(queue, RAW_DATA_FILE, write_mode))
    writer.start()

    start_time = time.time()
    total_extracted = 0

    logger.info(f"Extracting features with {NUM_WORKERS} workers...")
    
    for input_file in args.input:
        if not os.path.exists(input_file):
            logger.warning(f"[Skip] File not found: {input_file}")
            continue

        fname = os.path.basename(input_file)
        estimated_total = get_estimated_total(fname)

        skip_count = checkpoint['entries_done'] if checkpoint['last_file'] == fname else 0

        # Create stream generator
        entries_processed = 0
        stream = fast_xml_stream(input_file, BATCH_SIZE)
        
        # Map workers to the stream
        with tqdm(total=estimated_total, desc=f"Parsing {fname}", unit=" proteins") as pbar:
            for batch_results in pool.imap(parse_worker, stream, chunksize=1):
                batch_len = len(batch_results)
                entries_processed += batch_len
                pbar.update(batch_len)
                
                if entries_processed <= skip_count: 
                    continue
                    
                if batch_results:
                    total_extracted += batch_len
                    queue.put(batch_results)
                    # Save progress periodically
                    save_extraction_checkpoint(CHECKPOINT_FILE, fname, entries_processed)
    # Teardown
    queue.put(None)
    writer.join()
    pool.close()
    pool.join()

    elapsed = (time.time() - start_time) / 60
    logger.info(f"Extracted {total_extracted} entries in {elapsed:.1f} minutes.")
    logger.info(f"Output saved to: {RAW_DATA_FILE}")

if __name__ == "__main__":
    main()