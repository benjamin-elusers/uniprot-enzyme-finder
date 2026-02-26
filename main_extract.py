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

    # Load Checkpoint Dictionary
    checkpoint = load_extraction_checkpoint(CHECKPOINT_FILE)

    # Setup Writer Queue
    manager = multiprocessing.Manager()
    queue = manager.Queue(maxsize=50)
    pool = multiprocessing.Pool(NUM_WORKERS)
    
    # Use append mode if ANY file in the dictionary has processed entries
    write_mode = 'at' if any(info.get('entries_done', 0) > 0 for info in checkpoint.values()) else 'wt'
    writer = multiprocessing.Process(target=writer_listener, args=(queue, RAW_DATA_FILE, write_mode))
    writer.start()

    start_time = time.time()
    total_extracted_this_run = 0

    logger.info(f"Extracting features with {NUM_WORKERS} workers...")
    
    for input_file in args.input:
        if not os.path.exists(input_file):
            logger.warning(f"[Skip] File not found: {input_file}")
            continue

        fname = os.path.basename(input_file)
        
        # Initialize file tracking if it doesn't exist
        if fname not in checkpoint:
            estimated_total = get_estimated_total(fname)
            checkpoint[fname] = {
                "entries_done": 0,
                "total_entries": estimated_total,
                "remaining_entries": estimated_total,
                "completed": False
            }
            
        file_info = checkpoint[fname]

        if file_info.get("completed", False):
            logger.info(f"[Skip] {fname} is already fully processed.")
            continue

        skip_count = file_info["entries_done"]
        stream_entries = 0 # Tracks our position in the current XML stream
        stream = fast_xml_stream(input_file, BATCH_SIZE)
        
        with tqdm(total=file_info["total_entries"], desc=f"Parsing {fname}", unit=" proteins") as pbar:
            pbar.update(skip_count) # Advance the visual progress bar to where we left off
            
            for batch_results in pool.imap(parse_worker, stream, chunksize=1):
                batch_len = len(batch_results)
                stream_entries += batch_len
                
                # Fast-forward until we reach the skipped point
                if stream_entries <= skip_count: 
                    continue
                    
                if batch_results:
                    total_extracted_this_run += batch_len
                    queue.put(batch_results)
                    
                    # Update dictionary state
                    file_info["entries_done"] += batch_len
                    file_info["remaining_entries"] = max(0, file_info["total_entries"] - file_info["entries_done"])
                    save_extraction_checkpoint(CHECKPOINT_FILE, checkpoint)
                    
                    pbar.update(batch_len)

        # Mark file as complete once the loop finishes successfully
        logger.info(f"[Done] Finished processing {fname}.")
        file_info["completed"] = True
        file_info["remaining_entries"] = 0
        # Sync the total to the exact amount processed just in case the estimate was off
        file_info["total_entries"] = file_info["entries_done"] 
        save_extraction_checkpoint(CHECKPOINT_FILE, checkpoint)

    # Teardown
    queue.put(None)
    writer.join()
    pool.close()
    pool.join()

    elapsed = (time.time() - start_time) / 60
    logger.info(f"Extracted {total_extracted_this_run} entries this session in {elapsed:.1f} minutes.")
    logger.info(f"Output appended to: {RAW_DATA_FILE}")

if __name__ == "__main__":
    main()