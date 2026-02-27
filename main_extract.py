# main_extract.py
import argparse
import multiprocessing
import os
import time
import json
import gzip
import logging
from tqdm import tqdm

from core.parser import UniprotEntryParser
from utils.io import fast_xml_stream, load_checkpoint, save_checkpoint, get_estimated_total
from utils.logger import setup_logger
from config.settings import NUM_WORKERS, BATCH_SIZE, RAW_DATA_DIR, CHECKPOINT_FILE, VERSION_DIR, INPUT_FILES

def parse_and_write_worker(xml_batch):
    """
    Worker parses XML and appends directly to its own designated file.
    """
    # Get the unique name of the current process (e.g., 'ForkPoolWorker-1')
    worker_name = multiprocessing.current_process().name
    out_path = os.path.join(RAW_DATA_DIR, f"{worker_name}.jsonl.gz")
    
    logger = logging.getLogger("pipeline_logger")
    results = []
    
    for xml_str in xml_batch:
        try:
            parser = UniprotEntryParser(xml_str)
            data = parser.parse()
            if data:
                results.append(data)
        except Exception as e:
            logger.debug(f"Failed to parse XML: {e}")
            
    if results:
        # Append mode ('at') seamlessly adds a new compressed block to the worker's file
        with gzip.open(out_path, 'at', encoding='utf-8') as f:
            for res in results:
                f.write(json.dumps(res) + "\n")
                
    # Return only the integer count back to the main thread (Memory safe!)
    return len(results)

def main():
    parser = argparse.ArgumentParser(description="Stage 1: Extract Raw Features")
    parser.add_argument("--input", "-i", nargs='+', default=INPUT_FILES, help="Input XML.gz files")
    args = parser.parse_args()

    log_file = os.path.join(VERSION_DIR, "extraction.log")
    logger = setup_logger("pipeline_logger", log_file)
    logger.info("Starting UniProt XML extraction pipeline...")

    checkpoint = load_checkpoint(CHECKPOINT_FILE)
    os.makedirs(RAW_DATA_DIR, exist_ok=True)

    pool = multiprocessing.Pool(NUM_WORKERS)
    start_time = time.time()
    total_extracted_this_run = 0

    logger.info(f"Extracting features to {RAW_DATA_DIR} with {NUM_WORKERS} workers...")
    
    for input_file in args.input:
        if not os.path.exists(input_file):
            logger.warning(f"[Skip] File not found: {input_file}")
            continue

        fname = os.path.basename(input_file)
        
        # Initialize tracking
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
        stream = fast_xml_stream(input_file, BATCH_SIZE)
        
        # This generator acts as a gatekeeper, preventing skipped XML strings 
        # from ever being sent to the workers or taking up RAM
        def generate_payloads():
            stream_entries = 0
            for batch in stream:
                stream_entries += len(batch)
                if stream_entries <= skip_count:
                    continue
                yield batch

        with tqdm(total=file_info["total_entries"], desc=f"Parsing {fname}", unit=" proteins") as pbar:
            pbar.update(skip_count) 
            
            # imap_unordered ensures we only load as many batches into memory as we have workers
            for batch_len in pool.imap_unordered(parse_and_write_worker, generate_payloads(), chunksize=1):
                total_extracted_this_run += batch_len
                
                # State updates
                file_info["entries_done"] += batch_len
                file_info["remaining_entries"] = max(0, file_info["total_entries"] - file_info["entries_done"])
                save_checkpoint(CHECKPOINT_FILE, checkpoint)
                
                pbar.update(batch_len)

        logger.info(f"[Done] Finished processing {fname}.")
        file_info["completed"] = True
        file_info["remaining_entries"] = 0
        file_info["total_entries"] = file_info["entries_done"] 
        save_checkpoint(CHECKPOINT_FILE, checkpoint)

    pool.close()
    pool.join()

    elapsed = (time.time() - start_time) / 60
    logger.info(f"Extracted {total_extracted_this_run} entries this session in {elapsed:.1f} minutes.")
    logger.info(f"Output saved to Dataset Directory: {RAW_DATA_DIR}")

if __name__ == "__main__":
    main()