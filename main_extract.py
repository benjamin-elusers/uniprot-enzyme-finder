# main_extract.py

import argparse
import multiprocessing
import json
import os
import time
from tqdm import tqdm

from core.parser import UniprotEntryParser
from config.settings import FULL_VERSION
from utils.io import fast_xml_stream, smart_open
from config.settings import NUM_WORKERS, BATCH_SIZE, RAW_DATA_FILE, CHECKPOINT_FILE

def parse_worker(xml_batch):
    """
    Worker function: Takes a list of XML strings, returns a list of dicts.
    """
    results = []
    for xml_str in xml_batch:
        try:
            parser = UniprotEntryParser(xml_str)
            data = parser.parse()
            if data:
                results.append(data)
        except Exception as e:
            # In production, log this to a file instead of printing
            pass
    return results

def writer_listener(queue, output_file, mode='wt'):
    """
    Writes dictionaries to a compressed JSONL file.
    """
    with smart_open(output_file, mode) as f_out:
        while True:
            batch = queue.get()
            if batch is None: 
                break
            for entry in batch:
                f_out.write(json.dumps(entry) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Stage 1: Extract Raw Features")
    parser.add_argument("--input", "-i", nargs='+', required=True, help="Input XML.gz files")
    args = parser.parse_args()

    # Load Checkpoint
    checkpoint = {"last_file": None, "entries_done": 0}
    if os.path.exists(CHECKPOINT_FILE):
        with open(CHECKPOINT_FILE, 'r') as f:
            checkpoint = json.load(f)


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

    print(f"[Start] Extracting features with {NUM_WORKERS} workers...")
    for input_file in args.input:

        if not os.path.exists(input_file):
            print(f"[Skip] File not found: {input_file}")
        continue

        fname = os.path.basename(input_file)
        skip_count = checkpoint['entries_done'] if checkpoint['last_file'] == fname else 0

        # Create stream generator
        entries_processed = 0
        stream = fast_xml_stream(input_file, BATCH_SIZE)
        
        # Map workers to the stream
        # chunksize=1 ensures the pool doesn't buffer too much RAM
        for batch_results in tqdm(pool.imap(parse_worker, stream, chunksize=1), desc=f"Parsing {fname}"):
            entries_processed += len(batch_results)
            if entries_processed <= skip_count: continue
            if batch_results:
                total_extracted += len(batch_results)
                queue.put(batch_results)
                # Save progress periodically
                with open(CHECKPOINT_FILE, 'w') as f:  
                    json.dump({"last_file": fname, "entries_done": entries_processed}, f)

    queue.put(None)
    writer.join()
    pool.close()
    pool.join()

    print(f"\n[Done] Extracted {total_extracted} entries in {(time.time() - start_time) / 60:.1f} minutes.")
    print(f"[Output] {RAW_DATA_FILE}")

if __name__ == "__main__":
    main()
