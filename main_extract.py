# main_extract.py

import argparse
import multiprocessing
import json
import os
import time
from tqdm import tqdm

from core.parser import UniprotEntryParser
from utils.io import fast_xml_stream, smart_open

# Configuration
BATCH_SIZE = 5000
NUM_WORKERS = max(1, multiprocessing.cpu_count() - 1)

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

def writer_listener(queue, output_file):
    """
    Writes dictionaries to a compressed JSONL file.
    """
    with smart_open(output_file, 'wt') as f_out:
        while True:
            batch = queue.get()
            if batch is None: 
                break
            for entry in batch:
                f_out.write(json.dumps(entry) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Stage 1: Extract Raw Features")
    parser.add_argument("--input", "-i", nargs='+', required=True, help="Input XML.gz files")
    parser.add_argument("--output", "-o", default="uniprot_raw_features.jsonl.gz", help="Output JSONL.gz")
    args = parser.parse_args()

    # Setup Writer Queue
    manager = multiprocessing.Manager()
    queue = manager.Queue(maxsize=50)
    pool = multiprocessing.Pool(NUM_WORKERS)
    
    # Start Writer Process
    writer = multiprocessing.Process(target=writer_listener, args=(queue, args.output))
    writer.start()

    start_time = time.time()
    total_extracted = 0

    print(f"[Start] Extracting features with {NUM_WORKERS} workers...")

    try:
        for input_file in args.input:
            if not os.path.exists(input_file):
                print(f"[Skip] File not found: {input_file}")
                continue

            # Create stream generator
            stream = fast_xml_stream(input_file, BATCH_SIZE)
            
            # Map workers to the stream
            # chunksize=1 ensures the pool doesn't buffer too much RAM
            for batch_results in tqdm(pool.imap(parse_worker, stream, chunksize=1), desc=f"Parsing {os.path.basename(input_file)}"):
                if batch_results:
                    total_extracted += len(batch_results)
                    queue.put(batch_results)

    except KeyboardInterrupt:
        print("\n[Stop] Interrupted by user.")
    finally:
        queue.put(None)
        writer.join()
        pool.close()
        pool.join()

    print(f"\n[Done] Extracted {total_extracted} entries in {(time.time() - start_time) / 60:.1f} minutes.")
    print(f"[Output] {args.output}")

if __name__ == "__main__":
    main()