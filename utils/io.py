# utils/io.py
import subprocess
import shutil
import glob
import os
import json
import gzip
import re
import logging
from pathlib import Path

from config.settings import UNIPROT_RELEASE

def get_estimated_total(filename=None):
    """
    Finds the specific release notes to get exact total entries.
    If a filename is provided, it returns the specific count for that file.
    Otherwise, it returns the combined total.
    """
    sp_count = 573000
    tr_count = 200000000
    
    target_note = f"data/relnotes_{UNIPROT_RELEASE}.txt"
    fallback_note = f"relnotes_{UNIPROT_RELEASE}.txt"
    
    note_file = None
    if os.path.exists(target_note):
        note_file = target_note
    elif os.path.exists(fallback_note):
        note_file = fallback_note

    if note_file:
        try:
            with open(note_file, "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
                sp_match = re.search(r"Swiss-Prot:?\s*([\d,]+)\s+entries", content, re.IGNORECASE)
                tr_match = re.search(r"TrEMBL:?\s*([\d,]+)\s+entries", content, re.IGNORECASE)
                
                if sp_match: sp_count = int(sp_match.group(1).replace(",", ""))
                if tr_match: tr_count = int(tr_match.group(1).replace(",", ""))
        except Exception:
            pass

    if filename:
        if "sprot" in filename.lower():
            return sp_count
        elif "trembl" in filename.lower():
            return tr_count
            
    return sp_count + tr_count

def smart_open(filepath, mode='rt'):
    """
    Opens a file with gzip if it ends in .gz, otherwise normal open.
    """
    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode, encoding='utf-8')
    return open(filepath, mode, encoding='utf-8')

def fast_xml_stream(filepath, batch_size=5000):
    """
    Yields list of XML strings using pigz for parallel decompression 
    if available, otherwise gunzip.
    """
    # Detect fast decompressor
    if shutil.which("pigz"):
        cmd = ['pigz', '-dc', filepath]
    else:
        cmd = ['gunzip', '-c', filepath]

    print(f"[IO] Streaming {filepath} using {' '.join(cmd)}...")
    
    process = subprocess.Popen(
        cmd, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.DEVNULL, 
        encoding='utf-8', 
        bufsize=4*1024*1024
    )

    batch = []
    current_entry = []
    in_entry = False

    try:
        for line in process.stdout:
            if '<entry' in line:
                in_entry = True
                current_entry = [line]
            elif in_entry:
                current_entry.append(line)
                if '</entry>' in line:
                    in_entry = False
                    batch.append("".join(current_entry))
                    
                    if len(batch) >= batch_size:
                        yield batch
                        batch = []
        
        # Yield remaining
        if batch:
            yield batch

    finally:
        process.terminate()

def load_checkpoint(checkpoint_file):
    """Loads the extraction state dictionary."""
    if os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, 'r') as f:
                return json.load(f)
        except Exception:
            pass
    return {} # Return an empty dict if no checkpoint exists

def save_checkpoint(checkpoint_file, checkpoint_data):
    """Saves the full extraction progress dictionary with indentation for readability."""
    with open(checkpoint_file, 'w') as f:  
        json.dump(checkpoint_data, f, indent=4)

def get_dataset_files(filepath_or_dir, pattern="*.jsonl.gz"):
    """
    Resolves a file or directory path into a sorted list of file paths.
    """
    path_obj = Path(filepath_or_dir)
    
    if path_obj.is_dir():
        # .glob() returns Path objects, we convert them to strings
        return sorted(str(p) for p in path_obj.glob(pattern))
    elif path_obj.exists():
        return [str(path_obj)]
    
    return []

def chunk_reader(filepath_or_dir, chunk_size=2000, skip_to=0):
    """
    Yields lists of text lines, seamlessly handling single files or dataset directories.
    """
    from utils.io import smart_open 
    logger = logging.getLogger("pipeline_logger")
    
    # Delegate the file resolution to the helper function
    files = get_dataset_files(filepath_or_dir)
        
    if not files:
        logger.error(f"[Error] No files found for path: {filepath_or_dir}")
        return
        
    if skip_to > 0:
        logger.info(f"Resuming... Skipping first {skip_to} lines across the dataset.")
        
    lines_to_skip = skip_to
    chunk = []
    
    for fpath in files:
        with smart_open(fpath, 'rt') as f:
            for line in f:
                if lines_to_skip > 0:
                    lines_to_skip -= 1
                    continue
                    
                chunk.append(line)
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
    if chunk:
        yield chunk

def writer_listener(queue, output_file, mode='wt'):
    """Writes dictionaries from a multiprocessing queue to a compressed JSONL file."""
    # Ensure smart_open is imported or available in this file
    with smart_open(output_file, mode) as f_out:
        while True:
            batch = queue.get()
            if batch is None: 
                break
            for entry in batch:
                f_out.write(json.dumps(entry) + "\n")

