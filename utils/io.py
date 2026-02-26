# utils/io.py
import subprocess
import shutil
import glob
import os
import json
import gzip
import re

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
    if os.path.exists(checkpoint_file):
        with open(checkpoint_file, 'r') as f:
            return json.load(f).get("processed_lines", 0)
    return 0

def save_checkpoint(checkpoint_file, line_count, version):
    with open(checkpoint_file, 'w') as f:
        json.dump({"processed_lines": line_count, "version": version}, f)

def chunk_reader(filepath, chunk_size=2000, skip_to=0):
    """Yields lists of text lines, skipping already processed ones."""
    from utils.io import smart_open # Ensure this is imported
    with smart_open(filepath, 'rt') as f:
        if skip_to > 0:
            print(f"[Resume] Skipping first {skip_to} lines...")
            for _ in range(skip_to):
                next(f, None)
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk

def get_estimated_total():
    """Attempts to find release notes to get the exact total entries for accurate ETA."""
    total = 200573000 # defaul fallback 
    note_files = glob.glob("relnotes*.txt") + glob.glob("data/relnotes*.txt")
    if note_files:
        try:
            with open(note_files[0], "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
                sp_match = re.search(r"Swiss-Prot:?\s*([\d,]+)\s+entries", content, re.IGNORECASE)
                tr_match = re.search(r"TrEMBL:?\s*([\d,]+)\s+entries", content, re.IGNORECASE)
                sp_count = int(sp_match.group(1).replace(",", "")) if sp_match else 573000
                tr_count = int(tr_match.group(1).replace(",", "")) if tr_match else 200000000
                total = sp_count + tr_count
        except Exception:
            pass
    return total