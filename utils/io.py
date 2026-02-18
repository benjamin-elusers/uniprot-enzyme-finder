# utils/io.py
import subprocess
import shutil
import os
import gzip

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