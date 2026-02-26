# utils/download_xml.py
import os
import requests
from tqdm import tqdm
import sys

# Add root directory to sys.path so we can import config
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.settings import UNIPROT_RELEASE

FTP_BASE = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"

def download_file(url, dest_path):
    """Downloads a file with a progress bar, skipping if it already exists."""
    if os.path.exists(dest_path):
        print(f"[Skip] File already exists: {dest_path}")
        return

    print(f"[Download] Fetching {url}\n  -> {dest_path}")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(dest_path, 'wb') as f, tqdm(
            desc=os.path.basename(dest_path),
            total=total_size,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for data in response.iter_content(chunk_size=1024*1024):
                size = f.write(data)
                bar.update(size)
    except Exception as e:
        print(f"[Error] Failed to download {url}: {e}")

def main():
    os.makedirs("data", exist_ok=True)
    
    # Define exact URLs
    sprot_url = f"{FTP_BASE}/uniprot_sprot.xml.gz"
    trembl_url = f"{FTP_BASE}/uniprot_trembl.xml.gz"
    relnotes_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt"
    
    # Define expected destination paths
    sprot_dest = f"data/uniprot_sprot_{UNIPROT_RELEASE}.xml.gz"
    trembl_dest = f"data/uniprot_trembl_{UNIPROT_RELEASE}.xml.gz"
    relnotes_dest = f"data/relnotes_{UNIPROT_RELEASE}.txt"
    
    print(f"--- Uniprot Downloader for Version {UNIPROT_RELEASE} ---")
    
    # 1. Download Release Notes
    download_file(relnotes_url, relnotes_dest)
    
    # 2. Download Swiss-Prot
    download_file(sprot_url, sprot_dest)
    
    # 3. Download TrEMBL
    download_file(trembl_url, trembl_dest)
    
    print("\n[Done] Data is ready for extraction.")

if __name__ == "__main__":
    main()