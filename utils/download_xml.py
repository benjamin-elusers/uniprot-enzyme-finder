import os
import re
import requests
import shutil
from pathlib import Path
from tqdm import tqdm

# Configuration
BASE_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release"
DATA_DIR = Path("../data")  # Saves to project_root/data/

def get_remote_file_size(url):
    """Returns the size of the remote file in bytes, or None if failed."""
    try:
        response = requests.head(url, allow_redirects=True)
        return int(response.headers.get('content-length', 0))
    except Exception:
        return None

def download_file(url, dest_path):
    """
    Downloads a file with a progress bar. 
    Skips if local file exists and has the same size as remote.
    """
    dest_path = Path(dest_path)
    remote_size = get_remote_file_size(url)
    
    # 1. Check if file exists and matches size
    if dest_path.exists():
        local_size = dest_path.stat().st_size
        if remote_size and local_size == remote_size:
            print(f"[Skip] {dest_path.name} already exists and matches remote size ({local_size / 1e9:.2f} GB).")
            return
        else:
            print(f"[Update] Local file size ({local_size}) differs from remote ({remote_size}). Re-downloading...")

    # 2. Download with streaming
    print(f"[Download] Fetching {url} -> {dest_path}...")
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            total_size = int(r.headers.get('content-length', 0))
            
            with open(dest_path, 'wb') as f, tqdm(
                desc=dest_path.name,
                total=total_size,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
            ) as bar:
                for chunk in r.iter_content(chunk_size=8192):
                    size = f.write(chunk)
                    bar.update(size)
    except KeyboardInterrupt:
        print("\n[Stop] Download interrupted. Deleting incomplete file.")
        dest_path.unlink(missing_ok=True)
        raise
    except Exception as e:
        print(f"[Error] Failed to download {url}: {e}")

def get_latest_version():
    """Fetches release notes to determine the current UniProt version (e.g., 2024_01)."""
    relnotes_url = f"{BASE_URL}/relnotes.txt"
    print(f"[Init] Fetching version info from {relnotes_url}...")
    
    try:
        response = requests.get(relnotes_url)
        response.raise_for_status()
        content = response.text
        
        # Regex to find: "UniProtKB Release 2024_01"
        match = re.search(r"UniProtKB Release (\d{4}_\d{2})", content)
        if match:
            version = match.group(1)
            print(f"[Info] Detected UniProt Release: {version}")
            return version
        else:
            print("[Warn] Could not regex parse version. Defaulting to 'unknown'.")
            return "unknown"
    except Exception as e:
        print(f"[Error] Failed to fetch release notes: {e}")
        return "unknown"

def main():
    # Ensure data directory exists
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Get Version
    version = get_latest_version()

    # 2. Define Targets
    files_to_download = [
        {
            "url": f"{BASE_URL}/knowledgebase/complete/uniprot_sprot.xml.gz",
            "name": f"uniprot_sprot_{version}.xml.gz"
        },
        {
            "url": f"{BASE_URL}/knowledgebase/complete/uniprot_trembl.xml.gz",
            "name": f"uniprot_trembl_{version}.xml.gz"
        }
    ]

    # 3. Save Release Notes
    with open(DATA_DIR / f"relnotes_{version}.txt", "w") as f:
        f.write(f"Version detected: {version}\nDownloaded on: {os.getcwd()}\n")

    # 4. Download Files
    for item in files_to_download:
        dest = DATA_DIR / item["name"]
        download_file(item["url"], dest)

    print("\n[Done] Data directory is ready.")

if __name__ == "__main__":
    main()