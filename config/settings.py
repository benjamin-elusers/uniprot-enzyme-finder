# config/settings.py

import multiprocessing
import os

# VERSIONING
UNIPROT_RELEASE = "2026_01"
FINDER_VERSION = "v15"
FULL_VERSION = f"{UNIPROT_RELEASE}_{FINDER_VERSION}" # Result: 2026_01_v15

# DIRECTORY MANAGEMENT
# Establish a main output folder and a versioned subfolder
BASE_OUTPUT_DIR = "output"
VERSION_DIR = os.path.join(BASE_OUTPUT_DIR, FULL_VERSION)
os.makedirs(VERSION_DIR, exist_ok=True)


# FILE PATHS
INPUT_FILES = ["uniprot_sprot.xml.gz", "uniprot_trembl.xml.gz"]
# output/[FULL_VERSION]/
RAW_DATA_FILE = os.path.join(VERSION_DIR, "uniprot_raw_features.jsonl.gz")
STATS_FILE = os.path.join(VERSION_DIR, "uniprot_stats_summary.json")
FINAL_OUTPUT_FILE = os.path.join(VERSION_DIR, "uniprot_enzymes.parquet")

# CHECKPOINTING
CHECKPOINT_FILE = os.path.join(VERSION_DIR, f"uniprot_checkpoint_state.json")
CHECKPOINT_FREQ = 1000 # Save state every 1000 chunks

# RESOURCES
OBO_FILE = "data/go-basic.obo"

# PARSING
BATCH_SIZE = 5000
NUM_WORKERS = max(1, multiprocessing.cpu_count() - 2)


