# config/settings.py

import multiprocessing

# =============================================================================
#   VERSIONING
# =============================================================================
UNIPROT_RELEASE = "2026_01"
FINDER_VERSION = "v15"
FULL_VERSION = f"{UNIPROT_RELEASE}_{FINDER_VERSION}" # Result: 2026_01_v15

# =============================================================================
#   FILE PATHS
# =============================================================================
# Default Input/Output
INPUT_FILES = ["uniprot_sprot.xml.gz", "uniprot_trembl.xml.gz"]
RAW_DATA_FILE = f"uniprot_raw_features_{FULL_VERSION}.jsonl.gz"
FINAL_OUTPUT_FILE = f"uniprot_enzymes_{FULL_VERSION}.tsv.gz"
CHECKPOINT_FILE = f"uniprot_checkpoint_{FULL_VERSION}.json" #

# Resources
OBO_FILE = "data/go-basic.obo"
STATS_FILE = f"uniprot_stats_summary_{FULL_VERSION}.json"

# =============================================================================
#   TUNING
# =============================================================================
# Parsing
BATCH_SIZE = 5000
NUM_WORKERS = max(1, multiprocessing.cpu_count() - 2)

# Checkpointing
CHECKPOINT_FREQ = 1000 # Save state every 1000 chunks
