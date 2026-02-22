# config/settings.py
import multiprocessing

# =============================================================================
#   FILE PATHS
# =============================================================================

# Default Input/Output
INPUT_FILES = ["uniprot_sprot.xml.gz", "uniprot_trembl.xml.gz"]
RAW_DATA_FILE = "uniprot_raw_features.jsonl.gz"
FINAL_OUTPUT_FILE = "uniprot_enzymes_v14.tsv.gz"

# Resources
OBO_FILE = "data/go-basic.obo"
STATS_FILE = "uniprot_stats_summary.json"

# =============================================================================
#   TUNING
# =============================================================================

# Parsing
BATCH_SIZE = 5000
NUM_WORKERS = max(1, multiprocessing.cpu_count() - 2)

# Checkpointing
CHECKPOINT_FREQ = 100
CHECKPOINT_FILE = "uniprot_parse_state.json"