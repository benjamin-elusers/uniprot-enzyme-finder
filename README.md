# Uniprot Enzyme Scanner ![Version](https://img.shields.io/badge/version-v15-blue)

A robust, memory-safe, three-stage pipeline to extract, classify, and analyze enzymes from massive UniProt XML dumps (Swiss-Prot & TrEMBL). 

This pipeline features automatic state-checkpointing, native PyArrow Parquet dataset generation, and centralized logging to handle datasets well over 150GB with minimal RAM footprint.

## Project Structure
* **Stage 1 (Extract):** Streams massive XML files using multiprocessing, extracts raw biological features, and outputs a compressed JSONL file. Includes entry-level checkpointing to resume mid-file if interrupted.
* **Stage 2 (Classify):** Reads the JSONL stream, applies combinatorial logic to determine Enzyme Type and Confidence, and writes native Parquet datasets in chunks. 
* **Stage 3 (Analyze):** Utilizes `pyarrow.dataset` to instantly scan the generated Parquet directories, producing comprehensive TSV coverage reports and detailed combinatorial JSON statistics.

## Installation
1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

2. Optional (Recommended): 

Ensure you have pigz installed on your system for significantly faster parallel XML decompression:
```bash
sudo apt-get install pigz # Ubuntu/Debian
```

```bash
brew install pigz # macOS
```

## Usage

### Step 0: Download Data

Before running the pipeline, you need the UniProt XML datasets.

##### *Option A: Python Script (Recommended)*
This utility automatically detects the latest release, checks file sizes, and handles resumes.

```Bash
python -m utils.download_xml
```
Files will be downloaded to the data/ directory.

##### *Option B: Manual Download*
You can manually download the files from the [UniProt FTP Current Release](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/).

* [uniprot_sprot.xml.gz](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz)
* [uniprot_trembl.xml.gz](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.xml.gz)
* [Release Notes (relnotes.txt)](https://www.google.com/search?q=https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt)

### Step 1: Extract Features (Run once per data update)


This streams the raw XML files and saves the relevant data to output/[VERSION]/uniprot_raw_features.jsonl.gz.
Time estimate: ~2 hours for full Swiss-Prot + TrEMBL.

```bash
python main_extract.py -i data/uniprot_sprot.xml.gz data/uniprot_trembl.xml.gz
```

### Step 2: Classify Enzymes (Run frequently)

This applies the classification logic (Active Sites, Pseudo-enzymes, Confidence scores) and writes natively to a Parquet dataset directory.
Time estimate: ~5-10 minutes.

```bash
python main_classify_count.py -i output/2026_01_v15/uniprot_raw_features.jsonl.gz
```

*If you also want to save the massive dataset of Non-Enzymes, simply add the flag:*
> 
```bash
python main_classify_count.py -i output/2026_01_v15/uniprot_raw_features.jsonl.gz
```

### Step 3: Analyze and Generate Reports
This utilizes PyArrow to instantly read the Parquet datasets, generating the final field coverage report and combinatorial statistics.


```bash
python main_analyze.py -i output/2026_01_v15/uniprot_enzymes_dataset
```

### Configuration
* **Process optimization:**
   * Edit config/settings.py to change the release version, adjust RAM usage (Batch Sizes), or scale worker threads.
* **Biological definition:**
   * Edit config/keywords.py and config/patterns.py to tweak enzyme keyword mappings and regex text mining.
* **Classification Weights:** 
   * Edit core/classifier.py to fine-tune the decision tree for Enzyme/Pseudo-Enzyme detection.


## Project organization
Your project folder should look exactly like this:

```text
uniprot-enzyme-finder/
├── config/
│   ├── keywords.py              # Keyword mappings (CATALYTIC_KW_MAP, FUNCTIONAL_CATEGORIES)
│   ├── patterns.py              # Regex patterns (HARD_INACTIVE, BIOSYNTHESIS)
│   └── settings.py              # Centralized paths, versioning, and tuning (BATCH_SIZE)
├── core/
│   ├── classifier.py            # Main decision tree logic
│   ├── models.py                # Schema definitions and Parquet record formatting
│   └── parser.py                # Fast lxml-based UniProt XML parsing
├── utils/
│   ├── convert_to_parquet.py    # Legacy/Standalone Parquet converter
│   ├── download_xml.py          # Auto-downloader for UniProt releases
│   ├── extract_uniprot_entry.py # Helper to extract single accessions for debugging
│   ├── io.py                    # Streaming (pigz), checkpoints, and queue management
│   ├── logger.py                # Centralized dual-output logging 
│   ├── obo_loader.py            # Gene Ontology DAG loader
│   └── stats.py                 # Combinatorial statistics aggregator
│
├── data/*                       # Contain the XML files and the release notes and other resources
├── output/                      # Output folder
│   └── VERSION/*                # Subfolder for each version of the uniprot release
│
├── main_extract.py              # Stage 1: Memory-safe XML Extraction
├── main_classify_count.py       # Stage 2: Classification and Parquet dataset generation
├── main_analyze.py              # Stage 3: Fast Parquet scanning and reporting
├── requirements.txt             # Python dependencies
└── README.md                    # Project documentation


```
> (Note: output/ and data/ directories are generated automatically during runtime and are ignored by git).

You are now ready to run. 

Start Stage 1 (main_extract.py), and while that churns through the XML, you can open core/classifier.py and start fine-tuning your logic for Stage 2 and Stage 3
