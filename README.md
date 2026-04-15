# Uniprot Enzyme Scanner ![Version](https://img.shields.io/badge/version-v26-blue)

A robust, memory-safe, three-stage pipeline to extract, classify, and analyze enzymes from massive UniProt XML dumps (Swiss-Prot & TrEMBL). 

This pipeline features automatic state-checkpointing, native PyArrow Parquet dataset generation, and centralized logging to handle datasets well over 150GB with minimal RAM footprint.


## Key Features
* **High-Performance XML Parsing:** Uses `lxml.etree` iterative parsing and multiprocessing to gracefully handle massive (>150GB) XML dumps without blowing out memory.
* **Tiered Confidence Scoring:** Radically decouples a protein's predicted *Identity* (e.g., `Inferred_Enzyme`) from its Annotation *Confidence* (`Gold`, `Silver`, `Bronze`, and a quantitative `0-100+` score).
* **Precise Evidence Provenance:** Maps XML `<evidence>` keys to specific ECO ontology codes to accurately distinguish between strict low-throughput manual assays, high-throughput screens, manual curation, and automated electronic inferences.
* **Dynamic GO Clustering:** A dedicated, lightweight Pandas stage allows you to rapidly map enzymes into high-level functional groups based on GO Slim terms without re-running the heavy XML parser.
* **Kinetics Parsing:** Uses linguistic pivots to extract `KM`, `Vmax`, and `kcat` directly from free-text comment strings.


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
Time estimate: ~2-3 hours for full Swiss-Prot + TrEMBL.

```bash
python main_extract.py # use input files defined in settings.py

# For previous releases you can still give the correct path for each XML files
python main_extract.py -i data/uniprot_sprot_2026_01.xml.gz data/uniprot_trembl_2026_01.xml.gz

```

### Step 2: Classify Enzymes (Run frequently)

This applies the classification logic (Active Sites, Pseudo-enzymes, Confidence scores...) and writes natively to a Parquet dataset directory.
Time estimate: ~60-90 minutes.

```bash
python main_classify.py -i output/2026_01_v26/uniprot_raw_features.jsonl.gz
```

*If you also want to save the massive dataset of Non-Enzymes, simply add the flag:*
> 
```bash
python main_classify.py -i output/2026_01_v26/uniprot_raw_features.jsonl.gz --save_non_enzymes
```

### Step 3: Analyze and Generate Reports
This utilizes PyArrow to instantly read the Parquet datasets, generating the final field coverage report and combinatorial statistics.


```bash
python main_analyze.py -i output/2026_01_v26/uniprot_enzymes_dataset
```

### Configuration
* **Process optimization:**
   * Edit config/settings.py to change the release version, adjust RAM usage (Batch Sizes), or scale worker threads.
* **Biological definition:**
   * Edit config/keywords.py and config/patterns.py to tweak enzyme keyword mappings and regex text mining.
* **Classification Weights:** 
   * Edit core/classifier.py to fine-tune the decision tree for Enzyme/Pseudo-Enzyme detection.


## Project organization

The repository follows this modular structure:

```text
uniprot-enzyme-finder/
│
├── analysis/                    # R scripts for downstream plotting and analysis
│   ├── 01_prepare_data.r        # Base data preparation
│   ├── 02_plot_ec_donut.r       # Visualization scripts (Donuts, Waffles, Dashboards)
│   ├── ...                      # (Various plot generation scripts)
│   └── 13_plot_yeast.r          # Species-specific subset analysis
│
├── config/
│   ├── evidences.py             # ECO evidence code hierarchy mapping
│   ├── go_clusters.py           # GO Slim to Functional Cluster mapping
│   ├── keywords.py              # UniProt Keyword hierarchical mapping
│   ├── patterns.py              # Regex patterns for ECs, names, and text-mining
│   └── settings.py              # Global paths, versioning, and multiprocessing config
│
├── core/
│   ├── classifier.py            # Primary classification logic and weighted confidence scoring
│   ├── models.py                # Flat Parquet Schema definitions
│   └── parser.py                # Fast lxml-based UniProt XML parsing
│
├── utils/
│   ├── convert_to_parquet.py    # Legacy/Standalone Parquet converter
│   ├── download_xml.py          # Auto-downloader for UniProt releases
│   ├── extract_uniprot_entry.py # Helper to extract single accessions for debugging
│   ├── io.py                    # Streaming (pigz), checkpoints, and queue management
│   ├── logger.py                # Centralized dual-output logging 
│   ├── obo_loader.py            # Gene Ontology DAG loader
│   ├── parse_kinetics_jsonl.py  # Utility to parse kinetics from raw JSONL output
│   └── stats.py                 # Combinatorial statistics aggregator
│
├── data/*                       # Contain the XML files and the release notes and other resources
├── output/                      # Output folder
│   └── VERSION/*                # Subfolder for each version of the uniprot release
│
├── main_extract.py              # Stage 1: Memory-safe XML Extraction to Raw JSONL
├── main_classify.py             # Stage 2: Fast Parquet Classification (No Aggregation)
├── main_analyze.py              # Stage 4: Parquet scanning, counting, and reporting
├── requirements.txt             # Python dependencies
└── README.md                    # Project documentation

```
> (Note: output/ and data/ directories are generated automatically during runtime and are ignored by git).

You are now ready to run. 

Start Stage 1 (main_extract.py), and while that churns through the XML, you can open core/classifier.py and start fine-tuning your logic for Stage 2 and Stage 3
