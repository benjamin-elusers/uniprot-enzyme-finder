# Uniprot Enzyme Scanner (v14)

A two-stage pipeline to extract and classify enzymes from UniProt XML dumps (Swiss-Prot & TrEMBL).

## Project Structure
* **Stage 1 (Extract):** Parses massive XML files and extracts raw features into a compressed JSONL file. (I/O Bound)
* **Stage 2 (Classify):** Reads the JSONL, applies logic to determine Enzyme Type/Confidence, and outputs a TSV. (CPU Bound)

## Installation
1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

2. Optional (Recommended): 

Ensure you have pigz installed on your system for faster processing:

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

This reads the raw XML files and saves the relevant data to uniprot_raw_features.jsonl.gz.
Time estimate: ~2 hours for full Swiss-Prot + TrEMBL.

```bash
python main_extract.py -i data/uniprot_sprot.xml.gz data/uniprot_trembl.xml.gz
```

### Step 2: Classify Enzymes (Run frequently)

This applies your classification logic (Active Sites, Pseudo-enzymes, Confidence scores).
Time estimate: ~5-10 minutes.

```bash
python main_classify_count.py -i output/2026_01_v15/uniprot_raw_features.jsonl.gz
```

*If you also want to save the massive dataset of Non-Enzymes, simply add the flag:*
> 
```bash
python main_classify_count.py -i output/2026_01_v15/uniprot_raw_features.jsonl.gz
```

### Step 3: Generate reports
Generate the final coverage report and JSON statistics.

```bash
python main_analyze.py -i output/2026_01_v15/uniprot_enzymes_dataset
```

### Configuration:
* Edit config/keywords.py to change enzyme keyword mappings.
* Edit config/settings.py to adjust batch sizes or worker counts.
* Edit core/classifier.py to tweak the logic for "Gold/Silver" confidence or "Pseudo-Enzyme" detection.


## Project organization
Your project folder should look exactly like this:

```text
uniprot_enzyme_scanner/
├── config/
│   ├── __init__.py        
│   ├── keywords.py        # Keyword-based objects (CATALYTIC_KW_MAP, EC_CLASSES, FUNCTIONAL_CATEGORIES)
│   ├── patterns.py        # Comment-based patterns (HARD_INACTIVE_PATTERN, etc.)
│   └── settings.py        # Settings (Paths, Batch sizes)
├── core/
│   ├── __init__.py        
│   ├── classifier.py      # Classifier logic (determine_enzyme_type, calculate_confidence...)
│   ├── models.py          # Dataclasses (EnzymeEntry)
│   └── parser.py          # XML Parsing (UniprotEntryParser)
├── utils/
│   ├── __init__.py        
│   ├── io.py              # Helper (smart_open, fast_xml_stream)
│   ├── download_XML.py    # Download latest released Uniprot XML datasets
│   └── obo_loader.py      # Gene Ontology DAG loader
├── main_extract.py        # XML Extractor
├── main_classify.py       # Uniprot annotations classifier
├── requirements.txt       # Dependencies
└── README.md              # Instructions
```

You are now ready to run. Start Stage 1 (main_extract.py), and while that churns through the XML, you can open core/classifier.py and start fine-tuning your logic for Stage 2.