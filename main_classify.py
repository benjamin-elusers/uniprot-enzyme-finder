# main_classify.py

import argparse
import json
import csv
import multiprocessing
import time
from tqdm import tqdm

from core.classifier import EnzymeClassifier
from utils.io import smart_open

# Output Columns for the TSV
COLUMNS = [
    "accession", "id", "protein_name", 
    "enzyme_type", "confidence", "functional_category",
    "is_inactive", "ec_numbers", "ec_family", 
    "catalytic_activity", "cofactor", "activity_regulation",
    "active_sites", "binding_sites", 
    "keywords", "go_terms"
]

def classify_chunk(lines):
    """
    Worker: Decodes JSON lines and applies classification logic.
    """
    classifier = EnzymeClassifier() # Lightweight instantiation
    results = []

    for line in lines:
        entry = json.loads(line)
        
        # 1. Flatten Features for Classifier
        # We merge the boolean flags and the list existence into one dict
        feats = entry['features']
        feats['has_ec'] = bool(entry['lists'].get('ec_numbers'))
        feats['has_catalytic_kw'] = any(k in classifier.context_kw_map for k in entry['lists'].get('kw_ids', [])) 
        # Note: You might need to import CATALYTIC_KW_MAP here if classifier doesn't expose it publically, 
        # or rely on the classifier to check IDs. For now, let's assume the classifier handles raw lists if we pass them.
        
        # Actually, let's just pass the whole 'entry' to a helper or prep it here
        # PREP:
        kw_ids = entry['lists'].get('kw_ids', [])
        ec_list = entry['lists'].get('ec_numbers', [])
        
        # Check catalytic kws
        from config.keywords import CATALYTIC_KW_MAP
        feats['has_catalytic_kw'] = any(k in CATALYTIC_KW_MAP for k in kw_ids)
        
        # 2. Determine Type
        e_type = classifier.determine_enzyme_type(feats)
        
        # 3. Filter? (Optional: Skip "Non_Enzyme_Context" if you only want strict enzymes)
        if not e_type:
            continue

        # 4. Assign Category
        cat = classifier.assign_category(
            kw_ids, ec_list, e_type, feats.get('has_text_mining')
        )

        # 5. Calculate Confidence
        conf = classifier.calculate_confidence(
            e_type, ec_list, [], # ec_name isn't separated in raw yet, passing empty for now or fix in Parser
            [], # evidence_types not calculated in stage 2 yet, can add logic
            entry['meta']['reviewed'], 
            entry['meta']['protein_existence'],
            False, # is_fragment (need to add to parser!)
            False, # has_experimental (need to add to parser!)
            feats.get('has_kinetics'),
            False, # caution
            feats.get('has_active_site')
        )
        
        # 6. Build Row
        row = {
            "accession": entry['meta']['accession'],
            "id": entry['meta']['id'],
            "protein_name": entry['text'].get('protein_name', ''),
            "enzyme_type": e_type,
            "confidence": conf,
            "functional_category": cat,
            "is_inactive": feats.get('is_inactive_text'),
            "ec_numbers": ";".join(ec_list),
            "catalytic_activity": entry['text'].get('catalytic_activity', ''),
            "active_sites": ";".join(entry['lists'].get('active_sites', [])),
            "keywords": ";".join(kw_ids),
        }
        results.append(row)
        
    return results

def chunk_reader(filepath, chunk_size=10000):
    """Yields chunks of lines from the JSONL file."""
    with smart_open(filepath, 'rt') as f:
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk

def main():
    parser = argparse.ArgumentParser(description="Stage 2: Classify Enzymes")
    parser.add_argument("--input", "-i", required=True, help="Input raw JSONL.gz")
    parser.add_argument("--output", "-o", default="final_enzymes_v14.tsv.gz", help="Output TSV.gz")
    args = parser.parse_args()

    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    
    print(f"[Start] Classifying from {args.input}...")
    start = time.time()
    
    count = 0
    with smart_open(args.output, 'wt') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=COLUMNS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        
        stream = chunk_reader(args.input)
        
        for batch_rows in tqdm(pool.imap_unordered(classify_chunk, stream), desc="Classifying"):
            if batch_rows:
                count += len(batch_rows)
                writer.writerows(batch_rows)

    pool.close()
    pool.join()
    print(f"[Done] Written {count} enzymes to {args.output}")
    print(f"[Time] {(time.time() - start)/60:.2f} minutes")

if __name__ == "__main__":
    main()