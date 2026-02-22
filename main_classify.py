# main_classify.py

import argparse
import json
import csv
import io
import multiprocessing
import time
import itertools
import glob
import re
import os
from tqdm import tqdm

from core.classifier import EnzymeClassifier
from utils.io import smart_open
from utils.obo_loader import load_catalytic_go_ids
from config.settings import OBO_FILE
from config.keywords import CATALYTIC_KW_MAP

COLUMNS = [
    "accession", "id", "protein_name", 
    "reviewed", "protein_existence", "length", "fragment",
    "version", "date_modified", "date_created",
    "organism_id", "organism_name", "kingdom", "lineage", "gene_primary",
    "primary_class", "enzyme_type", "confidence", "functional_category",
    "ec_family", "ec_class", "ec_completeness",
    "flag_has_ec", "flag_has_active_site", "flag_has_kinetics", 
    "flag_has_reaction", "flag_has_experimental", "flag_is_pseudo",
    "ec", "ec_from_name", "all_ec",
    "rhea_ids", "keyword", "go_term",
    "kinetics_km", "kinetics_vmax", "kinetics_kcat",
    "reaction_text", "substrates", "products",
    "ft_act_site", "ft_binding", "ft_mutagenesis",
    "cc_catalytic_activity", "cc_cofactor", "cc_pathway", "cc_function",
    "cc_subunit", "cc_subcellular_location", "cc_induction", "cc_similarity",
    "cc_caution", "caution_activity_relevance", "cc_interaction",
    "cc_disease", "cc_phenotype",
    "xref_alphafolddb", "xref_brenda", "xref_reactome", 
    "xref_string", "xref_ensembl", "xref_kegg", "xref_unipathway"
]

worker_context = {}

def worker_init(valid_go_ids):
    """Initializes the classifier ONCE per worker process."""
    worker_context['classifier'] = EnzymeClassifier(valid_catalytic_go_ids=valid_go_ids)

def get_estimated_total():
    """Attempts to find release notes to get the exact total entries for accurate ETA."""
    total = 200573000 # Default fallback (approx 200M TrEMBL + 573k SProt)
    
    # Check current directory and data/ directory for relnotes
    note_files = glob.glob("relnotes*.txt") + glob.glob("data/relnotes*.txt")
    
    if note_files:
        try:
            with open(note_files[0], "r", encoding="utf-8", errors="ignore") as f:
                content = f.read()
                sp_match = re.search(r"Swiss-Prot:?\s*([\d,]+)\s+entries", content, re.IGNORECASE)
                tr_match = re.search(r"TrEMBL:?\s*([\d,]+)\s+entries", content, re.IGNORECASE)
                
                sp_count = int(sp_match.group(1).replace(",", "")) if sp_match else 573000
                tr_count = int(tr_match.group(1).replace(",", "")) if tr_match else 200000000
                total = sp_count + tr_count
        except Exception:
            pass
            
    return total

def classify_chunk(lines):
    """
    Worker: Decodes JSON, classifies, and formats the TSV text ITSELF.
    Returns a single string block to prevent IPC memory bloat.
    """
    classifier = worker_context['classifier']
    
    # Internal text buffer
    output_buffer = io.StringIO()
    writer = csv.DictWriter(output_buffer, fieldnames=COLUMNS, delimiter='\t', extrasaction='ignore')
    
    processed_count = 0

    for line in lines:
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            continue
            
        feats = entry.get('features', {})
        meta = entry.get('meta', {})
        text = entry.get('text', {})
        lists = entry.get('lists', {})
        
        kw_ids = lists.get('kw_ids', [])
        ec_list = lists.get('ec_numbers', [])
        feats['has_catalytic_kw'] = any(k in CATALYTIC_KW_MAP for k in kw_ids)
        
        primary_class, stats_flags = classifier.classify_entry(entry, feats)
        
        enzyme_type = primary_class
        conf = "None"
        if primary_class == "Enzyme": conf = "High"
        elif primary_class == "Suspected_Enzyme": conf = "Medium"
        elif primary_class == "Pseudo_Enzyme": conf = "Bronze"

        cat = classifier.assign_category(kw_ids, ec_list, primary_class, feats.get('has_text_mining', False))
        fam, cls, comp = classifier.get_ec_hierarchy(ec_list)

        row = {
            "accession": meta.get('accession', ''),
            "id": meta.get('id', ''),
            "protein_name": text.get('protein_name', ''),
            "reviewed": "true" if meta.get('reviewed') else "false",
            "protein_existence": meta.get('protein_existence', ''),
            "length": meta.get('length', '0'),
            "fragment": "true" if meta.get('is_fragment') else "false",
            "version": meta.get('version', ''),
            "date_modified": meta.get('date_modified', ''),
            "date_created": meta.get('date_created', ''),
            
            "organism_id": text.get('organism_id', ''),
            "organism_name": text.get('organism_name', ''),
            "kingdom": text.get('kingdom', ''),
            "lineage": text.get('lineage', ''),
            "gene_primary": text.get('gene_primary', ''),
            
            "primary_class": primary_class,
            "enzyme_type": enzyme_type,
            "confidence": conf,
            "functional_category": cat,
            "ec_family": fam,
            "ec_class": cls,
            "ec_completeness": comp,
            
            "flag_has_ec": str(stats_flags.get("flag_has_ec", False)).upper(),
            "flag_has_active_site": str(stats_flags.get("flag_has_active_site", False)).upper(),
            "flag_has_kinetics": str(stats_flags.get("flag_has_kinetics", False)).upper(),
            "flag_has_reaction": str(stats_flags.get("flag_has_reaction", False)).upper(),
            "flag_has_experimental": str(stats_flags.get("flag_has_experimental", False)).upper(),
            "flag_is_pseudo": str(stats_flags.get("flag_is_pseudo", False)).upper(),
            
            "ec": ";".join(lists.get('ec_db', [])),
            "ec_from_name": ";".join(lists.get('ec_from_name', [])),
            "all_ec": ";".join(ec_list),
            "rhea_ids": ";".join(lists.get('rhea_ids', [])),
            "keyword": ";".join(lists.get('kw_text', lists.get('kw_ids', []))),
            "go_term": ";".join(lists.get('go_text', lists.get('go_ids', []))),
            
            "kinetics_km": text.get('kinetics_km', ''),
            "kinetics_vmax": text.get('kinetics_vmax', ''),
            "kinetics_kcat": text.get('kinetics_kcat', ''),
            "reaction_text": text.get('catalytic_activity', ''),
            "substrates": ";".join(lists.get('substrates', [])),
            "products": ";".join(lists.get('products', [])),
            
            "ft_act_site": ";".join(lists.get('active_sites', [])),
            "ft_binding": ";".join(lists.get('binding_sites', [])),
            "ft_mutagenesis": ";".join(lists.get('mutagenesis', [])),
            
            "cc_catalytic_activity": text.get('catalytic_activity', ''),
            "cc_cofactor": text.get('cofactor', ''),
            "cc_pathway": text.get('pathway', ''),
            "cc_function": text.get('function', ''),
            "cc_subunit": text.get('subunit', ''),
            "cc_subcellular_location": text.get('subcellular_location', ''),
            "cc_induction": text.get('induction', ''),
            "cc_similarity": text.get('similarity', ''),
            "cc_caution": text.get('caution', ''),
            "cc_interaction": text.get('interaction', ''),
            "cc_disease": text.get('disease', ''),
            "cc_phenotype": text.get('phenotype', ''),
            "caution_activity_relevance": "true" if feats.get('caution_activity_relevance') else "false",
            
            "xref_alphafolddb": ";".join(lists.get('xref_af', [])),
            "xref_brenda": ";".join(lists.get('xref_brenda', [])),
            "xref_reactome": ";".join(lists.get('xref_reactome', [])),
            "xref_string": ";".join(lists.get('xref_string', [])),
            "xref_ensembl": ";".join(lists.get('xref_ensembl', [])),
            "xref_kegg": ";".join(lists.get('xref_kegg', [])),
            "xref_unipathway": ";".join(lists.get('xref_unipath', []))
        }
        
        writer.writerow(row)
        processed_count += 1
        
        # Free memory immediately
        del entry
        del line
        
    return (processed_count, output_buffer.getvalue())

def chunk_reader(filepath, chunk_size=2000):
    """Yields lists of text lines."""
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

    valid_go_ids = load_catalytic_go_ids(OBO_FILE)
    num_workers = max(1, multiprocessing.cpu_count() - 1)
    total_entries = get_estimated_total()
    
    # maxtasksperchild forces workers to restart, eliminating garbage collection leaks
    pool = multiprocessing.Pool(
        num_workers,
        initializer=worker_init,
        initargs=(valid_go_ids,),
        maxtasksperchild=100
    )
    
    print(f"[Start] Classifying from {args.input} with {num_workers} workers...")
    start_time = time.time()
    count = 0
    
    try:
        # We handle writing directly in the main process to eliminate Manager overhead
        with smart_open(args.output, 'wt') as f_out:
            f_out.write("\t".join(COLUMNS) + "\n")
            
            stream = chunk_reader(args.input, chunk_size=2000)
            
            # Super Batch size: Limit exactly how many chunks are pulled into memory at once
            # e.g., If 15 workers, read 45 chunks, process them fully, write to disk, THEN read more.
            super_batch_size = num_workers * 3  
            
            with tqdm(total=total_entries, desc="Classifying", unit=" proteins") as pbar:
                while True:
                    # Explicitly pull a limited number of chunks into memory
                    super_batch = list(itertools.islice(stream, super_batch_size))
                    
                    if not super_batch:
                        break # End of file
                        
                    # Calculate how many proteins are exactly in this super batch
                    entries_in_batch = sum(len(chunk) for chunk in super_batch)
                        
                    # Blocks here until all chunks in the super_batch are parsed
                    results = pool.map(classify_chunk, super_batch)
                    
                    # Write results straight to disk
                    for batch_count, text_block in results:
                        if batch_count > 0:
                            count += batch_count
                            f_out.write(text_block)
                            
                    pbar.update(entries_in_batch)
                    
                    # Free the variables so garbage collector can flush RAM
                    del super_batch
                    del results

    except KeyboardInterrupt:
        print("\n[Stop] Interrupted by user.")
    finally:
        pool.close()
        pool.join()

    print(f"\n[Done] Written {count} entries to {args.output}")
    print(f"[Time] {(time.time() - start_time)/60:.2f} minutes")

if __name__ == "__main__":
    main()