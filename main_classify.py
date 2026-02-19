# main_classify.py

import argparse
import json
import csv
import multiprocessing
import time
from tqdm import tqdm

from core.classifier import EnzymeClassifier
from utils.io import smart_open
from utils.obo_loader import load_catalytic_go_ids
from config.settings import OBO_FILE
from config.keywords import CATALYTIC_KW_MAP

# =============================================================================
#   OUTPUT COLUMNS (V13 Format + New Plotting Flags)
# =============================================================================
COLUMNS = [
    # Identifiers & Meta
    "accession", "id", "protein_name", 
    "reviewed", "protein_existence", "length", "fragment",
    "version", "date_modified", "date_created",
    
    # Taxonomy
    "organism_id", "organism_name", "kingdom", "lineage", "gene_primary",
    
    # Classification & Stats (New Plotting Data)
    "primary_class", "enzyme_type", "confidence", "functional_category",
    "ec_family", "ec_class", "ec_completeness",
    
    # Boolean Flags for Plotting
    "flag_has_ec", "flag_has_active_site", "flag_has_kinetics", 
    "flag_has_reaction", "flag_has_experimental", "flag_is_pseudo",
    
    # EC & External IDs
    "ec", "ec_from_name", "all_ec",
    "rhea_ids", "keyword", "go_term",
    
    # Kinetics & Reactions
    "kinetics_km", "kinetics_vmax", "kinetics_kcat",
    "reaction_text", "substrates", "products",
    
    # Sequence Features
    "ft_act_site", "ft_binding", "ft_mutagenesis",
    
    # Comments & Context
    "cc_catalytic_activity", "cc_cofactor", "cc_pathway", "cc_function",
    "cc_subunit", "cc_subcellular_location", "cc_induction", "cc_similarity",
    "cc_caution", "caution_activity_relevance", "cc_interaction",
    "cc_disease", "cc_phenotype",
    
    # Database Cross-References
    "xref_alphafolddb", "xref_brenda", "xref_reactome", 
    "xref_string", "xref_ensembl", "xref_kegg", "xref_unipathway"
]

def classify_chunk(args):
    """
    Worker: Decodes JSON lines, applies classification logic, and flattens 
    the nested dictionary into the final V13 TSV row format.
    """
    lines, valid_go_ids = args
    classifier = EnzymeClassifier(valid_catalytic_go_ids=valid_go_ids) 
    results = []

    for line in lines:
        try:
            entry = json.loads(line)
        except json.JSONDecodeError:
            continue
            
        feats = entry.get('features', {})
        meta = entry.get('meta', {})
        text = entry.get('text', {})
        lists = entry.get('lists', {})
        
        # 1. Prep Lists for Classification
        kw_ids = lists.get('kw_ids', [])
        ec_list = lists.get('ec_numbers', [])
        feats['has_catalytic_kw'] = any(k in CATALYTIC_KW_MAP for k in kw_ids)
        
        # 2. Apply Classification Logic
        primary_class, stats_flags = classifier.classify_entry(entry, feats)
        
        # Determine internal legacy type/confidence
        enzyme_type = primary_class  # Map to standard enum
        conf = "None"
        if primary_class == "Enzyme": conf = "High"
        elif primary_class == "Suspected_Enzyme": conf = "Medium"
        elif primary_class == "Pseudo_Enzyme": conf = "Bronze"

        # Determine Context Category
        cat = classifier.assign_category(
            kw_ids, ec_list, primary_class, feats.get('has_text_mining', False)
        )
        
        # Hierarchy
        fam, cls, comp = classifier.get_ec_hierarchy(ec_list)

        # 3. Flatten and Map to TSV Row
        row = {
            # Meta
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
            
            # Taxonomy
            "organism_id": text.get('organism_id', ''),
            "organism_name": text.get('organism_name', ''),
            "kingdom": text.get('kingdom', ''),
            "lineage": text.get('lineage', ''),
            "gene_primary": text.get('gene_primary', ''),
            
            # Classification
            "primary_class": primary_class,
            "enzyme_type": enzyme_type,
            "confidence": conf,
            "functional_category": cat,
            "ec_family": fam,
            "ec_class": cls,
            "ec_completeness": comp,
            
            # Plotting Flags
            "flag_has_ec": str(stats_flags.get("flag_has_ec", False)).upper(),
            "flag_has_active_site": str(stats_flags.get("flag_has_active_site", False)).upper(),
            "flag_has_kinetics": str(stats_flags.get("flag_has_kinetics", False)).upper(),
            "flag_has_reaction": str(stats_flags.get("flag_has_reaction", False)).upper(),
            "flag_has_experimental": str(stats_flags.get("flag_has_experimental", False)).upper(),
            "flag_is_pseudo": str(stats_flags.get("flag_is_pseudo", False)).upper(),
            
            # IDs & Lists
            "ec": ";".join(lists.get('ec_db', [])),
            "ec_from_name": ";".join(lists.get('ec_from_name', [])),
            "all_ec": ";".join(ec_list),
            "rhea_ids": ";".join(lists.get('rhea_ids', [])),
            "keyword": ";".join(lists.get('kw_text', lists.get('kw_ids', []))),
            "go_term": ";".join(lists.get('go_text', lists.get('go_ids', []))),
            
            # Kinetics
            "kinetics_km": text.get('kinetics_km', ''),
            "kinetics_vmax": text.get('kinetics_vmax', ''),
            "kinetics_kcat": text.get('kinetics_kcat', ''),
            "reaction_text": text.get('catalytic_activity', ''),
            "substrates": ";".join(lists.get('substrates', [])),
            "products": ";".join(lists.get('products', [])),
            
            # Sequence Features
            "ft_act_site": ";".join(lists.get('active_sites', [])),
            "ft_binding": ";".join(lists.get('binding_sites', [])),
            "ft_mutagenesis": ";".join(lists.get('mutagenesis', [])),
            
            # Comments
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
            
            # XRefs
            "xref_alphafolddb": ";".join(lists.get('xref_af', [])),
            "xref_brenda": ";".join(lists.get('xref_brenda', [])),
            "xref_reactome": ";".join(lists.get('xref_reactome', [])),
            "xref_string": ";".join(lists.get('xref_string', [])),
            "xref_ensembl": ";".join(lists.get('xref_ensembl', [])),
            "xref_kegg": ";".join(lists.get('xref_kegg', [])),
            "xref_unipathway": ";".join(lists.get('xref_unipath', []))
        }
        
        results.append(row)
        
    return results

def chunk_reader(filepath, valid_go_ids, chunk_size=10000):
    """Yields chunks of JSON lines and the GO validation set for multiprocessing."""
    with smart_open(filepath, 'rt') as f:
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield (chunk, valid_go_ids)
                chunk = []
        if chunk:
            yield (chunk, valid_go_ids)

def main():
    parser = argparse.ArgumentParser(description="Stage 2: Classify Enzymes")
    parser.add_argument("--input", "-i", required=True, help="Input raw JSONL.gz")
    parser.add_argument("--output", "-o", default="final_enzymes_v14.tsv.gz", help="Output TSV.gz")
    args = parser.parse_args()

    # Load GO Terms once to pass to workers to prevent massive memory duplication
    valid_go_ids = load_catalytic_go_ids(OBO_FILE)

    pool = multiprocessing.Pool(max(1, multiprocessing.cpu_count() - 1))
    
    print(f"[Start] Classifying from {args.input}...")
    start = time.time()
    
    count = 0
    with smart_open(args.output, 'wt') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=COLUMNS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        
        stream = chunk_reader(args.input, valid_go_ids)
        
        for batch_rows in tqdm(pool.imap_unordered(classify_chunk, stream), desc="Classifying"):
            if batch_rows:
                count += len(batch_rows)
                writer.writerows(batch_rows)

    pool.close()
    pool.join()
    print(f"[Done] Written {count} entries to {args.output}")
    print(f"[Time] {(time.time() - start)/60:.2f} minutes")

if __name__ == "__main__":
    main()