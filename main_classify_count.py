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
from collections import Counter
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

def worker_init(valid_go_ids, save_non_enzymes):
    """Initializes the classifier ONCE per worker process."""
    worker_context['classifier'] = EnzymeClassifier(valid_catalytic_go_ids=valid_go_ids)
    worker_context['save_non_enzymes'] = save_non_enzymes

def get_estimated_total():
    """Attempts to find release notes to get the exact total entries for accurate ETA."""
    total = 200573000
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
    Worker: Decodes JSON, classifies, writes enzyme rows to an internal buffer,
    and tallies strict combinatorial statistics.
    """
    classifier = worker_context['classifier']
    save_non_enzymes = worker_context['save_non_enzymes']
    
    enz_buffer = io.StringIO()
    non_enz_buffer = io.StringIO() if save_non_enzymes else None
    
    enz_writer = csv.DictWriter(enz_buffer, fieldnames=COLUMNS, delimiter='\t', extrasaction='ignore')
    non_enz_writer = csv.DictWriter(non_enz_buffer, fieldnames=COLUMNS, delimiter='\t', extrasaction='ignore') if save_non_enzymes else None
    
    stats = Counter()
    enz_count = 0
    non_enz_count = 0

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
        
        # -------------------------------------------------------------
        # TRACK COMBINATORIAL STATS (In Memory)
        # -------------------------------------------------------------
        fam, cls, comp = classifier.get_ec_hierarchy(ec_list)
        
        has_rxn = 1 if stats_flags.get("flag_has_reaction") else 0
        has_act = 1 if stats_flags.get("flag_has_active_site") else 0
        has_kin = 1 if stats_flags.get("flag_has_kinetics") else 0
        has_kw = 1 if feats.get("has_catalytic_kw") else 0
        has_exp = 1 if stats_flags.get("flag_has_experimental") else 0
        
        # Create a unique string signature for this exact combination of evidence
        combo_sig = f"{primary_class} | EC:{comp}, Rxn:{has_rxn}, ActSite:{has_act}, Kin:{has_kin}, KW:{has_kw}, Exp:{has_exp}"
        stats[combo_sig] += 1
        stats[f"TOTAL_{primary_class}"] += 1

        # -------------------------------------------------------------
        # FORMAT OUTPUT STRINGS
        # -------------------------------------------------------------
        if primary_class == "Non_Enzyme" and not save_non_enzymes:
            non_enz_count += 1
            del entry
            del line
            continue # SKIP formatting strings for non-enzymes! Massive speedup.

        conf = "None"
        if primary_class == "Enzyme": conf = "High"
        elif primary_class == "Suspected_Enzyme": conf = "Medium"
        elif primary_class == "Pseudo_Enzyme": conf = "Bronze"

        cat = classifier.assign_category(kw_ids, ec_list, primary_class, feats.get('has_text_mining', False))

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
            "enzyme_type": primary_class,
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
        
        if primary_class == "Non_Enzyme":
            non_enz_writer.writerow(row)
            non_enz_count += 1
        else:
            enz_writer.writerow(row)
            enz_count += 1
            
        del entry
        del line
        
    non_enz_text = non_enz_buffer.getvalue() if save_non_enzymes else ""
    return (enz_count, enz_buffer.getvalue(), non_enz_count, non_enz_text, stats)

def chunk_reader(filepath, chunk_size=2000):
    with smart_open(filepath, 'rt') as f:
        chunk = []
        for line in f:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        if chunk:
            yield chunk

def generate_stats_report(global_stats, out_file):
    """Processes the flat counters into a deep analytical JSON report with percentages."""
    report = { "Total_Entries": 0, "Classes": {} }
    
    # 1. Setup Base Classes
    for key, count in global_stats.items():
        if key.startswith("TOTAL_"):
            cls_name = key.split("TOTAL_")[1]
            report["Classes"][cls_name] = {
                "Total": count, 
                "Marginal_Features": Counter(), 
                "Combinations": {},
                "Combinations_Pct": {}
            }
            report["Total_Entries"] += count

    # 2. Parse Combinations & Calculate Marginals
    for key, count in global_stats.items():
        if " | " in key:
            cls_name, combo = key.split(" | ")
            report["Classes"][cls_name]["Combinations"][combo] = count
            
            # Tally Marginals (e.g., Total EC_Full regardless of other features)
            features = combo.split(", ")
            for feat in features:
                k, v = feat.split(":")
                if v != "0" and v != "None":
                    feat_name = f"Has_EC_{v}" if k == "EC" else f"Has_{k}"
                    report["Classes"][cls_name]["Marginal_Features"][feat_name] += count

    # 3. Calculate Percentages
    for cls_name, data in report["Classes"].items():
        total = data["Total"]
        if total == 0: continue
        
        # Marginal Percentages
        data["Marginal_Percentages"] = {
            k: f"{(v / total * 100):.2f}%" for k, v in data["Marginal_Features"].items()
        }
        
        # Sort Combinations by frequency and calculate percentage
        sorted_combos = sorted(data["Combinations"].items(), key=lambda x: x[1], reverse=True)
        data["Combinations"] = dict(sorted_combos)
        data["Combinations_Pct"] = {
            k: f"{(v / total * 100):.2f}%" for k, v in sorted_combos
        }

    with open(out_file, 'w') as f:
        json.dump(report, f, indent=4)
    print(f"\n[Stats] Generated detailed statistics report: {out_file}")


def main():
    parser = argparse.ArgumentParser(description="Stage 2: Classify Enzymes")
    parser.add_argument("--input", "-i", required=True, help="Input raw JSONL.gz")
    parser.add_argument("--out_enzymes", "-o", default="final_enzymes_v14.tsv.gz", help="Output enzymes TSV")
    parser.add_argument("--out_non_enzymes", "-n", default="final_non_enzymes_v14.tsv.gz", help="Output non-enzymes TSV")
    parser.add_argument("--out_stats", "-s", default="uniprot_combinatorial_stats_v14.json", help="Output Stats JSON")
    
    # NEW FLAG: By default, we skip writing Non-Enzymes to save massive time and RAM
    parser.add_argument("--save_non_enzymes", action="store_true", help="Pass this flag to actually save the non-enzymes to disk.")
    args = parser.parse_args()

    valid_go_ids = load_catalytic_go_ids(OBO_FILE)
    num_workers = max(1, multiprocessing.cpu_count() - 1)
    total_entries = get_estimated_total()
    
    pool = multiprocessing.Pool(
        num_workers,
        initializer=worker_init,
        initargs=(valid_go_ids, args.save_non_enzymes),
        maxtasksperchild=100
    )
    
    print(f"[Start] Classifying from {args.input} with {num_workers} workers...")
    if not args.save_non_enzymes:
        print("[Info] Non-Enzymes will be SKIPPED from TSV output to save I/O (Stats will still be calculated).")
        
    start_time = time.time()
    total_enz_written = 0
    total_non_enz_written = 0
    global_stats = Counter()
    
    try:
        f_enz = smart_open(args.out_enzymes, 'wt')
        f_enz.write("\t".join(COLUMNS) + "\n")
        
        f_non = None
        if args.save_non_enzymes:
            f_non = smart_open(args.out_non_enzymes, 'wt')
            f_non.write("\t".join(COLUMNS) + "\n")
            
        stream = chunk_reader(args.input, chunk_size=2000)
        super_batch_size = num_workers * 4  
        
        with tqdm(total=total_entries, desc="Classifying", unit=" proteins") as pbar:
            while True:
                super_batch = list(itertools.islice(stream, super_batch_size))
                if not super_batch: break
                    
                entries_in_batch = sum(len(chunk) for chunk in super_batch)
                results = pool.map(classify_chunk, super_batch)
                
                # Unpack and write
                for e_cnt, e_text, ne_cnt, ne_text, chunk_stats in results:
                    global_stats.update(chunk_stats)
                    
                    if e_cnt > 0:
                        total_enz_written += e_cnt
                        f_enz.write(e_text)
                    if ne_cnt > 0 and args.save_non_enzymes:
                        total_non_enz_written += ne_cnt
                        f_non.write(ne_text)
                        
                pbar.update(entries_in_batch)
                del super_batch
                del results

    except KeyboardInterrupt:
        print("\n[Stop] Interrupted by user.")
    finally:
        f_enz.close()
        if f_non: f_non.close()
        pool.close()
        pool.join()

    # Generate the beautiful stats report
    generate_stats_report(global_stats, args.out_stats)

    print(f"\n[Done] Time elapsed: {(time.time() - start_time)/60:.2f} minutes")
    print(f"[Result] Wrote {total_enz_written} Enzymes to {args.out_enzymes}")
    if args.save_non_enzymes:
        print(f"[Result] Wrote {total_non_enz_written} Non-Enzymes to {args.out_non_enzymes}")

if __name__ == "__main__":
    main()