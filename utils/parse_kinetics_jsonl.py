# utils/parse_kinetics_jsonl.py
import re
import argparse
import gzip
import json
import pandas as pd
import os
import glob
from multiprocessing import Pool, cpu_count, current_process, RLock
from tqdm import tqdm

def init_worker(tqdm_lock):
    tqdm.set_lock(tqdm_lock)

def normalize_value(val_str, param_type):
    """Normalizes KM to uM and kcat to sec(-1). Returns (value, unit)."""
    if not val_str: return "", ""
    
    match = re.match(r'^([~<>\s]*[\d\.,]+(?:[eE][\+\-]?\d+)?)\s*(.*)$', str(val_str).strip())
    if not match:
        return val_str, "" 
        
    try:
        clean_num = re.sub(r'[~<>\s]', '', match.group(1))
        num = float(clean_num)
        unit = match.group(2).strip()
        prefix = match.group(1).replace(clean_num, '').strip() 
        
        if param_type == "KM":
            if unit.lower() == "mm": return f"{prefix}{round(num * 1000, 4)}", "uM"
            if unit.lower() == "nm": return f"{prefix}{round(num / 1000, 4)}", "uM"
            if unit.lower() == "m": return f"{prefix}{round(num * 1000000, 4)}", "uM"
        elif param_type == "kcat":
            if unit.lower() == "min(-1)": return f"{prefix}{round(num / 60, 4)}", "sec(-1)"
            if unit.lower() == "h(-1)": return f"{prefix}{round(num / 3600, 4)}", "sec(-1)"
            
        return f"{prefix}{num}", unit
    except ValueError:
        return val_str, ""

def parse_kinetic_string(text, param_type):
    """Structurally parses kinetic strings using linguistic pivots and edge-case handlers."""
    results = []
    
    # Pre-processing: Clean PubMed and Linguistic Operators
    text = re.sub(r'\s*\(PubMed:\d+(?:,\s*PubMed:\d+)*\)', '', text, flags=re.IGNORECASE)
    text = re.sub(r'\babout\b', '~', text, flags=re.IGNORECASE)
    text = re.sub(r'\binferior to\b', '<', text, flags=re.IGNORECASE)
    
    if param_type == "kcat":
        sentences = [s.strip() for s in re.split(r'(?<=[a-zA-Z\)])\.\s+|(?<=[a-zA-Z\)])\.$', text) if s.strip()]
    else:
        sentences = [text.strip()]
        
    for sentence in sentences:
        sl = sentence.lower()
        if param_type == "kcat":
            if "kcat" not in sl or "kcat/km" in sl or "efficiency" in sl: 
                continue 

        # Ignore descriptive sentences lacking absolute kinetic parameters
        ignore_phrases = ["calculations are based", "fold increase", "orders of magnitude", "remains unchanged"]
        if any(phrase in sl for phrase in ignore_phrases):
            continue

        # 1a. Handle standard "respectively" (V1 and V2 with S1 and S2)
        resp_match = re.search(r"(?:(?:kcat|KM|Vmax)\s+(?:is|are|=)\s*)?(.*?)\s+and\s+(.*?)\s+(?:with|for|toward|against)\s+(.*?)\s+and\s+(.*?)\s+as substrates(?:,\s*respectively)?", sentence, flags=re.IGNORECASE)
        if resp_match:
            v1, v2, s1, s2 = resp_match.groups()
            results.append({"value": v1.strip(), "substrate": s1.strip(), "condition": "", "status": "Parsed"})
            results.append({"value": v2.strip(), "substrate": s2.strip(), "condition": "", "status": "Parsed"})
            continue
            
        # 1b. Handle apparent "respectively" (S1 and S2 is V1 and V2)
        resp_match_2 = re.search(r"(?:The\s+apparent\s+)?(?:kcat|KM|Vmax)\s+of\s+(.*?)\s+and\s+(.*?)(?:\s+degradation)?\s+(?:is|are|=)\s*(.*?)\s*,?\s*and\s*(.*?)\s*,?\s*respectively", sentence, flags=re.IGNORECASE)
        if resp_match_2:
            s1, s2, v1, v2 = resp_match_2.groups()
            results.append({"value": v1.strip(), "substrate": s1.strip(), "condition": "", "status": "Parsed"})
            results.append({"value": v2.strip(), "substrate": s2.strip(), "condition": "", "status": "Parsed"})
            continue
            
        # 1c. Handle multi-condition (at pH X is Y, at pH Z is W)
        multi_cond = re.search(r'^(?:kcat|KM|Vmax)\s+(at\s+pH\s+[\d\.]+)\s+(?:is|are|=)\s*([~<>\s]*[\d\.]+\s*[a-zA-Z\/\-\(\)\d]+)\s*,\s*(at\s+pH\s+[\d\.]+)\s+(?:is|are|=)\s*([~<>\s]*[\d\.]+\s*[a-zA-Z\/\-\(\)\d]+)', sentence, flags=re.IGNORECASE)
        if multi_cond:
            c1, v1, c2, v2 = multi_cond.groups()
            results.append({"value": v1.strip(), "substrate": "", "condition": c1.strip(), "status": "Parsed"})
            results.append({"value": v2.strip(), "substrate": "", "condition": c2.strip(), "status": "Parsed"})
            continue

        # Extract Condition explicitly at the end to simplify the string
        condition = ""
        cond_match = re.search(r'\(\s*(in the presence of|at |pH |with |when ).+?\)$', sentence, flags=re.IGNORECASE)
        if cond_match:
            condition = cond_match.group(0).strip('()')
            sentence = sentence[:cond_match.start()].strip() + " " + sentence[cond_match.end():].strip()

        # Clean leading prefixes
        clean_sentence = re.sub(r'^(?:The\s+)?(?:apparent\s+)?(?:kcat|KM|Vmax)\s+(?:is|are|=|of|:)\s+', '', sentence, flags=re.IGNORECASE)
        clean_sentence = re.sub(r'^(?:The\s+)?(?:apparent\s+)?(?:kcat|KM|Vmax)\s+', '', clean_sentence, flags=re.IGNORECASE)
        clean_sentence = clean_sentence.strip()

        # 2. Inverted Syntax (e.g., "for Tryptophan is 1.5 min(-1)")
        inv_match = re.search(r'^(?:for|toward|with|against)\s+(.+?)\s+(?:is|are|=)\s*(.+)$', clean_sentence, flags=re.IGNORECASE)
        if inv_match:
            results.append({"value": inv_match.group(2).strip(), "substrate": inv_match.group(1).strip(), "condition": condition, "status": "Parsed"})
            continue

        # 3. Standard Syntax (e.g., "1.5 uM for Tryptophan")
        std_match = re.search(r'^(.+?)\s+(for|toward|with|against)\s+(.+)$', clean_sentence, flags=re.IGNORECASE)
        if std_match:
            val_unit = std_match.group(1).strip()
            rest = std_match.group(3).strip()
            rest = re.sub(r',?\s*as substrates?$', '', rest, flags=re.IGNORECASE)
            results.append({"value": val_unit, "substrate": rest, "condition": condition, "status": "Parsed"})
            continue

        # 4. Fallback: No pivot word (e.g., just "1.5 min(-1)")
        val_match = re.match(r'^([~<>\s]*[\d\.,]+(?:[eE][\+\-]?\d+)?(?:x10\^\d+)?)\s*(.*)$', clean_sentence)
        if val_match:
            val = val_match.group(1) + (" " + val_match.group(2) if val_match.group(2) else "")
            results.append({"value": val.strip(), "substrate": "", "condition": condition, "status": "Parsed"})
        else:
            results.append({"value": sentence, "substrate": "", "condition": condition, "status": "Fallback"})
            
    return results

def process_file(filepath):
    parsed_rows = []
    filename = os.path.basename(filepath)
    worker_id = current_process()._identity[0] if current_process()._identity else 1
    
    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        total_lines = sum(1 for _ in f)

    with gzip.open(filepath, 'rt', encoding='utf-8') as f:
        for line in tqdm(f, total=total_lines, position=worker_id, leave=False, desc=filename[:18], mininterval=0.5):
            try:
                entry = json.loads(line)
            except json.JSONDecodeError:
                continue
            
            if not entry.get('features', {}).get('has_kinetics', False):
                continue
                
            accession = entry.get('meta', {}).get('accession', 'Unknown')
            text_blocks = entry.get('text', {})

            for p_type in ['KM', 'Vmax', 'kcat']:
                raw_data = text_blocks.get(f'kinetics_{p_type.lower()}', [])
                if not raw_data: continue
                
                raw_strings = [s.strip() for s in raw_data.split(' | ')] if isinstance(raw_data, str) else raw_data
                
                for raw_str in raw_strings:
                    if not raw_str.strip(): continue
                    
                    for data in parse_kinetic_string(raw_str.strip(), p_type):
                        norm_val, norm_unit = normalize_value(data.get('value'), p_type)
                        parsed_rows.append({
                            'accession': accession,
                            'param_type': p_type,
                            'raw_string': raw_str,
                            'parsed_value': data.get('value'),
                            'norm_value': norm_val,
                            'norm_unit': norm_unit,
                            'parsed_substrate': data.get('substrate'),
                            'parsed_condition': data.get('condition'),
                            'status': data.get('status')
                        })
    return parsed_rows

def main(input_path, output_csv, num_workers):
    files = glob.glob(os.path.join(input_path, '*.jsonl.gz')) if os.path.isdir(input_path) else [input_path]
    if not files: return print(f"No .jsonl.gz files found in {input_path}")
        
    print(f"Found {len(files)} files. Processing with {num_workers} workers...")
    print("\n") 
    
    all_rows = []
    tqdm_lock = RLock()
    
    with Pool(num_workers, initializer=init_worker, initargs=(tqdm_lock,)) as pool:
        for result in tqdm(pool.imap_unordered(process_file, files), total=len(files), desc="Total Progress", position=0, leave=True):
            all_rows.extend(result)

    print("\n\nParsing finished. Saving results to CSV...")
    out_df = pd.DataFrame(all_rows)
    out_df.to_csv(output_csv, index=False)
    
    total = len(out_df)
    failed = len(out_df[out_df['status'] == 'Fallback'])
    print(f"\nTotal parsed: {total} | Fallbacks: {failed} | Success Rate: {((total-failed)/total)*100:.4f}%")
    print(f"Detailed output saved to: {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', default="data/kinetics_report_detailed.csv")
    parser.add_argument('-w', '--workers', type=int, default=cpu_count())
    args = parser.parse_args()
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    main(args.input, args.output, args.workers)