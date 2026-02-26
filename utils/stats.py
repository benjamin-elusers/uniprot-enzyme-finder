# utils/stats.py
import json
from collections import Counter

def generate_stats_report(global_stats, out_file):
    """
    Processes the flat multiprocessing counters into a deep, 
    analytical JSON report with percentages.
    """
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
            # Catch any unexpected classes that weren't initialized
            if cls_name not in report["Classes"]:
                report["Classes"][cls_name] = {
                    "Total": 0, "Marginal_Features": Counter(), 
                    "Combinations": {}, "Combinations_Pct": {}
                }
            
            report["Classes"][cls_name]["Combinations"][combo] = count
            
            # Tally Marginals (e.g., Total EC_Full regardless of other features)
            features = combo.split(", ")
            for feat in features:
                k, v = feat.split(":")
                if v != "0" and v != "None":
                    feat_name = f"Has_EC_{v}" if k == "EC" else f"Has_{k}"
                    report["Classes"][cls_name]["Marginal_Features"][feat_name] += count

    # 3. Calculate Percentages and Sort
    for cls_name, data in report["Classes"].items():
        total = data["Total"]
        if total == 0: continue
        
        # Marginal Percentages
        data["Marginal_Percentages"] = {
            k: f"{(v / total * 100):.2f}%" for k, v in data["Marginal_Features"].items()
        }
        
        # Sort Combinations by absolute frequency and calculate percentage
        sorted_combos = sorted(data["Combinations"].items(), key=lambda x: x[1], reverse=True)
        data["Combinations"] = dict(sorted_combos)
        data["Combinations_Pct"] = {
            k: f"{(v / total * 100):.2f}%" for k, v in sorted_combos
        }

    # Write the formatted report to disk
    with open(out_file, 'w') as f:
        json.dump(report, f, indent=4)
        
    print(f"\n[Stats] Generated detailed statistics report: {out_file}")