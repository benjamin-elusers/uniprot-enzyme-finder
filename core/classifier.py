# core/classifier.py

import re
from config.keywords import FUNCTIONAL_CATEGORIES, EC_CLASSES
from config.patterns import ENZYME_NAME_ENDINGS

class EnzymeClassifier:
    def __init__(self, valid_catalytic_go_ids=None):
        self.valid_catalytic_go_ids = valid_catalytic_go_ids or set()
        
        # Flatten Context Map for O(1) lookups
        self.context_kw_map = {}
        for category, kws in FUNCTIONAL_CATEGORIES.items():
            for kw in kws:
                self.context_kw_map[kw] = category

    def classify_entry(self, entry_data, features):
        """
        Main decision tree. Returns the primary category and a dictionary 
        of exact boolean flags for plotting downstream.
        """
        # Unpack evidence tiers
        is_inactive_text = features.get('is_inactive_text', False)
        caution_text = entry_data['text'].get('caution', '').lower()
        if "inactive" in caution_text or "catalytically dead" in caution_text:
            is_inactive_text = True
            
        # TIER 1: High Probability (Obvious)
        has_ec = bool(entry_data['lists'].get('ec_numbers'))
        has_catalytic_activity = features.get('has_reaction', False)
        has_enzyme_dbs = bool(entry_data['lists'].get('enzyme_dbs'))
        has_rhea = bool(entry_data['lists'].get('rhea_ids'))
        
        # Check GO Catalytic descendants
        go_ids = entry_data['lists'].get('go_ids', [])
        has_catalytic_go = any(go in self.valid_catalytic_go_ids for go in go_ids) if self.valid_catalytic_go_ids else False
        
        high_prob = any([has_ec, has_catalytic_activity, has_enzyme_dbs, has_rhea, has_catalytic_go])

        # TIER 2: Medium Probability (Clues)
        has_active_site = features.get('has_active_site', False)
        has_catalytic_kw = features.get('has_catalytic_kw', False)
        has_biosynth_text = features.get('has_biosynth_text', False)
        
        protein_name = entry_data['text'].get('protein_name', '')
        has_enzyme_name = bool(ENZYME_NAME_ENDINGS.search(protein_name))
        
        med_prob = any([has_active_site, has_catalytic_kw, has_biosynth_text, has_enzyme_name])

        # TIER 3: Low Probability / Context (Binding, pathways, etc. - mostly tracked for context)
        has_binding_site = features.get('has_binding_site', False)

        # --- DECISION TREE ---
        primary_class = "Non_Enzyme"
        
        if is_inactive_text:
            # If it has enzyme structure/clues but is flagged inactive -> Pseudo
            if high_prob or med_prob or has_binding_site:
                primary_class = "Pseudo_Enzyme"
            else:
                primary_class = "Non_Enzyme"
        elif high_prob:
            primary_class = "Enzyme"
        elif med_prob:
            primary_class = "Suspected_Enzyme"
        else:
            primary_class = "Non_Enzyme"

        # --- BOOLEAN FLAGS FOR STATS ---
        stats_flags = {
            "flag_has_ec": has_ec,
            "flag_has_active_site": has_active_site,
            "flag_has_kinetics": features.get('has_kinetics', False),
            "flag_has_reaction": has_catalytic_activity,
            "flag_has_experimental": features.get('has_experimental_evidence', False),
            "flag_is_pseudo": (primary_class == "Pseudo_Enzyme")
        }

        return primary_class, stats_flags

    def assign_category(self, kw_ids, ec_numbers, primary_class, has_sec_met_text):
        """
        Determines Functional Category (Context).
        """
        # 1. Secondary Metabolism (Highest Priority)
        if has_sec_met_text: return "Secondary_Metabolism"
        for kw in kw_ids:
            if self.context_kw_map.get(kw) == "Secondary_Metabolism":
                return "Secondary_Metabolism"
        
        # 2. Info Processing (DNA/RNA/Translation)
        for ec in ec_numbers:
            # 2.7.7 (Nucleotidyltransferases), 6.1.1 (Ligases forming aminoacyl-tRNA)
            if ec.startswith("2.7.7") or ec.startswith("6.1.1"): return "Info_Processing"
            # 3.1.2x (Thioesterases) - Specific check for exonucleases/etc often in 3.1.x
            if ec.startswith("3.1.2") and len(ec) > 5 and ec[4] in '6789': return "Info_Processing"
        for kw in kw_ids:
            if kw in self.context_kw_map and self.context_kw_map[kw] == "Info_Processing": return "Info_Processing"

        # 3. Signaling (Kinases, Phosphatases)
        for ec in ec_numbers:
            if ec.startswith("3.1.3") or (ec.startswith("2.7") and not ec.startswith("2.7.7")): 
                return "Signaling"
        for kw in kw_ids:
            if self.context_kw_map.get(kw) == "Signaling": return "Signaling"

        # 4. Core Metabolism
        for kw in kw_ids:
            if self.context_kw_map.get(kw) == "Core_Metabolism":
                return "Core_Metabolism"
            
        return "Unclassified"

    def get_ec_hierarchy(self, ec_numbers):
        if not ec_numbers: return "Non-Enzyme", "", "None"
        primary_ec = ec_numbers[0]
        # Clean potential dirty strings
        clean_ec = primary_ec.replace("-", "").strip()
        parts = clean_ec.split(".")
        
        completeness = "None"
        if "-" not in primary_ec and primary_ec.count(".") == 3: completeness = "Full"
        elif primary_ec.count(".") >= 1: completeness = "Partial"
        
        if not parts or not parts[0].isdigit(): return "Unclassified", "", completeness
        
        family_name = EC_CLASSES.get(parts[0], "Unclassified")
        return family_name, parts[0], completeness