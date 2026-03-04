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

    def _determine_primary_class(self, is_inactive, high_prob, med_prob, has_binding_site):
        """Decision tree logic isolated for modularity."""
        if is_inactive:
            # If it has enzyme structure/clues but is flagged inactive -> Pseudo
            if high_prob or med_prob or has_binding_site:
                return "Pseudo_Enzyme"
            return "Non_Enzyme"
        
        if high_prob:
            return "Enzyme"
        elif med_prob:
            return "Suspected_Enzyme"
        
        return "Non_Enzyme"

    def _gather_evidence(self, entry_data, features, catalytic_gos, protein_name):
        """Evaluates entry features and returns explicitly matched evidence fields and values."""
        ev_fields, ev_vals = [], []
        
        if entry_data['lists'].get('ec_numbers'): 
            ev_fields.append("EC_Number")
            ev_vals.append(";".join(entry_data['lists']['ec_numbers']))
            
        if features.get('has_reaction'): 
            ev_fields.append("Reaction_Text")
            ev_vals.append("Reaction found")
            
        if entry_data['lists'].get('enzyme_dbs'): 
            ev_fields.append("Enzyme_DB")
            ev_vals.append(";".join(entry_data['lists']['enzyme_dbs']))
            
        if entry_data['lists'].get('rhea_ids'): 
            ev_fields.append("Rhea_ID")
            ev_vals.append(";".join(entry_data['lists']['rhea_ids']))
            
        if catalytic_gos: 
            ev_fields.append("GO_Catalytic")
            ev_vals.append(";".join(catalytic_gos))
            
        if features.get('has_active_site'): 
            ev_fields.append("Active_Site")
            ev_vals.append("Active site annotated")
            
        if features.get('has_catalytic_kw'): 
            ev_fields.append("Keyword_Catalytic")
            ev_vals.append("Catalytic keyword matched")
            
        if protein_name and ENZYME_NAME_ENDINGS.search(protein_name): 
            ev_fields.append("Protein_Name")
            ev_vals.append(protein_name)
            
        if features.get('has_biosynth_text'): 
            ev_fields.append("Text_Mining")
            ev_vals.append("Biosynthesis context")

        return ev_fields, ev_vals


    def classify_entry(self, entry_data, features):
        """
        Main decision tree. Returns the primary category and a dictionary 
        of exact boolean flags for plotting downstream.
        """
        # 1. Base Extractions & Setup
        caution_text = entry_data['text'].get('caution', '').lower()
        is_inactive_text = features.get('is_inactive_text', False) or "inactive" in caution_text or "catalytically dead" in caution_text
        
        protein_name = entry_data['text'].get('protein_name', '')
        has_enzyme_name = bool(ENZYME_NAME_ENDINGS.search(protein_name))
        
        go_ids = entry_data['lists'].get('go_ids', [])
        catalytic_gos = [go for go in go_ids if go in self.valid_catalytic_go_ids] if self.valid_catalytic_go_ids else []

        # 2. Probability Tiering
        high_prob = any([
            bool(entry_data['lists'].get('ec_numbers')),
            features.get('has_reaction', False),
            bool(entry_data['lists'].get('enzyme_dbs')),
            bool(entry_data['lists'].get('rhea_ids')),
            bool(catalytic_gos)
        ])

        med_prob = any([
            features.get('has_active_site', False),
            features.get('has_catalytic_kw', False),
            features.get('has_biosynth_text', False),
            has_enzyme_name
        ])
        
        has_binding_site = features.get('has_binding_site', False)

        # 3. Determine Class (Modular call)
        primary_class = self._determine_primary_class(
            is_inactive=is_inactive_text, 
            high_prob=high_prob, 
            med_prob=med_prob, 
            has_binding_site=has_binding_site
        )

        # 4. Gather Evidence Strings (Modular call)
        ev_fields, ev_vals = self._gather_evidence(
            entry_data=entry_data, 
            features=features, 
            catalytic_gos=catalytic_gos, 
            protein_name=protein_name
        )

        # 5. Build Stats Flags
        stats_flags = {
            "flag_has_ec": bool(entry_data['lists'].get('ec_numbers')),
            "flag_has_active_site": features.get('has_active_site', False),
            "flag_has_kinetics": features.get('has_kinetics', False),
            "flag_has_reaction": features.get('has_reaction', False),
            "flag_has_experimental": features.get('has_experimental_evidence', False),
            "flag_is_pseudo": (primary_class == "Pseudo_Enzyme"),
            "evidence_field": "|".join(ev_fields) if ev_fields else "None",
            "evidence_value": "|".join(ev_vals) if ev_vals else "None"
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