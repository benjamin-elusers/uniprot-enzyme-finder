# core/classifier.py
import re
from config.keywords import FUNCTIONAL_CATEGORIES, EC_CLASSES

class EnzymeClassifier:
    def __init__(self):
        # Flatten Context Map for O(1) lookups
        self.context_kw_map = {}
        for category, kws in FUNCTIONAL_CATEGORIES.items():
            for kw in kws:
                self.context_kw_map[kw] = category

    def determine_enzyme_type(self, features):
        """
        Central logic to determine if an entry is an Active Enzyme, 
        a Pseudo-Enzyme, or a Non-Enzyme Context match.
        
        Args:
            features (dict): Dictionary containing boolean flags and lists 
                             (e.g., has_ec, is_inactive, ft_act_site, etc.)
        """
        # Unpack features for readability
        is_inactive_text = features.get('is_inactive_text', False)
        has_ec = features.get('has_ec', False)
        has_catalytic_kw = features.get('has_catalytic_kw', False)
        has_catalytic_go = features.get('has_catalytic_go', False)
        has_kinetics = features.get('has_kinetics', False)
        has_active_site = features.get('has_active_site', False) # New Feature
        has_binding_site = features.get('has_binding_site', False) # New Feature
        has_text_mining = features.get('has_text_mining', False)
        has_context_clue = features.get('has_context_clue', False)
        
        # Strong indicators of "Enzymatic Nature" (Active or Inactive)
        is_structurally_enzyme = (
            has_ec or has_catalytic_kw or has_catalytic_go or 
            has_kinetics or has_active_site
        )

        # --- LOGIC BRANCH 1: INACTIVE HANDLING ---
        if is_inactive_text:
            # It says "Inactive", but is it a Pseudo-Enzyme or just irrelevant?
            # A Pseudo-Enzyme usually retains the fold (Active Site / Binding Site) 
            # or has a specific inactive EC (n.n.n.n).
            if is_structurally_enzyme or has_binding_site:
                return "Pseudo_Enzyme"
            else:
                return "Inactive_Non_Enzyme"

        # --- LOGIC BRANCH 2: ACTIVE ENZYMES ---
        if is_structurally_enzyme:
            if has_ec: 
                return "Explicit_EC"
            elif has_catalytic_kw: 
                return "Explicit_KW"
            elif has_kinetics:
                return "Explicit_Kinetics"
            elif has_active_site:
                # Active site but no EC/KW? Likely a known enzyme missing annotations.
                return "Implicit_Structure"
            elif has_text_mining and not has_ec: 
                return "Text_Inferred"
            else: 
                return "Implicit"

        # --- LOGIC BRANCH 3: CONTEXT MATCHES ---
        if has_context_clue:
            # Matches "Toxin" or "Antibiotic" but NO catalytic evidence
            return "Non_Enzyme_Context"

        return None

    def assign_category(self, kw_ids, ec_numbers, enzyme_type, has_sec_met_text):
        """
        Determines Functional Category (Context).
        """
        if "Inactive" in enzyme_type and enzyme_type != "Pseudo_Enzyme": 
            return "Inactive_Non_Enzyme"

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


    def calculate_confidence(self, enzyme_type, ec_db, ec_name, evidence_types, 
                             is_reviewed, pe_level, is_fragment, 
                             has_experimental, has_kinetics, caution_relevant,
                             has_active_site):
        """
        Assigns confidence level (Certified, Gold, Silver, Bronze).
        """
        if not enzyme_type or "Non_Enzyme" in enzyme_type: return "None"
        if "Inactive" in enzyme_type: return "Inactive"
        if enzyme_type == "Pseudo_Enzyme": return "Bronze"
        
        has_full_ec_db = any(ec.count('.') == 3 and '-' not in ec for ec in ec_db)
        has_full_ec_name = any(ec.count('.') == 3 and '-' not in ec for ec in ec_name)
        has_strong_context = any(e in ["KW_Catalytic", "Database", "GO_Catalytic"] for e in evidence_types)

        base = "Bronze"

        # --- TIER 1: CERTIFIED ---
        # Verified Structure AND Verified Activity
        if has_full_ec_db and (has_experimental or has_kinetics) and not caution_relevant:
            base = "Certified"
        
        # --- TIER 2: GOLD ---
        # High confidence prediction OR Reviewed with Name/Structure
        elif has_full_ec_db:
            base = "Gold"
        elif has_full_ec_name and is_reviewed:
            base = "Gold"
        elif is_reviewed and (has_strong_context or has_active_site):
            base = "Gold"
        
        # --- TIER 3: SILVER ---
        # Inference / Partially unreviewed
        elif (ec_db or ec_name): 
            base = "Silver"
        elif has_strong_context or enzyme_type == "Text_Inferred":
            base = "Silver"
        elif has_active_site:
            # Active site mapped >= Silver
            base = "Silver"
            
        # --- DOWNGRADES ---
        # 1. Caution Flag
        if caution_relevant: return "Bronze"
        
        # 2. Fragment / Uncertain Existence
        # Note: If it has an Active Site, we might be more lenient, 
        # but generally fragments shouldn't be Certified/Gold.
        if (pe_level == "5" or is_fragment):
            if base in ["Certified", "Gold"]:
                base = "Silver" # Downgrade to Silver instead of Bronze if strong evidence existed
            
        return base

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

