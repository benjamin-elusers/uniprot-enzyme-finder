# core/classifier.py

import re
import json
from config.keywords import FUNCTIONAL_CATEGORIES, EC_CLASSES
from config.patterns import ENZYME_NAME_ENDINGS
from config.evidences import ECO_STRICT_EXP, ECO_HTP_EXP, ECO_CURATED, ECO_AUTOMATED

class EnzymeClassifier:
    def __init__(self, valid_catalytic_go_ids=None):
        self.valid_catalytic_go_ids = valid_catalytic_go_ids or set()
        
        # Flatten Context Map for O(1) lookups
        self.context_kw_map = {}
        for category, kws in FUNCTIONAL_CATEGORIES.items():
            for kw in kws:
                self.context_kw_map[kw] = category

    def _parse_kinetic_string(self, text, param_type):
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

            ignore_phrases = ["calculations are based", "fold increase", "orders of magnitude", "remains unchanged"]
            if any(phrase in sl for phrase in ignore_phrases):
                continue

            # 1a. Handle standard "respectively" (V1 and V2 with S1 and S2)
            resp_match = re.search(r"(?:(?:kcat|KM|Vmax)\s+(?:is|are|=)\s*)?(.*?)\s+and\s+(.*?)\s+(?:with|for|toward|against)\s+(.*?)\s+and\s+(.*?)\s+as substrates(?:,\s*respectively)?", sentence, flags=re.IGNORECASE)
            if resp_match:
                v1, v2, s1, s2 = resp_match.groups()
                results.append({"value": v1.strip(), "substrate": s1.strip(), "condition": ""})
                results.append({"value": v2.strip(), "substrate": s2.strip(), "condition": ""})
                continue
                
            # 1b. Handle apparent "respectively" (S1 and S2 is V1 and V2)
            resp_match_2 = re.search(r"(?:The\s+apparent\s+)?(?:kcat|KM|Vmax)\s+of\s+(.*?)\s+and\s+(.*?)(?:\s+degradation)?\s+(?:is|are|=)\s*(.*?)\s*,?\s*and\s*(.*?)\s*,?\s*respectively", sentence, flags=re.IGNORECASE)
            if resp_match_2:
                s1, s2, v1, v2 = resp_match_2.groups()
                results.append({"value": v1.strip(), "substrate": s1.strip(), "condition": ""})
                results.append({"value": v2.strip(), "substrate": s2.strip(), "condition": ""})
                continue
                
            # 1c. Handle multi-condition (at pH X is Y, at pH Z is W)
            multi_cond = re.search(r'^(?:kcat|KM|Vmax)\s+(at\s+pH\s+[\d\.]+)\s+(?:is|are|=)\s*([~<>\s]*[\d\.]+\s*[a-zA-Z\/\-\(\)\d]+)\s*,\s*(at\s+pH\s+[\d\.]+)\s+(?:is|are|=)\s*([~<>\s]*[\d\.]+\s*[a-zA-Z\/\-\(\)\d]+)', sentence, flags=re.IGNORECASE)
            if multi_cond:
                c1, v1, c2, v2 = multi_cond.groups()
                results.append({"value": v1.strip(), "substrate": "", "condition": c1.strip()})
                results.append({"value": v2.strip(), "substrate": "", "condition": c2.strip()})
                continue

            # Extract Condition explicitly
            condition = ""
            cond_match = re.search(r'\(\s*(in the presence of|at |pH |with |when ).+?\)$', sentence, flags=re.IGNORECASE)
            if cond_match:
                condition = cond_match.group(0).strip('()')
                sentence = sentence[:cond_match.start()].strip() + " " + sentence[cond_match.end():].strip()

            clean_sentence = re.sub(r'^(?:The\s+)?(?:apparent\s+)?(?:kcat|KM|Vmax)\s+(?:is|are|=|of|:)\s+', '', sentence, flags=re.IGNORECASE)
            clean_sentence = re.sub(r'^(?:The\s+)?(?:apparent\s+)?(?:kcat|KM|Vmax)\s+', '', clean_sentence, flags=re.IGNORECASE).strip()

            # 2. Inverted Syntax (e.g., "for Tryptophan is 1.5 min(-1)")
            inv_match = re.search(r'^(?:for|toward|with|against)\s+(.+?)\s+(?:is|are|=)\s*(.+)$', clean_sentence, flags=re.IGNORECASE)
            if inv_match:
                results.append({"value": inv_match.group(2).strip(), "substrate": inv_match.group(1).strip(), "condition": condition})
                continue

            # 3. Standard Syntax (e.g., "1.5 uM for Tryptophan")
            std_match = re.search(r'^(.+?)\s+(for|toward|with|against)\s+(.+)$', clean_sentence, flags=re.IGNORECASE)
            if std_match:
                val_unit = std_match.group(1).strip()
                rest = std_match.group(3).strip()
                rest = re.sub(r',?\s*as substrates?$', '', rest, flags=re.IGNORECASE)
                results.append({"value": val_unit, "substrate": rest, "condition": condition})
                continue

            # 4. Fallback: No pivot word (e.g., just "1.5 min(-1)")
            val_match = re.match(r'^([~<>\s]*[\d\.,]+(?:[eE][\+\-]?\d+)?(?:x10\^\d+)?)\s*(.*)$', clean_sentence)
            if val_match:
                val = val_match.group(1) + (" " + val_match.group(2) if val_match.group(2) else "")
                results.append({"value": val.strip(), "substrate": "", "condition": condition})
            else:
                results.append({"value": sentence, "substrate": "", "condition": condition})
                
        return results

    def _determine_primary_class(self, is_inactive, high_prob, med_prob, is_passive_transport, is_putative, has_experimental):
        """Determines the functional identity label independent of confidence."""
        if is_inactive:
            return "Pseudo_Enzyme" if (high_prob or med_prob) else "Non_Enzyme"
        
        # Transport reaction (substrates = products)
        if is_passive_transport and not high_prob:
            return "Structural_or_Transport"

        if high_prob:
            if is_putative and not has_experimental:
                return "Inferred_Enzyme"
            return "Enzyme"
        
        elif med_prob: 
            return "Suspected_Enzyme"
        
        return "Non_Enzyme"

    def _calculate_confidence(self, features, entry_data):
        score = 0
        
        # --- PRIMARY ENZYMATIC EVIDENCE (High Weight) ---
        
        # 1. Evaluate Exact Catalytic Activity Evidence
        cat_act_ev = features.get('exp_level_catalytic_activity', 'NONE')
        if cat_act_ev == "STRICT_EXP": score += 60
        elif cat_act_ev == "HTP_EXP": score += 40
        elif cat_act_ev == "CURATED": score += 20
        elif cat_act_ev == "AUTOMATED": score += 5

        # 2. Evaluate Exact EC Number Evidence
        ec_ev = features.get('exp_level_ec', 'NONE')
        if ec_ev == "STRICT_EXP": score += 50
        elif ec_ev == "HTP_EXP": score += 35
        elif ec_ev == "CURATED": score += 15
        elif ec_ev == "AUTOMATED": score += 5
        
        # 3. Evaluate Exact Catalytic GO Term Evidence (NEW)
        go_ev = features.get('exp_level_go', 'NONE')
        if go_ev == "STRICT_EXP": score += 40
        elif go_ev == "HTP_EXP": score += 25
        elif go_ev == "CURATED": score += 15
        elif go_ev == "AUTOMATED": score += 5

        # 4. Structural & Strict Database Support
        if any(db in ["BRENDA", "SABIO-RK"] for db in entry_data['lists'].get('enzyme_dbs', [])): 
            score += 30
        if bool(entry_data['lists'].get('rhea_ids')): 
            score += 20
        if features.get('has_kinetics', False): 
            score += 25
        if features.get('has_active_site', False): 
            score += 15

        # --- HOLISTIC PROTEIN QUALITY BOOSTERS (Low Weight) ---
        pe = entry_data.get('meta', {}).get('protein_existence', '')
        if pe == '1': score += 10
        elif pe == '2': score += 5
        
        if features.get('has_experimental_evidence', False): score += 10
        if features.get('has_catalytic_kw', False): score += 10

        # --- TIER ASSIGNMENT ---
        tier = "Unrated"
        if score >= 100: tier = "Gold"      
        elif score >= 50: tier = "Silver" 
        elif score >= 20: tier = "Bronze"    

        return tier, score

    # def _determine_confidence_tier(self, primary_class, features, entry_data):
    #     """
    #     Evaluates the strength of evidence for the classification.
    #     Focuses on enzymatic evidence specifically.
    #     """
    #     # 1. Gold Tier: Direct Experimental Evidence or High-Quality Curation
    #     # Note: features['has_experimental_evidence'] should correspond to ECO:0000269 (Manual/Exp)
    #     strict_dbs = ["BRENDA", "SABIO-RK", "Rhea", "KEGG"]
    #     has_strict_db = any(db in entry_data['lists'].get('enzyme_dbs', []) for db in strict_dbs)
    #     has_rhea = bool(entry_data['lists'].get('rhea_ids'))
        
    #     if (features.get('has_experimental_evidence') or has_strict_db or has_rhea) and primary_class == "Enzyme":
    #         return "Gold"

    #     # 2. Silver Tier: Strongly Annotated but not explicitly experimental
    #     # Includes proteins with EC numbers or full reaction text but lacks the ECO code/Strict DB cross-ref
    #     if primary_class == "Enzyme":
    #         return "Silver"

    #     # 3. Bronze Tier: Homology, Rule-based, or Suspected
    #     if primary_class in ["Inferred_Enzyme", "Suspected_Enzyme"]:
    #         return "Bronze"

    #     # Default for Non_Enzyme, Pseudo_Enzyme, or Structural
    #     return "Low_or_None"

    def _evaluate_go_evidence(self, go_text_list):
        """
        Parses GO annotations and returns the highest evidence tier 
        specifically for valid catalytic GO terms.
        """
        highest_go_ev = "NONE"
        ranks = {"STRICT_EXP": 4, "HTP_EXP": 3, "CURATED": 2, "AUTOMATED": 1, "NONE": 0}
        
        for go_item in go_text_list:
            parts = go_item.split('|')
            if len(parts) >= 2:
                g_id, eco = parts[0], parts[1]
                if g_id in self.valid_catalytic_go_ids:
                    level = "NONE"
                    if eco in ECO_STRICT_EXP: level = "STRICT_EXP"
                    elif eco in ECO_HTP_EXP: level = "HTP_EXP"
                    elif eco in ECO_CURATED: level = "CURATED"
                    elif eco in ECO_AUTOMATED: level = "AUTOMATED"
                    
                    if ranks[level] > ranks[highest_go_ev]:
                        highest_go_ev = level
                        
        return highest_go_ev

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
        of exact boolean flags/metrics for downstream mapping.
        """

        # 1. Base Extractions & Setup
        caution_text = entry_data['text'].get('caution', '').lower()
        is_inactive_text = features.get('is_inactive_text', False) or "inactive" in caution_text or "catalytically dead" in caution_text
        
        protein_name = entry_data['text'].get('protein_name', '')
        has_enzyme_name = bool(ENZYME_NAME_ENDINGS.search(protein_name))
        
        # Putative terms
        putative_terms = ['putative', 'probable', 'uncharacterized', 'predicted']
        is_putative = any(term in protein_name.lower() for term in putative_terms)
        has_experimental = features.get('has_experimental_evidence', False)
        
        # Protein existence
        protein_existence = entry_data.get('meta', {}).get('protein_existence', '')
        
        go_ids = entry_data['lists'].get('go_ids', [])
        catalytic_gos = [go for go in go_ids if go in self.valid_catalytic_go_ids] if self.valid_catalytic_go_ids else []

        go_text_list = entry_data['lists'].get('go_text', [])
        highest_go_ev = self._evaluate_go_evidence(go_text_list)
        features['exp_level_go'] = highest_go_ev

        enz_dbs = entry_data['lists'].get('enzyme_dbs', [])
        strict_enz_db = any(db in ["BRENDA", "SABIO-RK"] for db in enz_dbs)

        # Parse Kinetics On-the-Fly
        parsed_kinetics = {'KM': [], 'Vmax': [], 'kcat': []}
        
        for p_type, raw_string in [('KM', entry_data['text'].get('kinetics_km', '')), 
                                   ('Vmax', entry_data['text'].get('kinetics_vmax', '')), 
                                   ('kcat', entry_data['text'].get('kinetics_kcat', ''))]:
            if raw_string:
                for raw_val in raw_string.split(" | "):
                    if raw_val.strip():
                        parsed_kinetics[p_type].extend(self._parse_kinetic_string(raw_val.strip(), p_type))

        # 2. Probability Tiering for Primary Class Label
        high_prob = any([
            bool(entry_data['lists'].get('ec_numbers')),
            (features.get('has_reaction') and not features.get('is_passive_transport', False)),
            strict_enz_db,
            (bool(entry_data['lists'].get('rhea_ids')) and not features.get('is_passive_transport', False)),
            (protein_existence == '1' and features.get('has_active_site', False))
        ])

        med_prob = any([
            features.get('has_active_site', False),
            features.get('has_catalytic_kw', False),
            features.get('has_biosynth_text', False),
            has_enzyme_name,
            bool(catalytic_gos)
        ])
        

        # 3.1 Determine Class 
        primary_class = self._determine_primary_class(
            is_inactive=is_inactive_text, 
            high_prob=high_prob,
            med_prob=med_prob,
            is_passive_transport=features.get('is_passive_transport', False),
            is_putative=is_putative,
            has_experimental=has_experimental
        )

        # 3.2 Calculate Decoupled Confidence Tier and Score
        confidence_tier, confidence_score = self._calculate_confidence(features, entry_data)

        # 3.2 Determine Confidence Tier ---
        #confidence_tier = self._determine_confidence_tier(primary_class, features, entry_data)


        # 4. Gather Evidence Strings
        ev_fields, ev_vals = self._gather_evidence(
            entry_data=entry_data, 
            features=features, 
            catalytic_gos=catalytic_gos, 
            protein_name=protein_name
        )


        # 5. Build Row Metadata 
        row_metadata = {
            "confidence_tier": confidence_tier,
            "confidence_score": confidence_score,
            "raw_evidence_map": json.dumps(entry_data.get('evidence_map', {})) if entry_data.get('evidence_map') else "{}",
            "exp_level_catalytic_activity": features.get("exp_level_catalytic_activity", "NONE"),
            "exp_level_ec": features.get("exp_level_ec", "NONE"),
            "exp_level_go": highest_go_ev,
            "kinetics_km": json.dumps(parsed_kinetics['KM']) if parsed_kinetics['KM'] else None,
            "kinetics_vmax": json.dumps(parsed_kinetics['Vmax']) if parsed_kinetics['Vmax'] else None,
            "kinetics_kcat": json.dumps(parsed_kinetics['kcat']) if parsed_kinetics['kcat'] else None,
            "flag_has_ec": bool(entry_data['lists'].get('ec_numbers')),
            "flag_has_active_site": features.get('has_active_site', False),
            "flag_has_kinetics": features.get('has_kinetics', False),
            "flag_has_reaction": features.get('has_reaction', False),
            "flag_has_experimental": has_experimental,
            "flag_is_pseudo": (primary_class == "Pseudo_Enzyme"),
            "flag_is_passive_transport": features.get('is_passive_transport', False), 
            "evidence_field": "|".join(ev_fields) if ev_fields else "None",
            "evidence_value": "|".join(ev_vals) if ev_vals else "None"
        }

        return primary_class, row_metadata

    def assign_category(self, kw_ids, ec_numbers, primary_class, has_sec_met_text):
        """Determines Functional Category (Context)."""
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
