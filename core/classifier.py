# core/classifier.py

import re
import json
from config.keywords import FUNCTIONAL_CATEGORIES, EC_CLASSES, CATALYTIC_KW_MAP
from config.patterns import (
    EC_NAME_PATTERN, HARD_INACTIVE_PATTERN, SOFT_INACTIVE_PATTERN, PUTATIVE_TERMS,
    BIOSYNTH_PATTERN, SEC_MET_PATTERNS, ENZYME_NAME_ENDINGS
)
from config.evidences import ECO_STRICT_EXP, ECO_HTP_EXP, ECO_CURATED, ECO_AUTOMATED

class EnzymeClassifier:
    def __init__(self, valid_catalytic_go_ids=None):
        self.valid_catalytic_go_ids = valid_catalytic_go_ids or set()
        
        # Flatten Context Map for O(1) lookups
        self.context_kw_map = {}
        for category, kws in FUNCTIONAL_CATEGORIES.items():
            for kw in kws:
                self.context_kw_map[kw] = category

    # =========================================================================
    #   BIOLOGICAL LOGIC & EVALUATORS 
    # =========================================================================
    
    def _build_evidence_summary(self, entry_data):
        """Inverts the XML evidence dictionary into a structured summary table."""
        evidence_map = entry_data.get('evidence_map', {})
        evidence_usages = entry_data.get('evidence_usages', {}) # <-- Load global usages
        if not evidence_map: return "[]"
            
        summary = []
        parsed_keys = {}
        
        # 1. Map explicitly tracked enzymatic keys from the parser
        for k in entry_data['text'].get('name_ec_evidence_keys', '').split():
            parsed_keys[k] = "EC_from_Name"
        for k in entry_data['text'].get('db_ec_evidence_keys', '').split():
            parsed_keys[k] = "EC_Database_Ref"
        for k in entry_data['text'].get('reaction_evidence_keys', '').split():
            parsed_keys[k] = "Catalytic_Activity"

        # 2. Add GO Terms explicitly
        for go_item in entry_data['lists'].get('go_text', []):
            parts = go_item.split('|')
            if len(parts) >= 2:
                g_id, eco = parts[0], parts[1]
                is_cat = (g_id in self.valid_catalytic_go_ids)
                summary.append({
                    "annotation": "GO_Term",
                    "evidence_key": "N/A",
                    "evidence_id": g_id,
                    "code": eco,
                    "scopes": [],
                    "is_parsed": is_cat
                })

        # 3. Process the dictionary with global fallback
        for ev_key, ev_data in evidence_map.items():
            if ev_key in parsed_keys:
                annotation = parsed_keys[ev_key]
                is_parsed = True
            else:
                # Dynamically grab the XML tag/type where this key was used!
                annotation = evidence_usages.get(ev_key, "Unknown_Feature")
                is_parsed = False
            
            summary.append({
                "annotation": annotation,
                "evidence_key": ev_key,
                "evidence_id": ev_data.get('pubmed', 'N/A'),
                "code": ev_data.get('eco', ''),
                "scopes": ev_data.get('scopes', []),
                "is_parsed": is_parsed
            })
            
        return json.dumps(summary)

    def _evaluate_evidence_level(self, evidence_str, evidence_map):
        """Translates raw evidence keys into ECO Tiers."""
        if not evidence_str: return "NONE"
            
        keys = evidence_str.split()
        found_levels = set()
        other_ecos = [] 
        
        for key in keys:
            ev_data = evidence_map.get(key)
            if not ev_data: continue
            
            eco = ev_data['eco']
            if eco in ECO_STRICT_EXP: return "STRICT_EXP" 
            elif eco in ECO_HTP_EXP: found_levels.add("HTP_EXP")
            elif eco in ECO_CURATED: found_levels.add("CURATED")
            elif eco in ECO_AUTOMATED: found_levels.add("AUTOMATED")
            else:
                other_ecos.append(eco) 

        if "HTP_EXP" in found_levels: return "HTP_EXP"
        if "CURATED" in found_levels: return "CURATED"
        if "AUTOMATED" in found_levels: return "AUTOMATED"
        if other_ecos: return f"OTHER_{other_ecos[0]}"
        return "NONE"

    def _max_evidence_level(self, level1, level2):
        ranks = {"STRICT_EXP": 4, "HTP_EXP": 3, "CURATED": 2, "AUTOMATED": 1, "NONE": 0}

        r1 = ranks.get(level1, 0)
        r2 = ranks.get(level2, 0)

        if r1 > r2: return level1
        elif r2 > r1: return level2
        else:
            if level1 != "NONE": return level1
            return level2

    def _clean_chemical_names(self, side_string):
        """ Removes transport/localization/compartment tags like (in), (out), or [cytoplasm]. """
        cleaned = re.sub(r'\((in|out|external|internal)\)', '', side_string, flags=re.IGNORECASE)
        cleaned = re.sub(r'\[.*?\]', '', cleaned)
        molecules = [m.strip().lower() for m in cleaned.split(' + ')]
        return " + ".join(sorted(filter(None, molecules)))

    def _process_reactions(self, catalytic_text):
        """Extracts molecules and identifies passive transport."""
        sub_rxns = [r.strip() for r in catalytic_text.split('|')]
        res = {'clean_equations': [], 'substrates': [], 'products': [], 'is_passive_transport': True}
        
        for rxn in sub_rxns:
            if " = " not in rxn: continue
            left, right = rxn.split(" = ", 1)
            
            if self._clean_chemical_names(left) != self._clean_chemical_names(right):
                res['is_passive_transport'] = False
            
            res['clean_equations'].append((left, right))
            res['substrates'].extend([s.strip() for s in left.split(' + ')])
            res['products'].extend([p.strip() for p in right.split(' + ')])

        return res

    def _parse_kinetic_string(self, text, param_type):
        """Structurally parses kinetic strings using linguistic pivots and edge-case handlers."""
        results = []
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
            if any(phrase in sl for phrase in ignore_phrases): continue

            resp_match = re.search(r"(?:(?:kcat|KM|Vmax)\s+(?:is|are|=)\s*)?(.*?)\s+and\s+(.*?)\s+(?:with|for|toward|against)\s+(.*?)\s+and\s+(.*?)\s+as substrates(?:,\s*respectively)?", sentence, flags=re.IGNORECASE)
            if resp_match:
                v1, v2, s1, s2 = resp_match.groups()
                results.append({"value": v1.strip(), "substrate": s1.strip(), "condition": ""})
                results.append({"value": v2.strip(), "substrate": s2.strip(), "condition": ""})
                continue
                
            resp_match_2 = re.search(r"(?:The\s+apparent\s+)?(?:kcat|KM|Vmax)\s+of\s+(.*?)\s+and\s+(.*?)(?:\s+degradation)?\s+(?:is|are|=)\s*(.*?)\s*,?\s*and\s*(.*?)\s*,?\s*respectively", sentence, flags=re.IGNORECASE)
            if resp_match_2:
                s1, s2, v1, v2 = resp_match_2.groups()
                results.append({"value": v1.strip(), "substrate": s1.strip(), "condition": ""})
                results.append({"value": v2.strip(), "substrate": s2.strip(), "condition": ""})
                continue
                
            multi_cond = re.search(r'^(?:kcat|KM|Vmax)\s+(at\s+pH\s+[\d\.]+)\s+(?:is|are|=)\s*([~<>\s]*[\d\.]+\s*[a-zA-Z\/\-\(\)\d]+)\s*,\s*(at\s+pH\s+[\d\.]+)\s+(?:is|are|=)\s*([~<>\s]*[\d\.]+\s*[a-zA-Z\/\-\(\)\d]+)', sentence, flags=re.IGNORECASE)
            if multi_cond:
                c1, v1, c2, v2 = multi_cond.groups()
                results.append({"value": v1.strip(), "substrate": "", "condition": c1.strip()})
                results.append({"value": v2.strip(), "substrate": "", "condition": c2.strip()})
                continue

            condition = ""
            cond_match = re.search(r'\(\s*(in the presence of|at |pH |with |when ).+?\)$', sentence, flags=re.IGNORECASE)
            if cond_match:
                condition = cond_match.group(0).strip('()')
                sentence = sentence[:cond_match.start()].strip() + " " + sentence[cond_match.end():].strip()

            clean_sentence = re.sub(r'^(?:The\s+)?(?:apparent\s+)?(?:kcat|KM|Vmax)\s+(?:is|are|=|of|:)\s+', '', sentence, flags=re.IGNORECASE)
            clean_sentence = re.sub(r'^(?:The\s+)?(?:apparent\s+)?(?:kcat|KM|Vmax)\s+', '', clean_sentence, flags=re.IGNORECASE).strip()

            inv_match = re.search(r'^(?:for|toward|with|against)\s+(.+?)\s+(?:is|are|=)\s*(.+)$', clean_sentence, flags=re.IGNORECASE)
            if inv_match:
                results.append({"value": inv_match.group(2).strip(), "substrate": inv_match.group(1).strip(), "condition": condition})
                continue

            std_match = re.search(r'^(.+?)\s+(for|toward|with|against)\s+(.+)$', clean_sentence, flags=re.IGNORECASE)
            if std_match:
                val_unit = std_match.group(1).strip()
                rest = std_match.group(3).strip()
                rest = re.sub(r',?\s*as substrates?$', '', rest, flags=re.IGNORECASE)
                results.append({"value": val_unit, "substrate": rest, "condition": condition})
                continue

            val_match = re.match(r'^([~<>\s]*[\d\.,]+(?:[eE][\+\-]?\d+)?(?:x10\^\d+)?)\s*(.*)$', clean_sentence)
            if val_match:
                val = val_match.group(1) + (" " + val_match.group(2) if val_match.group(2) else "")
                results.append({"value": val.strip(), "substrate": "", "condition": condition})
            else:
                results.append({"value": sentence, "substrate": "", "condition": condition})
                
        return results

    def _gather_evidence(self, entry_data, features, catalytic_gos, protein_name):
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

    def assign_category(self, kw_ids, ec_numbers, primary_class, has_sec_met_text):
        if has_sec_met_text: return "Secondary_Metabolism"
        for kw in kw_ids:
            if self.context_kw_map.get(kw) == "Secondary_Metabolism": return "Secondary_Metabolism"
        for ec in ec_numbers:
            if ec.startswith("2.7.7") or ec.startswith("6.1.1"): return "Info_Processing"
            if ec.startswith("3.1.2") and len(ec) > 5 and ec[4] in '6789': return "Info_Processing"
        for kw in kw_ids:
            if kw in self.context_kw_map and self.context_kw_map[kw] == "Info_Processing": return "Info_Processing"
        for ec in ec_numbers:
            if ec.startswith("3.1.3") or (ec.startswith("2.7") and not ec.startswith("2.7.7")): return "Signaling"
        for kw in kw_ids:
            if self.context_kw_map.get(kw) == "Signaling": return "Signaling"
        for kw in kw_ids:
            if self.context_kw_map.get(kw) == "Core_Metabolism": return "Core_Metabolism"
        return "Unclassified"

    def get_ec_hierarchy(self, ec_numbers):
        """Safely extracts EC family data and accurately flags partial ECs."""
        if not ec_numbers: return "Non-Enzyme", "", "None"
        primary_ec = ec_numbers[0].strip()
        
        # 1. Evaluate Completeness BEFORE stripping dashes
        completeness = "None"
        if "-" not in primary_ec and primary_ec.count(".") == 3: completeness = "Full"
        elif primary_ec.count(".") >= 1: completeness = "Partial"
        
        # 2. Clean and extract hierarchy
        clean_ec = primary_ec.replace("-", "")
        parts = clean_ec.split(".")
        
        if not parts or not parts[0].isdigit(): return "Unclassified", "", completeness
        
        family_name = EC_CLASSES.get(parts[0], "Unclassified")
        return family_name, parts[0], completeness

    # =========================================================================
    #   SCORING AND CLASSIFICATION
    # =========================================================================

    def _determine_primary_class(self, is_inactive, high_prob, med_prob, is_passive_transport, is_putative, has_catalytic_experiment):
        if is_inactive: return "Pseudo_Enzyme" if (high_prob or med_prob) else "Non_Enzyme"
        if is_passive_transport and not high_prob: return "Structural_or_Transport"
        if high_prob: return "Inferred_Enzyme" if (is_putative and not has_catalytic_experiment) else "Enzyme"
        if med_prob: return "Suspected_Enzyme"
        return "Non_Enzyme"

    def _evaluate_go_evidence(self, go_text_list):
        highest_go_ev = "NONE"
        highest_other = None
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
                    else:
                        highest_other = f"OTHER_{eco}"
                        continue
                        
                    if ranks[level] > ranks.get(highest_go_ev, 0): 
                        highest_go_ev = level
                        
        if highest_go_ev == "NONE" and highest_other:
            return highest_other
            
        return highest_go_ev


    def _calculate_confidence(self, features, entry_data, is_inactive, is_putative, is_fragment, pe):
        score = 0
        cat_act_ev = features.get('exp_level_catalytic_activity', 'NONE')
        if cat_act_ev == "STRICT_EXP": score += 60
        elif cat_act_ev == "HTP_EXP": score += 40
        elif cat_act_ev == "CURATED": score += 20
        elif cat_act_ev == "AUTOMATED": score += 5

        ec_ev = features.get('exp_level_ec', 'NONE')
        if ec_ev == "STRICT_EXP": score += 50
        elif ec_ev == "HTP_EXP": score += 35
        elif ec_ev == "CURATED": score += 15
        elif ec_ev == "AUTOMATED": score += 5
        
        go_ev = features.get('exp_level_go', 'NONE')
        if go_ev == "STRICT_EXP": score += 40
        elif go_ev == "HTP_EXP": score += 25
        elif go_ev == "CURATED": score += 15
        elif go_ev == "AUTOMATED": score += 5

        if any(db in ["BRENDA", "SABIO-RK"] for db in entry_data['lists'].get('enzyme_dbs', [])): score += 30
        if bool(entry_data['lists'].get('rhea_ids')): score += 30
        if features.get('has_kinetics', False): score += 50
        if features.get('has_active_site', False): score += 20
        if features.get('has_experimental_evidence', False): score += 5
        if features.get('has_catalytic_kw', False): score += 10

        if pe == '1': score += 10
        elif pe == '2': score += 5

        # PENALTIES
        if is_inactive: score -= 50
        if is_putative: score -= 20
        if pe in ['4', '5']: score -= 15 
        if is_fragment: score -= 10

        tier = "Unrated"
        if score >= 100: tier = "Gold"      
        elif score >= 50: tier = "Silver" 
        elif score >= 20: tier = "Bronze"    

        return tier, score

    def classify_entry(self, entry_data, features):
        evidence_map = entry_data.get('evidence_map', {})
        protein_name = entry_data['text'].get('protein_name', '')
        protein_existence = entry_data.get('meta', {}).get('protein_existence', '')
        is_fragment = entry_data.get('meta', {}).get('is_fragment', False)
        caution_text = entry_data['text'].get('caution', '').lower()
        combined_text = f"{entry_data['text'].get('pathway', '')} {entry_data['text'].get('function', '')} {entry_data['text'].get('catalytic_activity', '')} {protein_name}"
        ec_num = entry_data['lists'].get('ec_numbers', [])
        go_ids = entry_data['lists'].get('go_ids', [])
        go_text_list = entry_data['lists'].get('go_text', [])
        enz_dbs = entry_data['lists'].get('enzyme_dbs', [])
        rhea_ids = entry_data['lists'].get('rhea_ids', [])

        features['has_enzyme_name'] = bool(ENZYME_NAME_ENDINGS.search(protein_name))
        features['is_inactive_text'] = bool(HARD_INACTIVE_PATTERN.search(protein_name)) or "inactive" in caution_text or "catalytically dead" in caution_text
        features['is_putative'] = bool(PUTATIVE_TERMS.search(protein_name))
        features['has_sec_met_text'] = bool(SEC_MET_PATTERNS.search(combined_text))
        features['has_biosynth_text'] = bool(BIOSYNTH_PATTERN.search(combined_text))
        features['caution_activity_relevance'] = any(kw in caution_text for kw in ["affinity", "specificity", "activity", "kinetics", "catalysis", "dead", "inactive"])

        kw_ids = entry_data['lists'].get('kw_ids', [])
        features['has_catalytic_kw'] = any(kw in CATALYTIC_KW_MAP for kw in kw_ids)

        if EC_NAME_PATTERN.search(combined_text):
            found_ecs = EC_NAME_PATTERN.findall(combined_text)
            entry_data['lists']['ec_from_name'].extend([x.strip("()EC ") for x in found_ecs])
            entry_data['lists']['ec_numbers'] = list(set(entry_data['lists'].get('ec_db', []) + entry_data['lists'].get('ec_from_name', [])))

        catalytic_gos = [go for go in go_ids if go in self.valid_catalytic_go_ids] if self.valid_catalytic_go_ids else []

        rxn_text = entry_data['text'].get('catalytic_activity', '')
        features['has_reaction'] = bool(rxn_text)
        features['is_passive_transport'] = False
        if rxn_text:
            rxn_data = self._process_reactions(rxn_text)
            features['is_passive_transport'] = rxn_data['is_passive_transport']
            entry_data['lists']['substrates'] = list(set(rxn_data['substrates']))
            entry_data['lists']['products'] = list(set(rxn_data['products']))

        parsed_kinetics = {'KM': [], 'Vmax': [], 'kcat': []}
        for p_type, raw_string in [('KM', entry_data['text'].get('kinetics_km', '')), 
                                   ('Vmax', entry_data['text'].get('kinetics_vmax', '')), 
                                   ('kcat', entry_data['text'].get('kinetics_kcat', ''))]:
            if raw_string:
                for raw_val in raw_string.split(" | "):
                    if raw_val.strip():
                        parsed_kinetics[p_type].extend(self._parse_kinetic_string(raw_val.strip(), p_type))

        name_ec_ev = self._evaluate_evidence_level(entry_data['text'].get('name_ec_evidence_keys', ''), evidence_map)
        db_ec_ev = self._evaluate_evidence_level(entry_data['text'].get('db_ec_evidence_keys', ''), evidence_map)
        features['exp_level_ec'] = self._max_evidence_level(name_ec_ev, db_ec_ev)
        features['exp_level_catalytic_activity'] = self._evaluate_evidence_level(entry_data['text'].get('reaction_evidence_keys', ''), evidence_map)
        features['exp_level_go'] = self._evaluate_go_evidence(go_text_list)

        features['has_experimental_evidence'] = any(ev['eco'] in ECO_STRICT_EXP or ev['eco'] in ECO_HTP_EXP for ev in evidence_map.values())
        features['has_experimental_catalytic'] = (features.get("exp_level_catalytic_activity") == "STRICT_EXP" or features.get("exp_level_ec") == "STRICT_EXP")

        features['has_active_site'] = bool(entry_data['lists'].get('active_sites'))
        features['has_kinetics'] = bool(entry_data['text'].get('kinetics_km') or entry_data['text'].get('kinetics_vmax') or entry_data['text'].get('kinetics_kcat'))

        strict_enz_db = any(db in ["BRENDA", "SABIO-RK"] for db in enz_dbs)

        high_prob = any([
            bool(entry_data['lists'].get('ec_numbers')),
            (features['has_reaction'] and not features['is_passive_transport']),
            strict_enz_db,
            (bool(rhea_ids) and not features['is_passive_transport']),
            (protein_existence == '1' and features['has_active_site'])
        ])

        med_prob = any([
            features['has_active_site'],
            features['has_catalytic_kw'],
            features['has_biosynth_text'],
            bool(ENZYME_NAME_ENDINGS.search(protein_name)),
            bool(catalytic_gos)
        ])

        primary_class = self._determine_primary_class(
            is_inactive=features['is_inactive_text'], 
            high_prob=high_prob,
            med_prob=med_prob,
            is_passive_transport=features['is_passive_transport'],
            is_putative=features['is_putative'],
            has_catalytic_experiment=features['has_experimental_catalytic']
        )

        confidence_tier, confidence_score = self._calculate_confidence(
            features, entry_data, features['is_inactive_text'], features['is_putative'], is_fragment, protein_existence
        )

        ev_fields, ev_vals = self._gather_evidence(entry_data, features, catalytic_gos, protein_name)

        # 6. Build Row Metadata (Syntax Fixes Applied)
        row_metadata = {
            "confidence_tier": confidence_tier,
            "confidence_score": confidence_score,
            "raw_evidence_map": json.dumps(evidence_map) if evidence_map else "{}",
            "detailed_evidence_map": self._build_evidence_summary(entry_data),
            "exp_level_catalytic_activity": features.get("exp_level_catalytic_activity", "NONE"),
            "exp_level_ec": features.get("exp_level_ec", "NONE"),
            "exp_level_go": features.get("exp_level_go", "NONE"),
            "kinetics_km": json.dumps(parsed_kinetics['KM']) if parsed_kinetics['KM'] else None,
            "kinetics_vmax": json.dumps(parsed_kinetics['Vmax']) if parsed_kinetics['Vmax'] else None,
            "kinetics_kcat": json.dumps(parsed_kinetics['kcat']) if parsed_kinetics['kcat'] else None,
            "flag_has_ec": bool(ec_num),
            "flag_has_active_site": features['has_active_site'],
            "flag_has_kinetics": features['has_kinetics'],
            "flag_has_reaction": features['has_reaction'],
            "flag_has_experimental": features['has_experimental_evidence'],
            "flag_has_catalytic_experimental": features['has_experimental_catalytic'],
            "flag_is_pseudo": (primary_class == "Pseudo_Enzyme"),
            "flag_is_passive_transport": features['is_passive_transport'], 
            "evidence_field": "|".join(ev_fields) if ev_fields else "None",
            "evidence_value": "|".join(ev_vals) if ev_vals else "None"
        }

        return primary_class, row_metadata