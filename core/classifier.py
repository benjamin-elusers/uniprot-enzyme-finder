# core/classifier.py

import re
import json
from config.keywords import FUNCTIONAL_CATEGORIES, EC_CLASSES, CATALYTIC_KW_MAP
import config.patterns as pt
from config.evidences import ECO_STRICT_EXP, ECO_HTP_EXP, ECO_CURATED, ECO_AUTOMATED


class EnzymeClassifier:
    def __init__(self, valid_catalytic_go_ids=None, eco_defs=None):
        self.valid_catalytic_go_ids = valid_catalytic_go_ids or set()
        
        # Flatten Context Map for O(1) lookups
        self.context_kw_map = {}
        for category, kws in FUNCTIONAL_CATEGORIES.items():
            for kw in kws:
                self.context_kw_map[kw] = category

        # 2. Build O(1) ECO Tier Lookup Map
        self.eco_tier_map = {}
        for eco in ECO_STRICT_EXP: self.eco_tier_map[eco] = "STRICT_EXP"
        for eco in ECO_HTP_EXP: self.eco_tier_map[eco] = "HTP_EXP"
        for eco in ECO_CURATED: self.eco_tier_map[eco] = "CURATED"
        for eco in ECO_AUTOMATED: self.eco_tier_map[eco] = "AUTOMATED"
    # =========================================================================
    #   EVIDENCE MAPPING & OBO PARSING
    # =========================================================================

    def _rank_eco_codes(self, eco_codes):
        """Core logic to rank a list of raw ECO codes and map the OBO definition."""
        if not eco_codes:
            return "NONE", "None"

        highest_tier = "NONE"
        best_eco = None
        ranks = {"STRICT_EXP": 4, "HTP_EXP": 3, "CURATED": 2, "AUTOMATED": 1, "NONE": 0}

        for eco in eco_codes:
            level = self.eco_tier_map.get(eco, "OTHER_ECO")

            # OTHER_ECO ranks slightly above NONE (0) but below AUTOMATED (1)
            r_level = 0.5 if level == "OTHER_ECO" else ranks.get(level, 0)
            r_highest = 0.5 if highest_tier.startswith("OTHER_ECO") else ranks.get(highest_tier, 0)

            if r_level > r_highest:
                highest_tier = f"OTHER_{eco}" if level == "OTHER_ECO" else level
                best_eco = eco

        if highest_tier == "NONE":
            return "NONE", "None"

        def_str = self.eco_defs.get(best_eco, "No definition found in OBO")
        return highest_tier, f"{best_eco}: {def_str}"


    def _evaluate_keys_evidence(self, keys_list, evidence_map):
        """Extracts ECO codes from the evidence_map using keys and ranks them."""
        if not keys_list:
            return "NONE", "None"
            
        # Extract valid ECO codes using list comprehension
        eco_codes = [
            evidence_map[key]['eco'] 
            for key in keys_list if key in evidence_map and 'eco' in evidence_map[key]
        ]
        
        return self._rank_eco_codes(eco_codes)


    def _evaluate_go_evidence_extended(self, go_text_list):
        """Extracts ECO codes from valid catalytic GO terms and ranks them."""
        eco_codes = []
        for go_item in go_text_list:
            parts = go_item.split('|')
            if len(parts) >= 2:
                g_id, eco = parts[0], parts[1]
                if g_id in self.valid_catalytic_go_ids:
                    eco_codes.append(eco)
                    
        return self._rank_eco_codes(eco_codes)


    def _build_evidence_summary(self, entry_data):
        """Inverts the XML evidence dictionary into a structured summary JSON array."""
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

    # =========================================================================
    #   TEXT MINING & FEATURE PARSING
    # =========================================================================

    def _clean_chemical_names(self, side_string):
        """ Removes transport/localization/compartment tags like (in), (out), or [cytoplasm]. """
        cleaned = pt.TRANSPORT_PATTERN.sub('', side_string)
        cleaned = pt.COMPARTMENT_PATTERN.sub('', cleaned)
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
        results = []
        
        # Pre-cleaning
        text = pt.PUBMED_CLEAN_PATTERN.sub('', text)
        text = re.sub(r'\babout\b', '~', text, flags=re.IGNORECASE)
        text = re.sub(r'\binferior to\b', '<', text, flags=re.IGNORECASE)
        
        sentences = [s.strip() for s in re.split(r'(?<=[a-zA-Z\)])\.\s+|(?<=[a-zA-Z\)])\.$', text) if s.strip()] if param_type == "kcat" else [text.strip()]
        
        # Helper to format results cleanly
        def add_res(val, sub="", cond=""):
            results.append({"value": val.strip(), "substrate": sub.strip(), "condition": cond.strip()})

        for sentence in sentences:
            sl = sentence.lower()
            if param_type == "kcat" and ("kcat" not in sl or "kcat/km" in sl or "efficiency" in sl): 
                continue 
            if any(p in sl for p in ["calculations are based", "fold increase", "orders of magnitude", "remains unchanged"]): 
                continue

            # 1. Complex "Respectively" Syntax
            match = pt.KIN_RESP_1.search(sentence)
            if match:
                v1, v2, s1, s2 = match.groups()
                add_res(v1, s1); add_res(v2, s2)
                continue
                
            match = pt.KIN_RESP_2.search(sentence)
            if match:
                s1, s2, v1, v2 = match.groups()
                add_res(v1, s1); add_res(v2, s2)
                continue
                
            # 2. Multi-Condition Syntax
            match = pt.KIN_MULTI_COND.search(sentence)
            if match:
                c1, v1, c2, v2 = match.groups()
                add_res(v1, "", c1); add_res(v2, "", c2)
                continue

            # 3. Standard Extraction (Condition -> Prefix -> Value/Substrate)
            condition = ""
            cond_match = pt.KIN_COND_EXTRACT.search(sentence)
            if cond_match:
                condition = cond_match.group(0).strip('()')
                sentence = sentence[:cond_match.start()].strip() + " " + sentence[cond_match.end():].strip()

            clean_sentence = pt.KIN_PREFIX_1.sub('', sentence)
            clean_sentence = pt.KIN_PREFIX_2.sub('', clean_sentence).strip()

            match = pt.KIN_INV.search(clean_sentence)
            if match:
                add_res(match.group(2), match.group(1), condition)
                continue

            match = pt.KIN_STD.search(clean_sentence)
            if match:
                add_res(match.group(1), re.sub(r',?\s*as substrates?$', '', match.group(3), flags=re.IGNORECASE), condition)
                continue

            match = pt.KIN_VAL.search(clean_sentence)
            if match:
                val = f"{match.group(1)} {match.group(2)}".strip()
                add_res(val, "", condition)
            else:
                add_res(sentence, "", condition)
                
        return results

    def _gather_evidence(self, entry_data, features, catalytic_gos, protein_name):
        """Gathers evidence using efficient dictionary assignments."""
        lists = entry_data.get('lists', {})
        ev = {}

        if ecs := lists.get('ec_numbers'): ev["EC_Number"] = ";".join(ecs)
        if features.get('has_reaction'): ev["Reaction_Text"] = "Reaction found"
        if dbs := lists.get('enzyme_dbs'): ev["Enzyme_DB"] = ";".join(dbs)
        if rheas := lists.get('rhea_ids'): ev["Rhea_ID"] = ";".join(rheas)
        if catalytic_gos: ev["GO_Catalytic"] = ";".join(catalytic_gos)
        
        if features.get('has_active_site'): ev["Active_Site"] = "Active site annotated"
        if features.get('has_catalytic_kw'): ev["Keyword_Catalytic"] = "Catalytic keyword matched"
        if protein_name and pt.ENZYME_NAME_ENDINGS.search(protein_name): ev["Protein_Name"] = protein_name
        if features.get('has_biosynth_text'): ev["Text_Mining"] = "Biosynthesis context"

        return list(ev.keys()), list(ev.values())

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

    def _determine_primary_class(self, is_inactive, high_prob, med_prob, is_passive_transport, is_putative, has_enzyme_exp):
        """Determines the functional identity label based on probability and experimental validation."""
        if is_inactive: return "Pseudo_Enzyme" if (high_prob or med_prob) else "Non_Enzyme"
        if is_passive_transport and not high_prob: return "Structural_or_Transport"
        
        if high_prob: 
            # It is INFERRED if it explicitly has a putative name, 
            # OR if it completely lacks experimental evidence for its enzymatic features.
            if is_putative or not has_enzyme_exp:
                return "Inferred_Enzyme"
            return "Enzyme"
            
        if med_prob: return "Suspected_Enzyme"
        return "Non_Enzyme"

    def _calculate_confidence(self, features, entry_data, is_inactive, is_putative, is_fragment, pe):
        score = 0
        cat_act_ev = features.get('catalytic_activity_eco_rk', 'NONE')
        if cat_act_ev == "STRICT_EXP": score += 60
        elif cat_act_ev == "HTP_EXP": score += 40
        elif cat_act_ev == "CURATED": score += 30
        elif cat_act_ev == "AUTOMATED": score += 10

        ec_ev = features.get('ec_eco_rk', 'NONE')
        if ec_ev == "STRICT_EXP": score += 50
        elif ec_ev == "HTP_EXP": score += 30
        elif ec_ev == "CURATED": score += 20
        elif ec_ev == "AUTOMATED": score += 10
        
        # 1. Calculate GO Term Points
        go_ev = features.get('go_term_eco_rk', 'NONE')
        go_points = 0
        if go_ev == "STRICT_EXP": go_points = 40
        elif go_ev == "HTP_EXP": go_points = 25
        elif go_ev == "CURATED": go_points = 10
        elif go_ev == "AUTOMATED": go_points = 5

        # 2. Calculate Keyword Points
        kw_ev = features.get('keyword_eco_rk', 'NONE')
        kw_points = 0
        if kw_ev == "STRICT_EXP": kw_points = 40
        elif kw_ev == "HTP_EXP": kw_points = 25
        elif kw_ev == "CURATED": kw_points = 10
        elif kw_ev == "AUTOMATED": kw_points = 5

        # 3. Take the MAXIMUM of the two highly correlated features
        score += max(go_points, kw_points)

        # Optional: Corroboration Bonus
        # If they BOTH have experimental/curated evidence, give a tiny flat bonus 
        # to reward the density of the annotation without double-counting the main experiment.
        if go_points >= 10 and kw_points >= 10:
            score += 5
        
        if any(db in ["BRENDA", "SABIO-RK"] for db in entry_data['lists'].get('enzyme_dbs', [])): score += 30
        if bool(entry_data['lists'].get('rhea_ids')): score += 20
        if features.get('has_kinetics', False): score += 50
        if features.get('has_active_site', False): score += 20
        if features.get('has_experimental_evidence', False): score += 5

        if pe == '1': score += 10
        elif pe == '2': score += 5

        # PENALTIES
        if is_inactive: score -= 50
        if is_putative: score -= 20
        if pe in ['4', '5']: score -= 15 
        if is_fragment: score -= 10

        # Partial EC Penalty (Averaged for promiscuous/multifunctional enzymes)
        ec_list = entry_data['lists'].get('ec_numbers', [])
        if ec_list:
            total_ec_penalty = 0
            for ec in ec_list:
                missing_digits = ec.count("-")
                if missing_digits == 1:
                    total_ec_penalty += 5    # e.g., 1.1.1.-
                elif missing_digits == 2:
                    total_ec_penalty += 15   # e.g., 1.1.-.-
                elif missing_digits >= 3:
                    total_ec_penalty += 25   # e.g., 1.-.-.-
            
            # Apply the mean penalty across all functions, rounded to nearest integer
            average_penalty = round(total_ec_penalty / len(ec_list))
            score -= average_penalty

        tier = "Unrated"
        if score >= 100: tier = "Gold"      
        elif score >= 50: tier = "Silver" 
        elif score >= 20: tier = "Bronze"    

        return tier, score

    def classify_entry(self, entry_data, features):
        evidence_map = entry_data.get('evidence_map', {})
        evidence_usages = entry_data.get('evidence_usages', {})
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

        features['has_enzyme_name'] = bool(pt.ENZYME_NAME_ENDINGS.search(protein_name))
        features['is_inactive_text'] = bool(pt.HARD_INACTIVE_PATTERN.search(protein_name)) or "inactive" in caution_text or "catalytically dead" in caution_text
        features['is_putative'] = bool(pt.PUTATIVE_TERMS.search(protein_name))
        features['has_sec_met_text'] = bool(pt.SEC_MET_PATTERNS.search(combined_text))
        features['has_biosynth_text'] = bool(pt.BIOSYNTH_PATTERN.search(combined_text))
        features['caution_activity_relevance'] = any(kw in caution_text for kw in ["affinity", "specificity", "activity", "kinetics", "catalysis", "dead", "inactive"])

        kw_ids = entry_data['lists'].get('kw_ids', [])
        features['has_catalytic_kw'] = any(kw in CATALYTIC_KW_MAP for kw in kw_ids)

        if pt.EC_NAME_PATTERN.search(combined_text):
            found_ecs = pt.EC_NAME_PATTERN.findall(combined_text)
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

        # --- DYNAMIC EVIDENCE TRACKING ---
        def get_keys_for_usage(usage_pattern):
            return [k for k, usages in evidence_usages.items() if usage_pattern.lower() in usages.lower()]
            
        feature_key_mappings = {
            'ec': entry_data['text'].get('name_ec_evidence_keys', '').split() + entry_data['text'].get('db_ec_evidence_keys', '').split(),
            'catalytic_activity': entry_data['text'].get('reaction_evidence_keys', '').split(),
            'active_site': get_keys_for_usage("active site"),
            'cofactor': get_keys_for_usage("cofactor"),
            'rhea_id': get_keys_for_usage("rhea"),
            'keyword': get_keys_for_usage("keyword"),
            'kinetics': get_keys_for_usage("kinetics")
        }
        
        for feat_name, keys in feature_key_mappings.items():
            features[f'{feat_name}_eco_rk'], features[f'{feat_name}_eco_id'] = self._evaluate_keys_evidence(keys, evidence_map)
            
        features['go_term_eco_rk'], features['go_term_eco_id'] = self._evaluate_go_evidence_extended(go_text_list)

        features['has_experimental_evidence'] = any(ev['eco'] in ECO_STRICT_EXP or ev['eco'] in ECO_HTP_EXP for ev in evidence_map.values())
        enzymatic_ranks = [
            features.get('ec_eco_rk'), 
            features.get('catalytic_activity_eco_rk'), 
            features.get('go_term_eco_rk'), 
            features.get('active_site_eco_rk'),
            features.get('keyword_eco_rk')
        ]
        features['has_enzyme_experimental_evidence'] = any(rk in ["STRICT_EXP", "HTP_EXP"] for rk in enzymatic_ranks)


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
            bool(pt.ENZYME_NAME_ENDINGS.search(protein_name)),
            bool(catalytic_gos)
        ])

        primary_class = self._determine_primary_class(
            is_inactive=features['is_inactive_text'], 
            high_prob=high_prob,
            med_prob=med_prob,
            is_passive_transport=features['is_passive_transport'],
            is_putative=features['is_putative'],
            has_enzyme_exp=features['has_enzyme_experimental_evidence']
        )

        confidence_tier, confidence_score = self._calculate_confidence(
            features, entry_data, features['is_inactive_text'], features['is_putative'], is_fragment, protein_existence
        )

        ev_fields, ev_vals = self._gather_evidence(entry_data, features, catalytic_gos, protein_name)

        # 6. Build Row Metadata
        row_metadata = {
            "confidence_tier": confidence_tier,
            "confidence_score": confidence_score,
            "raw_evidence_map": json.dumps(evidence_map) if evidence_map else "{}",
            "detailed_evidence_map": self._build_evidence_summary(entry_data),
            
            "ec_eco_rk": features.get('ec_eco_rk', 'NONE'),
            "ec_eco_id": features.get('ec_eco_id', 'None'),
            "catalytic_activity_eco_rk": features.get('catalytic_activity_eco_rk', 'NONE'),
            "catalytic_activity_eco_id": features.get('catalytic_activity_eco_id', 'None'),
            "go_term_eco_rk": features.get('go_term_eco_rk', 'NONE'),
            "go_term_eco_id": features.get('go_term_eco_id', 'None'),
            "active_site_eco_rk": features.get('active_site_eco_rk', 'NONE'),
            "active_site_eco_id": features.get('active_site_eco_id', 'None'),
            "cofactor_eco_rk": features.get('cofactor_eco_rk', 'NONE'),
            "cofactor_eco_id": features.get('cofactor_eco_id', 'None'),
            "keyword_eco_rk": features.get('keyword_eco_rk', 'NONE'),
            "keyword_eco_id": features.get('keyword_eco_id', 'None'),
            "rhea_id_eco_rk": features.get('rhea_id_eco_rk', 'NONE'),
            "rhea_id_eco_id": features.get('rhea_id_eco_id', 'None'),
            "kinetics_eco_rk": features.get('kinetics_eco_rk', 'NONE'),
            "kinetics_eco_id": features.get('kinetics_eco_id', 'None'),
            
            "exp_level_catalytic_activity": features.get("catalytic_activity_eco_rk", "NONE"),
            "exp_level_ec": features.get("ec_eco_rk", "NONE"),
            "exp_level_go": features.get("go_term_eco_rk", "NONE"),
            
            "kinetics_km": json.dumps(parsed_kinetics['KM']) if parsed_kinetics['KM'] else None,
            "kinetics_vmax": json.dumps(parsed_kinetics['Vmax']) if parsed_kinetics['Vmax'] else None,
            "kinetics_kcat": json.dumps(parsed_kinetics['kcat']) if parsed_kinetics['kcat'] else None,
            "flag_has_ec": bool(ec_num),
            "flag_is_multifunctional": len(ec_list) > 1,
            "flag_is_promiscuous": any("-" in ec for ec in ec_list),
            "flag_has_active_site": features['has_active_site'],
            "flag_has_kinetics": features['has_kinetics'],
            "flag_has_reaction": features['has_reaction'],
            "flag_has_experimental": features['has_experimental_evidence'],
            "has_enzyme_experimental_evidence": features['has_enzyme_experimental_evidence'],
            "flag_is_pseudo": (primary_class == "Pseudo_Enzyme"),
            "flag_is_passive_transport": features['is_passive_transport'], 
            "evidence_field": "|".join(ev_fields) if ev_fields else "None",
            "evidence_value": "|".join(ev_vals) if ev_vals else "None"
        }

        return primary_class, row_metadata
