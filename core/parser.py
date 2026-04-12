# core/parser.py

import re
from lxml import etree
from config.patterns import (
    EC_NAME_PATTERN, HARD_INACTIVE_PATTERN, SOFT_INACTIVE_PATTERN,
    BIOSYNTH_PATTERN, SEC_MET_PATTERNS
)
from config.evidences import ECO_STRICT_EXP, ECO_HTP_EXP, ECO_CURATED, ECO_AUTOMATED
import logging

class UniprotEntryParser:
    """
    Parses a single UniProt XML entry into a comprehensive dictionary of raw features, 
    designed to capture all enzymatic, contextual, and structural evidence (v20).
    """
    def __init__(self, xml_str):
        self.parser = etree.XMLParser(recover=True)
        try:
            self.elem = etree.fromstring(xml_str.encode('utf-8'), parser=self.parser)
        except etree.XMLSyntaxError:
            if 'xmlns' not in xml_str:
                xml_str = xml_str.replace('<entry', '<entry xmlns="http://uniprot.org/uniprot"', 1)
                self.elem = etree.fromstring(xml_str.encode('utf-8'), parser=self.parser)
        
        m = re.match(r'\{(.*)\}entry', self.elem.tag)
        self.ns = f"{{{m.group(1)}}}" if m else ""
        
        self.pubmed_scopes = {}
        self.evidence_map = {}

    def _build_evidence_dictionaries(self):
        for ref in self.elem.findall(f"{self.ns}reference"):
            pubmed_id = None
            for dbref in ref.findall(f"{self.ns}citation/{self.ns}dbReference"):
                if dbref.get("type") == "PubMed":
                    pubmed_id = dbref.get("id")
                    break
            
            if pubmed_id:
                scopes = [s.text for s in ref.findall(f"{self.ns}scope") if s.text]
                self.pubmed_scopes[pubmed_id] = scopes

        for ev in self.elem.findall(f"{self.ns}evidence"):
            ev_key = ev.get("key")
            ev_type = ev.get("type") 
            
            pubmed_id = None
            source_dbref = ev.find(f"{self.ns}source/{self.ns}dbReference")
            if source_dbref is not None and source_dbref.get("type") == "PubMed":
                pubmed_id = source_dbref.get("id")

            if ev_key:
                self.evidence_map[ev_key] = {
                    'eco': ev_type,
                    'pubmed': pubmed_id,
                    'scopes': self.pubmed_scopes.get(pubmed_id, []) if pubmed_id else []
                }

    def _evaluate_evidence_level(self, evidence_str):
        if not evidence_str: return "NONE"
            
        keys = evidence_str.split()
        found_levels = set()
        
        for key in keys:
            ev_data = self.evidence_map.get(key)
            if not ev_data: continue
            
            eco = ev_data['eco']
            if eco in ECO_STRICT_EXP: return "STRICT_EXP" 
            elif eco in ECO_HTP_EXP: found_levels.add("HTP_EXP")
            elif eco in ECO_CURATED: found_levels.add("CURATED")
            elif eco in ECO_AUTOMATED: found_levels.add("AUTOMATED")

        if "HTP_EXP" in found_levels: return "HTP_EXP"
        if "CURATED" in found_levels: return "CURATED"
        if "AUTOMATED" in found_levels: return "AUTOMATED"
        return "NONE"

    def _max_evidence_level(self, level1, level2):
        ranks = {"STRICT_EXP": 4, "HTP_EXP": 3, "CURATED": 2, "AUTOMATED": 1, "NONE": 0}
        return level1 if ranks[level1] > ranks[level2] else level2

    def parse(self):
        self._build_evidence_dictionaries()
        
        data = {
            'meta': {},
            'features': {},   # Booleans/Flags
            'lists': {},      # Lists of IDs and Cross-references
            'text': {}        # Extracted text blocks
        }

        # --- 1. Header & Meta Info ---
        data['meta']['accession'] = self._get_text("accession")
        data['meta']['id'] = self._get_text("name")
        data['meta']['reviewed'] = (self.elem.get("dataset") == "Swiss-Prot")
        data['meta']['version'] = self.elem.get("version", "")
        data['meta']['date_modified'] = self.elem.get("modified", "")
        data['meta']['date_created'] = self.elem.get("created", "")
        
        # Protein Existence (1-5 scale)
        pe = self.elem.find(f"{self.ns}proteinExistence")
        pe_level = "1"
        if pe is not None:
            pt = pe.get("type")
            if pt == "evidence at transcript level": pe_level = "2"
            elif pt == "inferred from homology": pe_level = "3"
            elif pt == "predicted": pe_level = "4"
            elif pt == "uncertain": pe_level = "5"
        data['meta']['protein_existence'] = pe_level

        # Sequence Attributes
        seq = self.elem.find(f"{self.ns}sequence")
        data['meta']['length'] = seq.get("length") if seq is not None else "0"
        data['meta']['mass'] = seq.get("mass") if seq is not None else "0"
        data['meta']['is_fragment'] = False
        if seq is not None and 'fragment' in seq.attrib:
            data['meta']['is_fragment'] = (seq.get("fragment") != "single") 

        # --- 2. Taxonomy & Gene ---
        tax_info = self._extract_organism_and_gene()
        data['text'].update(tax_info)
        
        # 3. Name & Inactive Status
        full_name, name_ecs, name_ec_level = self._extract_name()
        data['text']['protein_name'] = full_name
        data['lists']['ec_from_name'] = name_ecs
        
        is_hard_inactive = bool(full_name and HARD_INACTIVE_PATTERN.search(full_name))
        is_soft_inactive = bool(full_name and SOFT_INACTIVE_PATTERN.search(full_name))
        data['features']['is_inactive_text'] = is_hard_inactive or is_soft_inactive

        # --- 4. Comments, Reactions & Kinetics ---
        comments = self._extract_comments()
        data['text'].update(comments['text'])
        data['lists'].update(comments['lists'])
        data['features']['has_kinetics'] = comments['has_kinetics']
        data['features']['has_reaction'] = bool(comments['text']['catalytic_activity'])
        data['features']['is_passive_transport'] = comments.get('is_passive_transport', False)
        data['features']['exp_level_catalytic_activity'] = comments['exp_level_catalytic_activity']

        # --- 5. Evidence Scanning (DB Refs, GO, Keywords, Experimental) ---
        refs = self._scan_references()
        for k, v in refs.items():
            if k in data['lists']: data['lists'][k].extend(v)
            elif isinstance(v, list): data['lists'][k] = v
            
        data['features']['has_experimental_evidence'] = refs.pop('has_experimental_evidence', False)
        
        # Resolve best EC evidence level between names and DBs
        data['features']['exp_level_ec'] = self._max_evidence_level(name_ec_level, refs.pop('exp_level_db_ec', "NONE"))

        # --- 6. Sequence Features (Active Sites, Mutagenesis) ---
        feats = self._extract_seq_features()
        data['lists']['active_sites'] = feats['active_sites']
        data['lists']['binding_sites'] = feats['binding_sites']
        data['lists']['mutagenesis'] = feats['mutagenesis']
        data['features']['has_active_site'] = bool(feats['active_sites'])
        data['features']['has_binding_site'] = bool(feats['binding_sites'])

        # Check mutagenesis for inactive flags
        if any(kw in m.lower() for m in feats['mutagenesis'] for kw in ["inactive", "loss", "abolish"]):
            data['features']['is_inactive_text'] = True

        # Check caution for inactive relevance
        caution_text = data['text'].get('caution', '').lower()
        data['features']['cauti1on_activity_relevance'] = any(kw in caution_text for kw in ["affinity", "specificity", "activity", "kinetics", "catalysis", "dead", "inactive"])

        # --- 7. Text Mining Flags ---
        combined_text = (
            f"{data['text'].get('pathway', '')} {data['text'].get('function', '')} "
            f"{data['text'].get('catalytic_activity', '')} {full_name} "
        )
        data['features']['has_text_mining'] = bool(SEC_MET_PATTERNS.search(combined_text))
        data['features']['has_biosynth_text'] = bool(BIOSYNTH_PATTERN.search(combined_text))
        
        # Catch stray ECs in text
        if EC_NAME_PATTERN.search(combined_text):
            found_ecs = EC_NAME_PATTERN.findall(combined_text)
            data['lists']['ec_from_name'].extend([x.strip("()EC ") for x in found_ecs])
            
        data['lists']['ec_numbers'] = list(set(data['lists'].get('ec_db', []) + data['lists'].get('ec_from_name', [])))

        return data

    # =========================================================================
    #   EXTRACTION HELPERS
    # =========================================================================

    def _get_text(self, path):
        found = self.elem.find(f"{self.ns}{path}")
        return found.text if found is not None else ""

    def _extract_name(self):
        protein = self.elem.find(f"{self.ns}protein")
        name = ""
        ecs = []
        highest_ec_level = "NONE"
        
        if protein is not None:
            for n_type in ["recommendedName", "submittedName", "alternativeName"]:
                for group in protein.findall(f"{self.ns}{n_type}"):
                    fn = group.find(f"{self.ns}fullName")
                    if fn is not None and not name: name = fn.text
                    # Extract EC numbers embedded in names
                    for ec in group.findall(f"{self.ns}ecNumber"):
                        if ec.text: 
                            ecs.append(ec.text)
                            level = self._evaluate_evidence_level(ec.get("evidence", ""))
                            highest_ec_level = self._max_evidence_level(highest_ec_level, level)
                                
        return name, ecs, highest_ec_level

    def _extract_organism_and_gene(self):
        out = {'organism_name': '', 'organism_id': '', 'lineage': '', 'kingdom': 'Unknown', 'gene_primary': ''}
        
        # Organism
        org = self.elem.find(f"{self.ns}organism")
        if org is not None:
            for n in org.findall(f"{self.ns}name"):
                if n.get("type") == "scientific": out['organism_name'] = n.text
            tax = org.find(f"{self.ns}dbReference[@type='NCBI Taxonomy']")
            if tax is not None: out['organism_id'] = tax.get("id")
            
            lineage = org.find(f"{self.ns}lineage")
            if lineage is not None:
                taxa = [t.text for t in lineage.findall(f"{self.ns}taxon") if t.text]
                out['lineage'] = ", ".join(taxa)
                if taxa:
                    if "Bacteria" in taxa: out['kingdom'] = "Bacteria"
                    elif "Archaea" in taxa: out['kingdom'] = "Archaea"
                    elif "Viruses" in taxa: out['kingdom'] = "Viruses"
                    elif "Eukaryota" in taxa:
                        out['kingdom'] = "Fungi" if "Fungi" in taxa else "Plants" if "Viridiplantae" in taxa else "Metazoa" if "Metazoa" in taxa else "Eukaryota (Other)"
                    else: out['kingdom'] = taxa[0]

        # Gene
        gene = self.elem.find(f"{self.ns}gene")
        if gene is not None:
            for n in gene.findall(f"{self.ns}name"):
                if n.get("type") == "primary": out['gene_primary'] = n.text
                
        return out

    def _scan_references(self):
        out = {
            'ec_db': [], 'kw_ids': [], 'kw_text': [], 'go_ids': [], 'go_text': [], 'rhea_ids': [],
            'enzyme_dbs': [], 'xref_ensembl': [], 'xref_af': [], 'xref_kegg': [], 
            'xref_brenda': [], 'xref_biocyc': [], 'xref_pathcom': [], 'xref_reactome': [], 'xref_string': [], 'xref_unipath': [],
            'has_experimental_evidence': False, 'exp_level_db_ec': "NONE"
        }
        
        # Global Evidence for ECO:0000269
        for ev in self.elem.findall(f"{self.ns}evidence"):
             if ev.get("type") == "ECO:0000269":
                 out['has_experimental'] = True
        
        for key, ev_data in self.evidence_map.items():
             if ev_data['eco'] in ECO_STRICT_EXP or ev_data['eco'] in ECO_HTP_EXP:
                 out['has_experimental_evidence'] = True
                 break

        # Database Cross-References
        for child in self.elem.findall(f"{self.ns}dbReference"):
            d_type = child.get("type")
            d_id = child.get("id")
            if not d_id: continue
            
            if d_type in ["EC", "ENZYME"]: 
                out['ec_db'].append(d_id)
                level = self._evaluate_evidence_level(child.get("evidence", ""))
                out['exp_level_db_ec'] = self._max_evidence_level(out['exp_level_db_ec'], level)
                    
            elif d_type == "Rhea": out['rhea_ids'].append(d_id)
            elif d_type in ["BRENDA", "SABIO-RK"]: out['enzyme_dbs'].append(d_id) # "Reactome", "BioCyc", "PathwayCommons" => NOT ONLY ENZYMES
            elif d_type == "GO":
                term = child.find(f"{self.ns}property[@type='term']")
                val = term.get("value", "") if term is not None else ""
                
                # Evidence Code (ECO)
                ev_prop = child.find(f"{self.ns}property[@type='evidence']")
                ev_val = ev_prop.get("value", "") if ev_prop is not None else "Unknown"
                eco_code = ev_val.split("|")[0] if ev_val else "Unknown"

                out['go_ids'].append(d_id)
                if val: out['go_text'].append(f"{d_id}|{eco_code}|{val}")
            # Other external references (xrefs)    
            elif d_type == "Ensembl": out['xref_ensembl'].append(d_id)
            elif d_type == "AlphaFoldDB": out['xref_af'].append(d_id)
            elif d_type == "BRENDA": out['xref_brenda'].append(d_id)
            elif d_type == "BioCyc": out['xref_biocyc'].append(d_id)
            elif d_type == "PathwayCommons": out['xref_pathcom'].append(d_id)
            elif d_type == "KEGG": out['xref_kegg'].append(d_id)
            elif d_type == "Reactome": out['xref_reactome'].append(d_id)
            elif d_type == "STRING": out['xref_string'].append(d_id)
            elif d_type == "UniPathway": out['xref_unipath'].append(d_id)

        # Keywords
        for child in self.elem.findall(f"{self.ns}keyword"):
            kw_id, kw_text = child.get("id"), child.text
            if kw_id: out['kw_ids'].append(kw_id)
            if kw_text: out['kw_text'].append(kw_text)

        return out

    def _clean_chemical_names(self, side_string):
        """
        Removes transport/localization markers like (in), (out), or [cytoplasm].
        """
        # Remove markers like (in), (out), (external)
        cleaned = re.sub(r'\((in|out|external|internal)\)', '', side_string, flags=re.IGNORECASE)
        # Remove compartment markers in brackets
        cleaned = re.sub(r'\[.*?\]', '', cleaned)
        
        # Split, sort, and rejoin to allow stable comparison of chemical identity
        molecules = [m.strip().lower() for m in cleaned.split(' + ')]
        return " + ".join(sorted(filter(None, molecules)))

    def _process_sub_reactions(self, catalytic_text):
        """
        Splits alternative reactions and identifies passive transport.
        """
        sub_rxns = [r.strip() for r in catalytic_text.split('|')]
        processed_data = {'clean_equations': [], 'is_passive_transport': True}
        for rxn in sub_rxns:
            if " = " not in rxn: continue
            left, right = rxn.split(" = ", 1)
            
            # Compare cleaned versions to see if chemistry actually changes
            if self._clean_chemical_names(left) != self._clean_chemical_names(right):
                processed_data['is_passive_transport'] = False
            
            processed_data['clean_equations'].append((left, right))

        return processed_data

    def _extract_comments(self):
        """Extracts text blocks, reaction data, context, and kinetics."""
        res = {
            'text': {
                'function': [], 'catalytic_activity': [], 'pathway': [], 'caution': [],
                'cofactor': [], 'subunit': [], 'subcellular_location': [], 'induction': [],
                'similarity': [], 'interaction': [], 'disease': [], 'phenotype': [],
                'kinetics_km': [], 'kinetics_vmax': [], 'kinetics_kcat': []
            },
            'lists': {'substrates': [], 'products': [], 'rhea_ids': []},
            'has_kinetics': False,
            'is_passive_transport': None,
            'exp_level_catalytic_activity': "NONE"
        }

        for c in self.elem.findall(f"{self.ns}comment"):
            ctype = c.get("type")
            text = c.findtext(f"{self.ns}text") or ""
            
            # Standard Text Fields
            if ctype in ["function", "pathway", "cofactor", "subunit", "induction", "similarity", "caution"]:
                res['text'][ctype].append(text)
            
            # Subcellular Location
            elif ctype == "subcellular location":
                for loc_group in c.findall(f"{self.ns}subcellularLocation"):
                    loc = loc_group.find(f"{self.ns}location")
                    if loc is not None and loc.text: res['text']['subcellular_location'].append(loc.text)
            
            # Disease
            elif ctype == "disease":
                d_node = c.find(f"{self.ns}disease")
                if d_node is not None:
                    d_name = d_node.findtext(f"{self.ns}name") or ""
                    d_desc = d_node.findtext(f"{self.ns}description") or ""
                    res['text']['disease'].append(f"{d_name}: {d_desc}")
            
            # Phenotype
            elif ctype == "disruption phenotype":
                res['text']['phenotype'].append(text)

            # Interactions
            elif ctype == "interaction":
                pair = []
                for i_node in c.findall(f"{self.ns}interactant"):
                    i_id = i_node.findtext(f"{self.ns}id")
                    if i_id: pair.append(i_id)
                if pair: res['text']['interaction'].append(" <-> ".join(pair))

            # Catalytic Activity & Reactions
            elif ctype == "catalytic activity":
                # 1. Handle Standard Reaction and Sub-reactions
                reaction = c.find(f"{self.ns}reaction")
                if reaction is not None:
                    text_eq = reaction.findtext(f"{self.ns}text")
                    level = self._evaluate_evidence_level(reaction.get("evidence", ""))
                    res['exp_level_catalytic_activity'] = self._max_evidence_level(res['exp_level_catalytic_activity'], level)

                    if text_eq:
                        res['text']['catalytic_activity'].append(text_eq)

                        # Process sub-reactions for substrates/products and passive transport flags
                        rxn_info = self._process_sub_reactions(text_eq)
                        
                        if rxn_info['is_passive_transport']:
                            if res['is_passive_transport'] is None: res['is_passive_transport'] = True
                        else: res['is_passive_transport'] = False # Any active reaction
                        
                        # Extract individual molecules from all sub-reactions
                        for left, right in rxn_info['clean_equations']:
                            res['lists']['substrates'].extend([s.strip() for s in left.split(' + ')])
                            res['lists']['products'].extend([p.strip() for p in right.split(' + ')])

                    # Preserved: Extract Rhea IDs from the reaction
                    for dbref in reaction.findall(f"{self.ns}dbReference"):
                        if dbref.get("type") == "Rhea": res['lists']['rhea_ids'].append(dbref.get("id"))
                
                # 2. Preserved: Handle Physiological Reactions
                phys_rxn = c.find(f"{self.ns}physiologicalReaction")
                if phys_rxn is not None:
                    d_dir = phys_rxn.get("direction", "")
                    for dbref in phys_rxn.findall(f"{self.ns}dbReference"):
                        if dbref.get("type") == "Rhea": res['lists']['rhea_ids'].append(f"{dbref.get('id')}({d_dir})")

            # Kinetics & Biophysicochemical
            elif ctype == "biophysicochemical properties":
                kinetics = c.find(f"{self.ns}kinetics")
                if kinetics is not None:
                    res['has_kinetics'] = True
                    
                    # 1. Extract Raw KM
                    for km in kinetics.findall(f"{self.ns}KM"):
                        if km.text: res['text']['kinetics_km'].append(km.text.strip())
                    for vmax in kinetics.findall(f"{self.ns}Vmax"):
                        if vmax.text: res['text']['kinetics_vmax'].append(vmax.text.strip())
                    for text_node in kinetics.findall(f"{self.ns}text"):
                        if text_node.text and "kcat" in text_node.text: res['text']['kinetics_kcat'].append(text_node.text.strip())             

        if res['is_passive_transport'] is None: res['is_passive_transport'] = False
        for k in res['text']: res['text'][k] = " | ".join(filter(None, res['text'][k]))
        # Deduplicate lists
        for k in res['lists']: res['lists'][k] = list(set(filter(None, res['lists'][k])))
            
        return res

    def _extract_seq_features(self):
        """Extracts Active Sites, Binding Sites, and Mutagenesis data."""
        out = {'active_sites': [], 'binding_sites': [], 'mutagenesis': []}
        
        for f in self.elem.findall(f"{self.ns}feature"):
            ftype = f.get("type")
            if ftype in ["active site", "binding site", "mutagenesis site"]:
                desc = f.get("description", "")
                loc = f.find(f"{self.ns}location")
                pos_str = ""
                if loc is not None:
                    p = loc.find(f"{self.ns}position")
                    if p is not None: pos_str = p.get("position", "")
                    else:
                        b, e = loc.find(f"{self.ns}begin"), loc.find(f"{self.ns}end")
                        if b is not None and e is not None: pos_str = f"{b.get('position', '')}..{e.get('position', '')}"
                entry = f"{pos_str}({desc})" if desc else pos_str
                if ftype == "active site": out['active_sites'].append(entry)
                elif ftype == "binding site": out['binding_sites'].append(entry)
                elif ftype == "mutagenesis site": out['mutagenesis'].append(entry)
        
        return out


def parse_worker(xml_batch):
    """
    Worker function: Takes a list of XML strings, returns a list of dicts.
    Logs exceptions instead of failing silently.
    """
    logger = logging.getLogger("pipeline_logger")
    results = []
    
    for xml_str in xml_batch:
        try:
            parser = UniprotEntryParser(xml_str)
            data = parser.parse()
            if data:
                results.append(data)
        except Exception as e:
            logger.error(f"Failed to parse XML entry: {e}")
            
    return results
