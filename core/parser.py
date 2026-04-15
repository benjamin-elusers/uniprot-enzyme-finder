# core/parser.py

import re
from lxml import etree
import logging

class UniprotEntryParser:
    """
    Parses a single UniProt XML entry into a comprehensive dictionary of raw features.
    Designed purely for ETL extraction (v22)
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
    
    def _build_evidence_usages(self):
        """Sweeps the XML tree to find where every evidence key was utilized."""
        usages = {}
        for node in self.elem.iter():
            ev_str = node.get("evidence")
            if ev_str:
                tag_name = node.tag.replace(self.ns, '')
                node_type = node.get('type')
                
                # Use the specific 'type' (e.g. 'subcellular location'), fallback to generic tag name
                annot_name = node_type if node_type else tag_name 
                
                for key in ev_str.split():
                    if key not in usages:
                        usages[key] = set()
                    usages[key].add(annot_name)
                    
        return {k: " | ".join(sorted(list(v))) for k, v in usages.items()}
    
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


    def _get_text(self, path):
        found = self.elem.find(f"{self.ns}{path}")
        return found.text if found is not None else ""

    def _extract_name(self):
        protein = self.elem.find(f"{self.ns}protein")
        name = ""
        ecs = []
        ev_keys = []
        
        if protein is not None:
            for n_type in ["recommendedName", "submittedName", "alternativeName"]:
                for group in protein.findall(f"{self.ns}{n_type}"):
                    fn = group.find(f"{self.ns}fullName")
                    if fn is not None and not name: name = fn.text
                    for ec in group.findall(f"{self.ns}ecNumber"):
                        if ec.text: 
                            ecs.append(ec.text)
                            if ec.get("evidence"): ev_keys.append(ec.get("evidence"))
                                
        return name, ecs, " ".join(ev_keys)

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
            'db_ec_evidence_keys': ""
        }
        
        # Database Cross-References        
        db_evs = []
        for child in self.elem.findall(f"{self.ns}dbReference"):
            d_type = child.get("type")
            d_id = child.get("id")
            if not d_id: continue
            
            if d_type in ["EC", "ENZYME"]: 
                out['ec_db'].append(d_id)
                if child.get("evidence"): db_evs.append(child.get("evidence"))
                    
            elif d_type == "Rhea": out['rhea_ids'].append(d_id)
            elif d_type in ["BRENDA", "SABIO-RK"]: out['enzyme_dbs'].append(d_id) 
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

        out['db_ec_evidence_keys'] = " ".join(db_evs)

        # Keywords
        for child in self.elem.findall(f"{self.ns}keyword"):
            kw_id, kw_text = child.get("id"), child.text
            if kw_id: out['kw_ids'].append(kw_id)
            if kw_text: out['kw_text'].append(kw_text)

        return out

    def _extract_comments(self):
        """Extracts text blocks, reaction data, context, and kinetics."""
        res = {
            'text': {
                'function': [], 'catalytic_activity': [], 'pathway': [], 'caution': [],
                'cofactor': [], 'subunit': [], 'subcellular_location': [], 'induction': [],
                'similarity': [], 'interaction': [], 'disease': [], 'phenotype': [],
                'kinetics_km': [], 'kinetics_vmax': [], 'kinetics_kcat': [],
                'reaction_evidence_keys': []
            },
            'lists': {'rhea_ids': []}
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
                    if text_eq: 
                        res['text']['catalytic_activity'].append(text_eq)

                    if reaction.get("evidence"): res['text']['reaction_evidence_keys'].append(reaction.get("evidence"))
                    for dbref in reaction.findall(f"{self.ns}dbReference"):
                        if dbref.get("type") == "Rhea": res['lists']['rhea_ids'].append(dbref.get("id"))
                
                # 2. Handle Physiological Reactions
                phys_rxn = c.find(f"{self.ns}physiologicalReaction")
                if phys_rxn is not None:
                    d_dir = phys_rxn.get("direction", "")
                    for dbref in phys_rxn.findall(f"{self.ns}dbReference"):
                        if dbref.get("type") == "Rhea": res['lists']['rhea_ids'].append(f"{dbref.get('id')}({d_dir})")

            # Kinetics & Biophysicochemical
            elif ctype == "biophysicochemical properties":
                kinetics = c.find(f"{self.ns}kinetics")
                if kinetics is not None:
                    # 1. Extract Raw KM                    
                    for km in kinetics.findall(f"{self.ns}KM"):
                        if km.text: res['text']['kinetics_km'].append(km.text.strip())
                    # 2. Extract Raw VMAX
                    for vmax in kinetics.findall(f"{self.ns}Vmax"):
                        if vmax.text: res['text']['kinetics_vmax'].append(vmax.text.strip())
                    # 3. Extract Raw kinetics free text (contains conditions, kcat, etc...)
                    for text_node in kinetics.findall(f"{self.ns}text"):
                        if text_node.text and "kcat" in text_node.text: res['text']['kinetics_kcat'].append(text_node.text.strip())             

        for k in res['text']: 
            res['text'][k] = " | ".join(filter(None, res['text'][k])) if isinstance(res['text'][k], list) else res['text'][k]

        # Deduplicate lists
        for k in res['lists']: 
            res['lists'][k] = list(set(filter(None, res['lists'][k])))
            
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

    def parse(self):
        self._build_evidence_dictionaries()
        
        data = {
            'meta': {},
            'features': {},   # Booleans/Flags
            'lists': {},      # Lists of IDs and Cross-references
            'text': {},       # Extracted text blocks
            'evidence_map': self.evidence_map, # Dictionary of evidences
            'evidence_usages': self._build_evidence_usages() # All tags linked to evidences
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
        
        # --- 3. Name ---
        full_name, name_ecs, name_ec_ev_keys = self._extract_name()
        data['text']['protein_name'] = full_name
        data['lists']['ec_from_name'] = name_ecs
        data['text']['name_ec_evidence_keys'] = name_ec_ev_keys

        # --- 4. Comments, Reactions & Kinetics ---
        comments = self._extract_comments()
        data['text'].update(comments['text'])
        data['lists'].update(comments['lists'])

        # --- 5. Database Cross-References ---
        refs = self._scan_references()
        for k, v in refs.items():
            if k in data['lists']: data['lists'][k].extend(v)
            elif isinstance(v, list): data['lists'][k] = v
            else: data['text'][k] = v # Catch raw evidence strings

        # --- 6. Sequence Features (Active Sites, Mutagenesis) ---
        feats = self._extract_seq_features()
        data['lists']['active_sites'] = feats['active_sites']
        data['lists']['binding_sites'] = feats['binding_sites']
        data['lists']['mutagenesis'] = feats['mutagenesis']
        data['lists']['ec_numbers'] = list(set(data['lists'].get('ec_db', []) + data['lists'].get('ec_from_name', [])))

        return data


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