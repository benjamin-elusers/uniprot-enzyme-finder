# core/parser.py

import re
from lxml import etree
from config.patterns import (
    EC_NAME_PATTERN, HARD_INACTIVE_PATTERN, SOFT_INACTIVE_PATTERN,
    BIOSYNTH_PATTERN, SEC_MET_PATTERNS
)

class UniprotEntryParser:
    """
    Parses a single UniProt XML entry into a dictionary of raw features.
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

    def parse(self):
        data = {
            'meta': {},
            'features': {},   # Booleans/Flags
            'lists': {},      # Lists of IDs
            'text': {}        # Extracted text blocks
        }

        # --- 1. Header Info ---
        data['meta']['accession'] = self._get_text("accession")
        data['meta']['id'] = self._get_text("name")
        data['meta']['reviewed'] = (self.elem.get("dataset") == "Swiss-Prot")
        
        # Protein Existence
        pe = self.elem.find(f"{self.ns}proteinExistence")
        data['meta']['protein_existence'] = pe.get("type") if pe is not None else "uncertain"

        # Sequence Attributes (Fragment / Length)
        seq = self.elem.find(f"{self.ns}sequence")
        data['meta']['is_fragment'] = False
        if seq is not None:
            if 'fragment' in seq.attrib:
                data['meta']['is_fragment'] = (seq.get("fragment") != "single") 
        
        # --- 2. Name & Inactive Status ---
        full_name = self._extract_name()
        data['text']['protein_name'] = full_name
        
        is_hard_inactive = bool(full_name and HARD_INACTIVE_PATTERN.search(full_name))
        is_soft_inactive = bool(full_name and SOFT_INACTIVE_PATTERN.search(full_name))
        data['features']['is_inactive_text'] = is_hard_inactive or is_soft_inactive

        # --- 3. Comments & Kinetics ---
        comments = self._extract_comments()
        data['text'].update(comments['text'])
        data['features']['has_kinetics'] = comments['has_kinetics']
        
        # --- 4. Evidence Scanning (EC, GO, Keywords, Experimental) ---
        refs = self._scan_references()
        data['lists'].update(refs)
        data['features']['has_experimental_evidence'] = refs['has_experimental']
        
        # --- 5. Sequence Features (Active Sites) ---
        feats = self._extract_seq_features()
        data['lists']['active_sites'] = feats['active_sites']
        data['lists']['binding_sites'] = feats['binding_sites']
        data['features']['has_active_site'] = bool(feats['active_sites'])
        data['features']['has_binding_site'] = bool(feats['binding_sites'])

        # --- 6. Text Mining Flags ---
        combined_text = (
            f"{data['text'].get('pathway', '')} {data['text'].get('function', '')} "
            f"{data['text'].get('catalytic_activity', '')} {full_name} "
        )
        
        data['features']['has_text_mining'] = bool(SEC_MET_PATTERNS.search(combined_text))
        
        # Check for EC in name (if not found in DB)
        if EC_NAME_PATTERN.search(combined_text):
            found_ecs = EC_NAME_PATTERN.findall(combined_text)
            cleaned_ecs = [x.strip("()EC ") for x in found_ecs]
            data['lists']['ec_numbers'].extend(cleaned_ecs)
            data['lists']['ec_numbers'] = list(set(data['lists']['ec_numbers']))

        return data

    # --- HELPERS ---

    def _get_text(self, path):
        found = self.elem.find(f"{self.ns}{path}")
        return found.text if found is not None else ""

    def _extract_name(self):
        protein = self.elem.find(f"{self.ns}protein")
        if protein is not None:
            for n_type in ["recommendedName", "submittedName"]:
                group = protein.find(f"{self.ns}{n_type}")
                if group is not None:
                    fn = group.find(f"{self.ns}fullName")
                    if fn is not None: return fn.text
        return ""

    def _scan_references(self):
        out = {
            'ec_numbers': [], 'kw_ids': [], 'go_ids': [], 'rhea_ids': [],
            'has_experimental': False
        }
        
        # Evidence Tags (ECO codes)
        # ECO:0000269 = Experimental evidence
        for ev in self.elem.findall(f"{self.ns}evidence"):
             if ev.get("type") == "ECO:0000269":
                 out['has_experimental'] = True

        # DB Refs
        for child in self.elem.findall(f"{self.ns}dbReference"):
            d_type = child.get("type")
            d_id = child.get("id")
            if not d_id: continue
            
            if d_type == "EC": out['ec_numbers'].append(d_id)
            elif d_type == "Rhea": out['rhea_ids'].append(d_id)
            elif d_type == "GO": out['go_ids'].append(d_id)

        # Keywords
        for child in self.elem.findall(f"{self.ns}keyword"):
            kw_id = child.get("id")
            if kw_id: out['kw_ids'].append(kw_id)

        return out

    def _extract_comments(self):
        res = {
            'text': {'function': [], 'catalytic_activity': [], 'pathway': []},
            'has_kinetics': False
        }
        
        for c in self.elem.findall(f"{self.ns}comment"):
            ctype = c.get("type")
            text = c.findtext(f"{self.ns}text") or ""
            
            if ctype == "function": res['text']['function'].append(text)
            elif ctype == "pathway": res['text']['pathway'].append(text)
            elif ctype == "catalytic activity":
                reaction = c.find(f"{self.ns}reaction")
                if reaction is not None:
                    eq = reaction.findtext(f"{self.ns}text")
                    if eq: res['text']['catalytic_activity'].append(eq)
            elif ctype == "biophysicochemical properties":
                kin = c.find(f"{self.ns}kinetics")
                if kin is not None:
                    if (kin.find(f"{self.ns}KM") is not None or 
                        kin.find(f"{self.ns}Vmax") is not None or
                        kin.find(f"{self.ns}text") is not None):
                        res['has_kinetics'] = True

        for k in res['text']:
            res['text'][k] = " | ".join(filter(None, res['text'][k]))
            
        return res

    def _extract_seq_features(self):
        out = {'active_sites': [], 'binding_sites': []}
        for f in self.elem.findall(f"{self.ns}feature"):
            ftype = f.get("type")
            if ftype in ["active site", "binding site"]:
                desc = f.get("description", "")
                loc = f.find(f"{self.ns}location")
                pos_str = ""
                if loc is not None:
                    p = loc.find(f"{self.ns}position")
                    if p is not None: pos_str = p.get("position")
                    else:
                        b, e = loc.find(f"{self.ns}begin"), loc.find(f"{self.ns}end")
                        if b is not None and e is not None: 
                            pos_str = f"{b.get('position')}..{e.get('position')}"
                entry = f"{pos_str}({desc})" if desc else pos_str
                if ftype == "active site": out['active_sites'].append(entry)
                elif ftype == "binding site": out['binding_sites'].append(entry)
        return out