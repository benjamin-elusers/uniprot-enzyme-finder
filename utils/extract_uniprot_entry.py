import argparse
import gzip
import os
import sys
from lxml import etree

def get_accessions(arg):
    """
    Parses the accession argument. 
    If it's a file, reads lines. If it's a string, creates a set of one.
    """
    accessions = set()
    if os.path.exists(arg):
        print(f"[Info] Reading accessions from file: {arg}")
        with open(arg, 'r') as f:
            for line in f:
                stripped = line.strip()
                if stripped:
                    accessions.add(stripped)
    else:
        print(f"[Info] Searching for single accession: {arg}")
        accessions.add(arg.strip())
    return accessions

def extract_entries(xml_path, target_accessions, output_path):
    """
    Streams the XML file and writes matching entries to the output.
    """
    # Open input stream (gzip aware)
    if xml_path.endswith('.gz'):
        source = gzip.open(xml_path, 'rb')
    else:
        source = open(xml_path, 'rb')

    # Prepare output
    with open(output_path, 'w', encoding='utf-8') as out_f:
        # Write XML header and root tag to ensure valid XML output
        out_f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        out_f.write('<uniprot xmlns="http://uniprot.org/uniprot">\n')
        
        # Determine namespace from the file content usually, but we'll strip it for tag matching
        # or use the local name.
        
        # Iterparse context
        context = etree.iterparse(source, events=("end",), tag="{*}entry")
        
        found_count = 0
        
        print("[Info] Scanning XML stream...")
        for event, elem in context:
            # Check accessions within this entry
            # UniProt entries can have multiple accession tags (primary + secondaries)
            # We look for any match.
            
            # Find all accession tags regardless of namespace
            # using local-name() is safer if namespace varies (http vs https)
            entry_accessions = set()
            for child in elem:
                if child.tag.endswith('accession') and child.text:
                    entry_accessions.add(child.text.strip())
            
            # Intersection check
            if not entry_accessions.isdisjoint(target_accessions):
                # Write the XML string for this element
                xml_str = etree.tostring(elem, encoding='unicode', pretty_print=True)
                out_f.write(xml_str)
                found_count += 1
                
                # Optional: Stop early if we only looked for 1 item and found it
                if len(target_accessions) == 1 and found_count == 1:
                    print("[Info] Single target found. Stopping early.")
                    elem.clear()
                    break

            # Memory Management: Clear the element to keep RAM usage low
            elem.clear()
            # Also clear references from the root to this element
            while elem.getprevious() is not None:
                del elem.getparent()[0]

        # Close root tag
        out_f.write('</uniprot>\n')
    
    source.close()
    print(f"[Done] Extracted {found_count} entries to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract specific UniProt entries by Accession.")
    parser.add_argument("-i", "--input", required=True, help="Path to UniProt XML file (.xml or .xml.gz)")
    parser.add_argument("-a", "--accession", required=True, help="Single Accession ID OR Path to a text file with list of Accessions")
    parser.add_argument("-o", "--output", required=True, help="Output file path (e.g., extracted.xml)")
    
    args = parser.parse_args()
    
    targets = get_accessions(args.accession)
    extract_entries(args.input, targets, args.output)

#python extract_uniprot_entries.py -i uniprot_sprot.xml.gz -a P0C9F0 -o output.xml