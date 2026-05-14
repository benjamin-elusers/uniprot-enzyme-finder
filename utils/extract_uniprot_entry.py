import argparse
import gzip
import hashlib
import os
import sys
from lxml import etree

def get_accessions(arg):
    """
    Supports:
      - single accession
      - comma-separated accessions
      - text file with one accession per line
    """
    accessions = set()

    # File input
    if os.path.exists(arg):
        print(f"[Info] Reading accessions from file: {arg}")
        with open(arg, "r") as f:
            for line in f:
                stripped = line.strip()
                if stripped:
                    accessions.add(stripped)
    else:
        # Comma-separated OR single
        parts = [x.strip() for x in arg.split(",") if x.strip()]
        accessions.update(parts)
        print(f"[Info] Loaded {len(accessions)} accession(s)")

        #print(f"[Info] Searching for single accession: {arg}")
        #accessions.add(arg.strip())
    return accessions

def make_output_name(accessions, output=None):
    """ Generate deterministic output filename if output not supplied. """

    if output:
        return output

    sorted_ids = sorted(accessions)
    # Short readable name for small sets
    if len(sorted_ids) <= 5:
        base = "_".join(sorted_ids)
        # Prevent absurdly long filenames
        if len(base) > 120:
            base = base[:120]

        return f"{base}.xml"

    # Deterministic hash-based name for larger sets
    joined = ",".join(sorted_ids)
    signature = hashlib.md5(joined.encode()).hexdigest()[:10]

    return f"uniprot_extract_n{len(sorted_ids)}_{signature}.xml"



def open_stream(path):
    """ Opens gzipped or plain XML file. """
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def scan_single_file(xml_path, remaining, out_f, found_accessions):
    """ Scan one XML file and write matching entries. """

    print(f"[Info] Scanning: {xml_path}")
    source = open_stream(xml_path)
    context = etree.iterparse(source, events=("end",), tag="{*}entry")
    entries_found = 0

    for _, elem in context:

        entry_accessions = {acc.text.strip() for acc in elem.findall("{*}accession") if acc.text }
        matched = entry_accessions & remaining

        if matched:

            xml_str = etree.tostring(elem,encoding="unicode",pretty_print=True)
            out_f.write(xml_str)
            found_accessions.update(matched)
            remaining -= matched
            entries_found += 1

            print(f"[Match] {', '.join(sorted(matched))} remaining: {len(remaining)})")

            # Stop early if everything found
            if not remaining:
                print("[Info] All requested accessions found.")
                elem.clear()
                break

        # Memory cleanup
        elem.clear()

        while elem.getprevious() is not None:
            del elem.getparent()[0]
    source.close()

    return entries_found


def extract_entries(xml_paths, target_accessions, output_path):
    """
    Scan multiple UniProt XML files sequentially.
    Stops automatically once all IDs are found.
    """

    remaining = set(target_accessions)
    found_accessions = set()
    total_entries = 0

    with open(output_path, "w", encoding="utf-8") as out_f:

        # XML header
        out_f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        out_f.write('<uniprot xmlns="http://uniprot.org/uniprot">\n')
        for xml_path in xml_paths:

            if not remaining:
                break
            count = scan_single_file(xml_path, remaining, out_f, found_accessions)
            total_entries += count

        out_f.write("</uniprot>\n")

    print("\n========== SUMMARY ==========")
    print(f"Requested accessions : {len(target_accessions)}")
    print(f"Found accessions     : {len(found_accessions)}")
    print(f"Entries written      : {total_entries}")
    print(f"Output file          : {output_path}")

    missing = target_accessions - found_accessions
    if missing:
        print("\n[Warning] Missing accessions:")
        for acc in sorted(missing):
            print(f"  - {acc}")
    else:
        print("\n[Success] All accessions found.")

#def extract_entries(xml_path, target_accessions, output_path):
#    """
#    Streams the XML file and writes matching entries to the output.
#    """
#    # Open input stream (gzip aware)
#    if xml_path.endswith('.gz'):
#        source = gzip.open(xml_path, 'rb')
#    else:
#        source = open(xml_path, 'rb')
#
#    # Prepare output
#    with open(output_path, 'w', encoding='utf-8') as out_f:
#        # Write XML header and root tag to ensure valid XML output
#        out_f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
#        out_f.write('<uniprot xmlns="http://uniprot.org/uniprot">\n')
#        
#        # Determine namespace from the file content usually, but we'll strip it for tag matching
#        # or use the local name.
#        
#        # Iterparse context
#        context = etree.iterparse(source, events=("end",), tag="{*}entry")
#        
#        found_count = 0
#        
#        print("[Info] Scanning XML stream...")
#        for event, elem in context:
#            # Check accessions within this entry
#            # UniProt entries can have multiple accession tags (primary + secondaries)
#            # We look for any match.
#            
#            # Find all accession tags regardless of namespace
#            # using local-name() is safer if namespace varies (http vs https)
#            entry_accessions = set()
#            for child in elem:
#                if child.tag.endswith('accession') and child.text:
#                    entry_accessions.add(child.text.strip())
#            
#            # Intersection check
#            if not entry_accessions.isdisjoint(target_accessions):
#                # Write the XML string for this element
#                xml_str = etree.tostring(elem, encoding='unicode', pretty_print=True)
#                out_f.write(xml_str)
#                found_count += 1
#                
#                # Optional: Stop early if we only looked for 1 item and found it
#                if len(target_accessions) == 1 and found_count == 1:
#                    print("[Info] Single target found. Stopping early.")
#                    elem.clear()
#                    break
#
#            # Memory Management: Clear the element to keep RAM usage low
#            elem.clear()
#            # Also clear references from the root to this element
#            while elem.getprevious() is not None:
#                del elem.getparent()[0]
#
#        # Close root tag
#        out_f.write('</uniprot>\n')
#    
#    source.close()
#    print(f"[Done] Extracted {found_count} entries to {output_path}")



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Extract UniProt entries from one or more XML/XML.GZ files.")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="One or more UniProt XML/XML.GZ files")
    parser.add_argument("-o", "--output", default=None, help="Output XML file (optional)")
    parser.add_argument("-a", "--accession", required=True, 
                        help=("Single accession, comma-separated accessions, or text file with accessions"))

    args = parser.parse_args()
    targets = get_accessions(args.accession)
    output_path = make_output_name(targets,args.output)
    extract_entries(xml_paths=args.input, target_accessions=targets, output_path=output_path)

"""
Examples
--------

Single accession:
python -m utils.extract_uniprot_entry \
    -i data/uniprot_sprot.xml.gz \
    -a Q12084

Multiple accessions:
python -m utils.extract_uniprot_entry \
    -i data/uniprot_sprot.xml.gz data/uniprot_trembl.xml.gz \
    -a Q12084,P32178,Q9BUT1

File of accessions:
python -m utils.extract_uniprot_entry \
    -i data/uniprot_sprot.xml.gz data/uniprot_trembl.xml.gz \
    -a ids.txt

Custom output:
python -m utils.extract_uniprot_entry \
    -i data/uniprot_sprot.xml.gz \
    -a Q12084,P32178 \
    -o my_subset.xml
"""

