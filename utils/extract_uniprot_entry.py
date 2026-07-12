import argparse
import gzip
import os
from lxml import etree

def get_accessions(arg):
    """
    Supports:
      - single accession
      - comma-separated accessions
      - text file with one accession per line
    """
    accessions = set()

    if os.path.exists(arg):
        print(f"[Info] Reading accessions from file: {arg}")
        with open(arg, "r") as f:
            for line in f:
                stripped = line.strip()
                if stripped:
                    accessions.add(stripped)
    else:
        parts = [x.strip() for x in arg.split(",") if x.strip()]
        accessions.update(parts)
        print(f"[Info] Loaded {len(accessions)} accession(s)")

    return accessions


def open_stream(path):
    """ Opens gzipped or plain XML file. """
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def scan_single_file(xml_path, remaining, found_accessions, outdir):
    """ Scan one XML file and write matching entries to their own files. """
    print(f"[Info] Scanning: {xml_path}")
    source = open_stream(xml_path)
    context = etree.iterparse(source, events=("end",), tag="{*}entry")
    files_written = 0

    for _, elem in context:
        entry_accessions = {acc.text.strip() for acc in elem.findall("{*}accession") if acc.text}
        matched = entry_accessions & remaining

        if matched:
            # Generate the XML string for this entry exactly once
            xml_str = etree.tostring(elem, encoding="unicode", pretty_print=True)

            # Write an individual XML file for each matched accession
            for acc in matched:
                out_path = os.path.join(outdir, f"{acc}.xml")
                with open(out_path, "w", encoding="utf-8") as out_f:
                    out_f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
                    out_f.write('<uniprot xmlns="http://uniprot.org/uniprot">\n')
                    out_f.write(xml_str)
                    out_f.write("</uniprot>\n")
                
                files_written += 1

            found_accessions.update(matched)
            remaining -= matched

            print(f"[Match] {', '.join(sorted(matched))} (remaining: {len(remaining)})")

            # Stop early if everything is found
            if not remaining:
                print("[Info] All requested accessions found.")
                elem.clear()
                break

        # Memory cleanup (keeps RAM footprint very low)
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]
            
    source.close()
    return files_written


def extract_entries(xml_paths, target_accessions, outdir):
    """
    Scan multiple UniProt XML files sequentially.
    Stops automatically once all IDs are found.
    """
    remaining = set(target_accessions)
    found_accessions = set()
    total_written = 0

    # Ensure output directory exists
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        print(f"[Info] Created output directory: {outdir}")

    for xml_path in xml_paths:
        if not remaining:
            break
        count = scan_single_file(xml_path, remaining, found_accessions, outdir)
        total_written += count

    print("\n========== SUMMARY ==========")
    print(f"Requested accessions : {len(target_accessions)}")
    print(f"Found accessions     : {len(found_accessions)}")
    print(f"Files written        : {total_written}")
    print(f"Output directory     : {os.path.abspath(outdir)}")

    missing = target_accessions - found_accessions
    if missing:
        print("\n[Warning] Missing accessions:")
        for acc in sorted(missing):
            print(f"  - {acc}")
    else:
        print("\n[Success] All accessions found.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract UniProt entries into individual XML files.")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="One or more UniProt XML/XML.GZ files")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory for the XML files (default: current dir)")
    parser.add_argument("-a", "--accession", required=True, 
                        help=("Single accession, comma-separated accessions, or text file with accessions"))

    args = parser.parse_args()
    targets = get_accessions(args.accession)
    
    extract_entries(xml_paths=args.input, target_accessions=targets, outdir=args.outdir)

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

