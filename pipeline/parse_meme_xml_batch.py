#!/usr/bin/env python3
import argparse
import os
import xml.etree.ElementTree as ET

def parse_one(xml_file, out_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    motifs = {}

    for motif in root.findall(".//motif"):
        mid = motif.attrib.get("id", motif.attrib.get("name", "NA"))
        alt = motif.attrib.get("alt", mid)
        width = int(motif.attrib.get("width", "0"))
        motifs[mid] = (alt, width)

    seq_name = {}
    for seq in root.findall(".//sequence"):
        sid = seq.attrib.get("id", seq.attrib.get("name", "NA"))
        name = seq.attrib.get("name", sid)
        seq_name[sid] = name

    hits = []

    for scanned_sites in root.findall(".//scanned_sites"):
        seq_id = scanned_sites.attrib.get("sequence_id", scanned_sites.attrib.get("id", "NA"))
        gene = seq_name.get(seq_id, seq_id)

        for hit in scanned_sites.findall(".//scanned_site"):
            motif_id = hit.attrib.get("motif_id", "NA")
            pos = int(hit.attrib.get("position", "0")) + 1
            pvalue = hit.attrib.get("pvalue", "NA")
            motif_name, width = motifs.get(motif_id, (motif_id, 0))
            hits.append((gene, motif_name, pos, pos + width - 1, pvalue))

    # 兼容部分 MEME XML 旧格式
    if not hits:
        for seq in root.findall(".//sequence"):
            gene = seq.attrib.get("name", seq.attrib.get("id", "NA"))

            for hit in seq.findall(".//contributing_site"):
                motif_id = hit.attrib.get("motif_id", "NA")
                pos = int(hit.attrib.get("position", "0")) + 1
                pvalue = hit.attrib.get("pvalue", "NA")
                motif_name, width = motifs.get(motif_id, (motif_id, 0))
                hits.append((gene, motif_name, pos, pos + width - 1, pvalue))

    with open(out_file, "w") as out:
        out.write("gene_id\tmotif\tstart\tend\tpvalue\n")
        for row in hits:
            out.write("\t".join(map(str, row)) + "\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root_dir", required=True)
    args = parser.parse_args()

    parsed = 0
    failed = []

    for og in sorted(os.listdir(args.root_dir)):
        og_dir = os.path.join(args.root_dir, og)
        if not os.path.isdir(og_dir):
            continue

        xml = os.path.join(og_dir, "meme", "meme.xml")
        out = os.path.join(og_dir, f"{og}.motif.tsv")

        if not os.path.exists(xml):
            continue

        try:
            parse_one(xml, out)
            parsed += 1
        except Exception as e:
            failed.append((og, str(e)))

    print(f"Parsed MEME XML for OGs: {parsed}")

    if failed:
        print("Failed OGs:")
        for og, err in failed:
            print(og, err)

if __name__ == "__main__":
    main()
