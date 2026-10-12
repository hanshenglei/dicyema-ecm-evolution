#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict

def read_fasta(path):
    records = []
    header = None
    seqs = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seqs)))
                header = line
                seqs = []
            else:
                seqs.append(line.strip())

    if header is not None:
        records.append((header, "".join(seqs)))

    return records

def write_fasta(records, path, width=60):
    with open(path, "w") as out:
        for header, seq in records:
            out.write(header + "\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")

def extract_gene_id(header):
    m = re.search(r"\[gene=([^\]]+)\]", header)
    if m:
        return m.group(1)

    m = re.search(r"\[locus_tag=([^\]]+)\]", header)
    if m:
        return m.group(1)

    m = re.search(r"\[db_xref=GeneID:([^\]]+)\]", header)
    if m:
        return "GeneID_" + m.group(1)

    return header[1:].split()[0]

def extract_protein_id(header):
    m = re.search(r"\[protein_id=([^\]]+)\]", header)
    if m:
        return m.group(1)
    return header[1:].split()[0]

def is_partial(header):
    return "[partial=" in header

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--prefer_complete", action="store_true")
    args = parser.parse_args()

    records = read_fasta(args.input)
    groups = defaultdict(list)

    for header, seq in records:
        gene_id = extract_gene_id(header)
        groups[gene_id].append({
            "header": header,
            "seq": seq.upper(),
            "length": len(seq),
            "protein_id": extract_protein_id(header),
            "partial": is_partial(header)
        })

    kept = []

    with open(args.map, "w") as out:
        out.write("gene_id\tkept_protein_id\tkept_length\tisoform_count\tpartial_status\tkept_header\n")

        for gene_id in sorted(groups):
            items = groups[gene_id]

            if args.prefer_complete:
                complete = [x for x in items if not x["partial"]]
                if complete:
                    best = max(complete, key=lambda x: x["length"])
                else:
                    best = max(items, key=lambda x: x["length"])
            else:
                best = max(items, key=lambda x: x["length"])

            kept.append((best["header"], best["seq"]))

            out.write(
                f"{gene_id}\t{best['protein_id']}\t{best['length']}\t"
                f"{len(items)}\t"
                f"{'partial' if best['partial'] else 'complete'}\t"
                f"{best['header']}\n"
            )

    write_fasta(kept, args.output)

    print(f"Input records: {len(records)}")
    print(f"Genes: {len(groups)}")
    print(f"Output records: {len(kept)}")
    print(f"Output: {args.output}")

if __name__ == "__main__":
    main()
