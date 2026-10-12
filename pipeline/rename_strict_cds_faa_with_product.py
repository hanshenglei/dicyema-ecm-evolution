#!/usr/bin/env python3
import argparse
import re

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
                header = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line.strip())

    if header is not None:
        records.append((header, "".join(seqs)))

    return records

def write_fasta(records, path, width=60):
    with open(path, "w") as out:
        for seq_id, seq in records:
            out.write(f">{seq_id}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")

def get_first_id_from_header(header):
    if header.startswith(">"):
        return header[1:].split()[0]
    return header.split()[0]

def extract_tag(header, tag):
    m = re.search(rf"\[{tag}=([^\]]+)\]", header)
    if m:
        return m.group(1)
    return "NA"

def simplify_product_name(name):
    if name == "NA":
        return "NA"

    # 适度清理，保留生物学含义
    name = name.replace(",", "")
    name = re.sub(r"\s+", " ", name).strip()

    # 去掉过长的 isoform 信息，避免画图太乱
    name = re.sub(r"\s+isoform\s+X\d+$", "", name)
    name = re.sub(r"\s+isoform\s+\S+$", "", name)

    return name

def read_longest_map(path):
    """
    读取 keep_longest_cds_by_gene.py 输出的 map:
    gene_id kept_protein_id kept_length isoform_count partial_status kept_header
    """
    oldid_to_info = {}

    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")

            gene_id = parts[col.get("gene_id")]
            kept_protein_id = parts[col.get("kept_protein_id")]
            kept_header = parts[col.get("kept_header")]

            old_id = get_first_id_from_header(kept_header)

            protein_name = extract_tag(kept_header, "protein")
            protein_id = extract_tag(kept_header, "protein_id")
            gene_from_header = extract_tag(kept_header, "gene")
            locus_tag = extract_tag(kept_header, "locus_tag")

            if protein_id == "NA":
                protein_id = kept_protein_id

            if gene_from_header != "NA":
                final_gene_id = gene_from_header
            elif locus_tag != "NA":
                final_gene_id = locus_tag
            else:
                final_gene_id = gene_id

            oldid_to_info[old_id] = {
                "gene_id": final_gene_id,
                "protein_id": protein_id,
                "protein_name": simplify_product_name(protein_name),
                "full_header": kept_header
            }

    return oldid_to_info

def main():
    parser = argparse.ArgumentParser(
        description="Rename strict CDS/FAA and keep original protein product annotation."
    )
    parser.add_argument("--cds", required=True, help="Strict complete CDS fasta")
    parser.add_argument("--faa", required=True, help="Strict complete FAA fasta")
    parser.add_argument("--longest_map", required=True, help="Map from keep_longest_cds_by_gene.py")
    parser.add_argument("--species", required=True, help="Full species name, e.g. Octopus_sinensis")
    parser.add_argument("--prefix", required=True, help="New ID prefix, e.g. O_sinensis")
    parser.add_argument("--out_cds", required=True)
    parser.add_argument("--out_faa", required=True)
    parser.add_argument("--out_map", required=True)
    args = parser.parse_args()

    cds_records = read_fasta(args.cds)
    faa_records = read_fasta(args.faa)

    if len(cds_records) != len(faa_records):
        raise RuntimeError(
            f"CDS and FAA record numbers differ: {len(cds_records)} vs {len(faa_records)}"
        )

    longest_info = read_longest_map(args.longest_map)

    renamed_cds = []
    renamed_faa = []

    with open(args.out_map, "w") as m:
        m.write(
            "species\told_id\tnew_id\tgene_id\tprotein_id\tprotein_name\tdisplay_label\tfull_header\n"
        )

        for idx, ((old_cds_id, cds_seq), (old_faa_id, faa_seq)) in enumerate(zip(cds_records, faa_records), 1):
            if old_cds_id != old_faa_id:
                raise RuntimeError(
                    f"CDS and FAA IDs not matched at record {idx}: {old_cds_id} vs {old_faa_id}"
                )

            old_id = old_cds_id
            new_id = f"{args.prefix}{idx:05d}"

            info = longest_info.get(old_id, None)

            if info is None:
                gene_id = "NA"
                protein_id = "NA"
                protein_name = "NA"
                full_header = "NA"
            else:
                gene_id = info["gene_id"]
                protein_id = info["protein_id"]
                protein_name = info["protein_name"]
                full_header = info["full_header"]

            if protein_name != "NA":
                display_label = f"{new_id} | {protein_name}"
            else:
                display_label = new_id

            renamed_cds.append((new_id, cds_seq))
            renamed_faa.append((new_id, faa_seq))

            m.write(
                f"{args.species}\t{old_id}\t{new_id}\t{gene_id}\t{protein_id}\t"
                f"{protein_name}\t{display_label}\t{full_header}\n"
            )

    write_fasta(renamed_cds, args.out_cds)
    write_fasta(renamed_faa, args.out_faa)

    print("Done.")
    print(f"Species: {args.species}")
    print(f"Input records: {len(cds_records)}")
    print(f"Output CDS: {args.out_cds}")
    print(f"Output FAA: {args.out_faa}")
    print(f"Map: {args.out_map}")

if __name__ == "__main__":
    main()
