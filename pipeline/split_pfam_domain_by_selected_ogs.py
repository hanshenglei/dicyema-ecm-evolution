#!/usr/bin/env python3
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root_dir", required=True)
    parser.add_argument("--pfam_domain", required=True)
    args = parser.parse_args()

    domain_by_gene = {}

    with open(args.pfam_domain) as f:
        header = f.readline().rstrip("\n")
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            gene_id = line.split("\t", 1)[0]
            domain_by_gene.setdefault(gene_id, []).append(line)

    n_og = 0

    for og in sorted(os.listdir(args.root_dir)):
        og_dir = os.path.join(args.root_dir, og)
        if not os.path.isdir(og_dir):
            continue

        ids_file = os.path.join(og_dir, f"{og}.ids")
        if not os.path.exists(ids_file):
            continue

        ids = [x.strip() for x in open(ids_file) if x.strip()]
        out_file = os.path.join(og_dir, f"{og}.pfam.domain.tsv")

        with open(out_file, "w") as out:
            out.write(header + "\n")
            for gid in ids:
                for line in domain_by_gene.get(gid, []):
                    out.write(line + "\n")

        n_og += 1

    print(f"Processed OGs: {n_og}")

if __name__ == "__main__":
    main()
