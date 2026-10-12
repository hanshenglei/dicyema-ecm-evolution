#!/usr/bin/env python3
import argparse
import os

def nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root_dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.out, "w") as out:
        out.write(
            "Orthogroup\tmember_count\tfaa\ttreefile\tmotif\tdomain\t"
            "annotation\tseed_summary\tready_for_plot\n"
        )

        for og in sorted(os.listdir(args.root_dir)):
            og_dir = os.path.join(args.root_dir, og)
            if not os.path.isdir(og_dir):
                continue

            faa = os.path.join(og_dir, f"{og}.faa")
            tree = os.path.join(og_dir, f"{og}.mafft.trim.fa.treefile")
            motif = os.path.join(og_dir, f"{og}.motif.tsv")
            domain = os.path.join(og_dir, f"{og}.pfam.domain.tsv")
            ann = os.path.join(og_dir, f"{og}.annotation.tsv")
            seed = os.path.join(og_dir, f"{og}.seed_summary.tsv")

            member_count = 0
            if os.path.exists(faa):
                with open(faa) as f:
                    member_count = sum(1 for line in f if line.startswith(">"))

            flags = {
                "faa": nonempty(faa),
                "tree": nonempty(tree),
                "motif": nonempty(motif),
                "domain": nonempty(domain),
                "annotation": nonempty(ann),
                "seed": nonempty(seed)
            }

            ready = flags["faa"] and flags["tree"] and flags["domain"] and flags["annotation"]

            out.write(
                f"{og}\t{member_count}\t{flags['faa']}\t{flags['tree']}\t"
                f"{flags['motif']}\t{flags['domain']}\t{flags['annotation']}\t"
                f"{flags['seed']}\t{ready}\n"
            )

if __name__ == "__main__":
    main()
