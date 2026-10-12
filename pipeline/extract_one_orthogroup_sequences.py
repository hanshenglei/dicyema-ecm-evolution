#!/usr/bin/env python3
import argparse

def read_fasta(path):
    seq = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seq[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name:
        seq[name] = "".join(chunks)
    return seq

def split_genes(cell):
    if not cell.strip():
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--og", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    seqs = read_fasta(args.fasta)
    ids = []

    with open(args.orthogroups) as f:
        header = f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == args.og:
                for cell in fields[1:]:
                    ids.extend(split_genes(cell))
                break

    with open(args.out, "w") as out:
        for gid in ids:
            if gid in seqs:
                out.write(f">{gid}\n")
                s = seqs[gid]
                for i in range(0, len(s), 60):
                    out.write(s[i:i+60] + "\n")

    print(f"{args.og}: {len(ids)} IDs")
    print(f"Output: {args.out}")

if __name__ == "__main__":
    main()
