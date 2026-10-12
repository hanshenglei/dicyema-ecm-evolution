#!/usr/bin/env python3
import argparse

def split_genes(cell):
    if not cell.strip():
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--species_col", required=True)
    parser.add_argument("--out_ids", required=True)
    parser.add_argument("--out_map", required=True)
    args = parser.parse_args()

    ids = []

    with open(args.orthogroups) as f, open(args.out_map, "w") as m:
        header = f.readline().rstrip("\n").split("\t")

        if args.species_col not in header:
            raise RuntimeError(f"{args.species_col} not found in header")

        sp_idx = header.index(args.species_col)

        m.write("Orthogroup\tgene_id\n")

        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))

            og = fields[0]
            genes = split_genes(fields[sp_idx])

            for g in genes:
                ids.append(g)
                m.write(f"{og}\t{g}\n")

    with open(args.out_ids, "w") as out:
        for g in sorted(set(ids)):
            out.write(g + "\n")

    print(f"{args.species_col} genes: {len(set(ids))}")

if __name__ == "__main__":
    main()
