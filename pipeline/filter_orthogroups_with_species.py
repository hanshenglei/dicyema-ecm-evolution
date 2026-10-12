#!/usr/bin/env python3
import argparse

def split_genes(cell):
    if not cell.strip():
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--species_col", required=True, help="Column name in Orthogroups.tsv, e.g. D_japonicum")
    parser.add_argument("--out_og", required=True)
    parser.add_argument("--out_ids", required=True)
    args = parser.parse_args()

    kept_ogs = 0
    kept_ids = set()

    with open(args.orthogroups) as f, open(args.out_og, "w") as out:
        header = f.readline().rstrip("\n").split("\t")

        if args.species_col not in header:
            raise RuntimeError(f"{args.species_col} not found in header: {header}")

        sp_idx = header.index(args.species_col)

        out.write("\t".join(header) + "\n")

        for line in f:
            row = line.rstrip("\n")
            if not row:
                continue

            fields = row.split("\t")

            # 补齐空列，防止尾部空白被 split 丢掉
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))

            target_genes = split_genes(fields[sp_idx])

            if not target_genes:
                continue

            kept_ogs += 1
            out.write("\t".join(fields) + "\n")

            for cell in fields[1:]:
                for g in split_genes(cell):
                    kept_ids.add(g)

    with open(args.out_ids, "w") as out:
        for g in sorted(kept_ids):
            out.write(g + "\n")

    print(f"Kept orthogroups with {args.species_col}: {kept_ogs}")
    print(f"All member IDs in kept orthogroups: {len(kept_ids)}")

if __name__ == "__main__":
    main()
