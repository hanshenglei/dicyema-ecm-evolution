#!/usr/bin/env python3
import argparse

def split_genes(cell):
    if not cell.strip():
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--category", required=True)
    parser.add_argument("--target_species", default="D_japonicum")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    with open(args.orthogroups) as f, open(args.out, "w") as out:
        header = f.readline().rstrip("\n").split("\t")

        if args.target_species not in header:
            raise RuntimeError(f"{args.target_species} not found in Orthogroups header")

        target_idx = header.index(args.target_species)

        out.write(
            "category\tOrthogroup\thas_D_japonicum\tD_japonicum_copy_number\t"
            "total_species_present\ttotal_gene_number\n"
        )

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))

            og = fields[0]
            target_genes = split_genes(fields[target_idx])

            species_present = 0
            total_genes = 0

            for cell in fields[1:]:
                genes = split_genes(cell)
                if genes:
                    species_present += 1
                    total_genes += len(genes)

            has_target = "yes" if target_genes else "no"

            out.write(
                f"{args.category}\t{og}\t{has_target}\t{len(target_genes)}\t"
                f"{species_present}\t{total_genes}\n"
            )

if __name__ == "__main__":
    main()
