#!/usr/bin/env python3
import argparse

def split_genes(cell):
    if not cell.strip():
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]

def main():
    parser = argparse.ArgumentParser(
        description="Select ECM/Muscle orthogroups containing target species for tree/motif/domain analysis."
    )
    parser.add_argument("--ecm_og", required=True, help="ECM orthogroups with D_japonicum")
    parser.add_argument("--muscle_og", required=True, help="Muscle orthogroups with D_japonicum")
    parser.add_argument("--target_species", default="D_japonicum")
    parser.add_argument("--min_species", type=int, default=3)
    parser.add_argument("--max_total_genes", type=int, default=150)
    parser.add_argument("--min_target_genes", type=int, default=1)
    parser.add_argument("--out_list", required=True)
    parser.add_argument("--out_info", required=True)
    args = parser.parse_args()

    og_records = {}

    def read_og(path, category):
        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            if args.target_species not in header:
                raise RuntimeError(f"{args.target_species} not found in {path}")

            target_idx = header.index(args.target_species)

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

                if og not in og_records:
                    og_records[og] = {
                        "categories": set(),
                        "species_present": species_present,
                        "total_genes": total_genes,
                        "target_gene_count": len(target_genes),
                        "target_genes": target_genes
                    }

                og_records[og]["categories"].add(category)

    read_og(args.ecm_og, "ECM")
    read_og(args.muscle_og, "Muscle")

    selected = []

    with open(args.out_info, "w") as info:
        info.write(
            "Orthogroup\tcategory\tspecies_present\ttotal_genes\t"
            f"{args.target_species}_copy_number\t{args.target_species}_genes\tselected\treason\n"
        )

        for og in sorted(og_records):
            r = og_records[og]
            category = ",".join(sorted(r["categories"]))

            reasons = []

            if r["species_present"] < args.min_species:
                reasons.append("too_few_species")

            if r["total_genes"] > args.max_total_genes:
                reasons.append("too_many_genes")

            if r["target_gene_count"] < args.min_target_genes:
                reasons.append("too_few_target_genes")

            is_selected = len(reasons) == 0

            if is_selected:
                selected.append(og)

            info.write(
                f"{og}\t{category}\t{r['species_present']}\t{r['total_genes']}\t"
                f"{r['target_gene_count']}\t{','.join(r['target_genes'])}\t"
                f"{'yes' if is_selected else 'no'}\t"
                f"{';'.join(reasons) if reasons else 'PASS'}\n"
            )

    with open(args.out_list, "w") as out:
        for og in selected:
            out.write(og + "\n")

    print(f"Total unique OGs: {len(og_records)}")
    print(f"Selected OGs: {len(selected)}")
    print(f"Output list: {args.out_list}")
    print(f"Output info: {args.out_info}")

if __name__ == "__main__":
    main()
