#!/usr/bin/env python3
import argparse
from collections import defaultdict

def read_annotation_map(path):
    ann = {}
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {x:i for i,x in enumerate(header)}

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            new_id = parts[col["new_id"]]
            species = parts[col["species"]]
            pname = parts[col["protein_name"]]
            ann[new_id] = (species, pname)
    return ann

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed_map", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--category", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ann = read_annotation_map(args.annotation)

    with open(args.seed_map) as f, open(args.out, "w") as out:
        header = f.readline()
        out.write("category\tOrthogroup\tseed_count\tseed_species\tseed_names\tseed_ids\n")

        for line in f:
            og, seed_str = line.rstrip("\n").split("\t")
            seeds = [x for x in seed_str.split(",") if x]

            species_set = []
            name_set = []

            for sid in seeds:
                sp, pname = ann.get(sid, ("NA", "NA"))
                species_set.append(sp)
                name_set.append(pname)

            species_u = sorted(set(species_set))
            names_u = sorted(set(name_set))

            out.write(
                f"{args.category}\t{og}\t{len(seeds)}\t"
                f"{';'.join(species_u)}\t"
                f"{';'.join(names_u)}\t"
                f"{';'.join(seeds)}\n"
            )

if __name__ == "__main__":
    main()
