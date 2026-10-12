#!/usr/bin/env python3
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--seed_ids", required=True)
    parser.add_argument("--out_og", required=True)
    parser.add_argument("--out_ids", required=True)
    parser.add_argument("--out_seed_map", required=True)
    args = parser.parse_args()

    seed = set(x.strip() for x in open(args.seed_ids) if x.strip())

    selected_ogs = []
    all_member_ids = set()

    with open(args.orthogroups) as f, \
         open(args.out_og, "w") as out_og, \
         open(args.out_seed_map, "w") as out_map:

        header = f.readline().rstrip("\n")
        out_og.write(header + "\n")
        out_map.write("Orthogroup\tseed_ids\n")

        for line in f:
            row = line.rstrip("\n")
            if not row:
                continue

            fields = row.split("\t")
            og = fields[0]

            members = []
            for cell in fields[1:]:
                if cell.strip():
                    members.extend([x.strip() for x in cell.split(",") if x.strip()])

            hit = sorted(seed.intersection(members))

            if hit:
                out_og.write(row + "\n")
                out_map.write(og + "\t" + ",".join(hit) + "\n")
                selected_ogs.append(og)
                all_member_ids.update(members)

    with open(args.out_ids, "w") as out:
        for x in sorted(all_member_ids):
            out.write(x + "\n")

    print(f"Seed IDs: {len(seed)}")
    print(f"Selected orthogroups: {len(selected_ogs)}")
    print(f"All member IDs in selected orthogroups: {len(all_member_ids)}")

if __name__ == "__main__":
    main()
