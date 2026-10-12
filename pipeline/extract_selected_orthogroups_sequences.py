#!/usr/bin/env python3
import argparse
import os

def read_fasta(path):
    seqs = {}
    name = None
    chunks = []

    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())

    if name is not None:
        seqs[name] = "".join(chunks)

    return seqs

def split_genes(cell):
    if not cell.strip():
        return []
    return [x.strip() for x in cell.split(",") if x.strip()]

def read_annotation(path):
    ann = {}

    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {x: i for i, x in enumerate(header)}

        if "new_id" not in col:
            raise RuntimeError("annotation file must contain column: new_id")

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")

            if len(parts) < len(header):
                parts += ["NA"] * (len(header) - len(parts))

            new_id = parts[col["new_id"]]
            ann[new_id] = parts

    return header, ann

def read_seed_summary(paths):
    summary = {}

    for path in paths:
        if not path:
            continue
        if not os.path.exists(path):
            continue

        with open(path) as f:
            first = f.readline()

            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) < 2:
                    continue

                # format:
                # category Orthogroup seed_count seed_species seed_names seed_ids
                og = parts[1]
                summary.setdefault(og, []).append(line)

    return summary

def read_selected_info(path):
    info = {}

    with open(path) as f:
        header = f.readline().rstrip("\n")

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            og = parts[0]
            info[og] = line

    return header, info

def main():
    parser = argparse.ArgumentParser(
        description="Batch extract FASTA, IDs, annotation and metadata for selected orthogroups."
    )
    parser.add_argument("--orthogroups", required=True, help="Orthogroups.tsv from OrthoFinder")
    parser.add_argument("--og_list", required=True, help="Selected orthogroup list")
    parser.add_argument("--fasta", required=True, help="All species FAA")
    parser.add_argument("--annotation", required=True, help="all_species.annotation.map.tsv")
    parser.add_argument("--selected_info", required=True, help="selected_orthogroups.info.tsv")
    parser.add_argument("--ecm_seed_summary", default=None)
    parser.add_argument("--muscle_seed_summary", default=None)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    selected = [x.strip() for x in open(args.og_list) if x.strip()]
    selected_set = set(selected)

    seqs = read_fasta(args.fasta)
    ann_header, ann = read_annotation(args.annotation)
    info_header, selected_info = read_selected_info(args.selected_info)

    seed_summary = read_seed_summary([
        args.ecm_seed_summary,
        args.muscle_seed_summary
    ])

    og_to_members = {}

    with open(args.orthogroups) as f:
        header = f.readline().rstrip("\n").split("\t")
        species_names = header[1:]

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) < len(header):
                fields += [""] * (len(header) - len(fields))

            og = fields[0]
            if og not in selected_set:
                continue

            members = []
            species_map = []

            for sp, cell in zip(species_names, fields[1:]):
                genes = split_genes(cell)
                for g in genes:
                    members.append(g)
                    species_map.append((g, sp))

            og_to_members[og] = {
                "members": members,
                "species_map": species_map
            }

    summary_path = os.path.join(args.outdir, "batch_extract.summary.tsv")

    with open(summary_path, "w") as summary:
        summary.write("Orthogroup\tmember_count\tmissing_sequence_count\tmissing_annotation_count\toutdir\n")

        for og in selected:
            if og not in og_to_members:
                continue

            members = og_to_members[og]["members"]
            species_map = og_to_members[og]["species_map"]

            og_dir = os.path.join(args.outdir, og)
            os.makedirs(og_dir, exist_ok=True)

            faa_out = os.path.join(og_dir, f"{og}.faa")
            ids_out = os.path.join(og_dir, f"{og}.ids")
            ann_out = os.path.join(og_dir, f"{og}.annotation.tsv")
            info_out = os.path.join(og_dir, f"{og}.info.tsv")
            seed_out = os.path.join(og_dir, f"{og}.seed_summary.tsv")
            spmap_out = os.path.join(og_dir, f"{og}.species_map.tsv")

            missing_seq = 0
            missing_ann = 0

            with open(faa_out, "w") as fa, \
                 open(ids_out, "w") as ids, \
                 open(ann_out, "w") as aout, \
                 open(info_out, "w") as iout, \
                 open(seed_out, "w") as sout, \
                 open(spmap_out, "w") as smap:

                aout.write("\t".join(ann_header) + "\n")
                smap.write("gene_id\tspecies\n")

                for gid, sp in species_map:
                    ids.write(gid + "\n")
                    smap.write(f"{gid}\t{sp}\n")

                    if gid in seqs:
                        fa.write(f">{gid}\n")
                        s = seqs[gid]
                        for i in range(0, len(s), 60):
                            fa.write(s[i:i+60] + "\n")
                    else:
                        missing_seq += 1

                    if gid in ann:
                        aout.write("\t".join(ann[gid]) + "\n")
                    else:
                        missing_ann += 1

                iout.write(info_header + "\n")
                if og in selected_info:
                    iout.write(selected_info[og] + "\n")

                sout.write("category\tOrthogroup\tseed_count\tseed_species\tseed_names\tseed_ids\n")
                for line in seed_summary.get(og, []):
                    sout.write(line + "\n")

            summary.write(
                f"{og}\t{len(members)}\t{missing_seq}\t{missing_ann}\t{og_dir}\n"
            )

    print(f"Selected OGs in list: {len(selected)}")
    print(f"Extracted OGs: {len(og_to_members)}")
    print(f"Summary: {summary_path}")

if __name__ == "__main__":
    main()
