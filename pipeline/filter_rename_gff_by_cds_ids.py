#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict

def read_fasta_ids(fasta):
    ids = []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                seq_id = line[1:].strip().split()[0]
                ids.append(seq_id)
    return ids

def parse_attr(attr_str):
    attrs = {}
    if attr_str == ".":
        return attrs
    for item in attr_str.strip().split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs

def format_attr(attrs):
    return ";".join([f"{k}={v}" for k, v in attrs.items()])

def get_id(attrs):
    return attrs.get("ID")

def get_parents(attrs):
    if "Parent" not in attrs:
        return []
    return attrs["Parent"].split(",")

def main():
    parser = argparse.ArgumentParser(
        description="Filter GFF3 by CDS fasta IDs, keep corresponding transcripts, rename gene/mRNA IDs."
    )
    parser.add_argument("--gff", required=True, help="Input longest isoform GFF3")
    parser.add_argument("--cds", required=True, help="Strict complete CDS fasta")
    parser.add_argument("--prefix", default="D_japonicum", help="New gene prefix, default D_japonicum")
    parser.add_argument("--out_gff", required=True, help="Output renamed GFF3")
    parser.add_argument("--map", required=True, help="Old-to-new ID map")
    args = parser.parse_args()

    keep_transcript_ids = set(read_fasta_ids(args.cds))

    features = []
    id_to_index = {}
    parent_to_children = defaultdict(list)
    child_to_parents = defaultdict(list)

    comments = []

    with open(args.gff) as f:
        for line in f:
            line = line.rstrip("\n")

            if not line:
                continue

            if line.startswith("#"):
                comments.append(line)
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            attrs = parse_attr(parts[8])
            fid = get_id(attrs)
            parents = get_parents(attrs)

            idx = len(features)
            feature = {
                "parts": parts,
                "attrs": attrs,
                "id": fid,
                "parents": parents,
                "type": parts[2],
                "index": idx
            }

            features.append(feature)

            if fid:
                id_to_index[fid] = idx

            for p in parents:
                parent_to_children[p].append(fid if fid else f"__noid__{idx}")
                child_to_parents[fid if fid else f"__noid__{idx}"].append(p)

    # 1. 确认 CDS fasta 里的 ID 是否能在 GFF3 的 mRNA/transcript 里找到
    gff_ids = set(id_to_index.keys())
    missing = sorted([x for x in keep_transcript_ids if x not in gff_ids])

    if missing:
        print("WARNING: Some CDS IDs were not found as GFF feature IDs.")
        print("First 20 missing IDs:")
        for x in missing[:20]:
            print(x)
        print("If many IDs are missing, check whether CDS fasta headers match mRNA IDs in GFF3.")

    # 2. 向上找到 gene，向下找到 exon/CDS/UTR 等子特征
    keep_ids = set()

    for tid in keep_transcript_ids:
        if tid not in id_to_index:
            continue

        # keep transcript itself
        keep_ids.add(tid)

        # upward: keep parents, usually gene
        stack = list(features[id_to_index[tid]]["parents"])
        while stack:
            pid = stack.pop()
            if pid in keep_ids:
                continue
            keep_ids.add(pid)
            if pid in id_to_index:
                stack.extend(features[id_to_index[pid]]["parents"])

        # downward: keep children, usually exon/CDS/UTR
        stack = [tid]
        while stack:
            pid = stack.pop()
            for child in parent_to_children.get(pid, []):
                keep_ids.add(child)
                if child in id_to_index:
                    stack.append(child)

    # 3. 按 GFF 中 mRNA 出现顺序分配新名字
    transcript_old_to_new = {}
    gene_old_to_new = {}

    count = 0
    for ft in features:
        if ft["id"] in keep_transcript_ids and ft["id"] in keep_ids:
            count += 1
            new_id = f"{args.prefix}{count:05d}"
            transcript_old_to_new[ft["id"]] = new_id

            # 找它的 gene parent
            for p in ft["parents"]:
                if p in keep_ids:
                    gene_old_to_new[p] = new_id + ".gene"

    # 4. 输出 map
    with open(args.map, "w") as mapout:
        mapout.write("old_id\tnew_id\tfeature_level\n")
        for old, new in transcript_old_to_new.items():
            mapout.write(f"{old}\t{new}\tmRNA\n")
        for old, new in gene_old_to_new.items():
            mapout.write(f"{old}\t{new}\tgene\n")

    # 5. 输出新的 GFF3
    feature_count_by_parent_type = defaultdict(lambda: defaultdict(int))

    with open(args.out_gff, "w") as out:
        out.write("##gff-version 3\n")

        for ft in features:
            fid = ft["id"]
            pseudo_id = fid if fid else f"__noid__{ft['index']}"

            if pseudo_id not in keep_ids and fid not in keep_ids:
                continue

            parts = ft["parts"][:]
            attrs = dict(ft["attrs"])
            ftype = parts[2]

            old_id = attrs.get("ID")
            old_parents = get_parents(attrs)

            # gene 层
            if ftype == "gene":
                if old_id in gene_old_to_new:
                    new_gene_id = gene_old_to_new[old_id]
                    base = new_gene_id.replace(".gene", "")
                    attrs = {
                        "ID": new_gene_id,
                        "Name": base
                    }
                else:
                    continue

            # mRNA/transcript 层
            elif old_id in transcript_old_to_new:
                new_transcript_id = transcript_old_to_new[old_id]
                new_attrs = {
                    "ID": new_transcript_id,
                    "Name": new_transcript_id
                }

                new_parents = []
                for p in old_parents:
                    if p in gene_old_to_new:
                        new_parents.append(gene_old_to_new[p])

                if new_parents:
                    new_attrs["Parent"] = ",".join(new_parents)

                attrs = new_attrs

            # exon/CDS/UTR 等子特征
            else:
                new_parents = []

                for p in old_parents:
                    if p in transcript_old_to_new:
                        new_parents.append(transcript_old_to_new[p])
                    elif p in gene_old_to_new:
                        new_parents.append(gene_old_to_new[p])

                if not new_parents:
                    continue

                parent_id = new_parents[0]
                feature_count_by_parent_type[parent_id][ftype] += 1
                idx = feature_count_by_parent_type[parent_id][ftype]

                new_child_id = f"{parent_id}.{ftype}{idx}"

                attrs = {
                    "ID": new_child_id,
                    "Parent": ",".join(new_parents)
                }

                if ftype == "CDS":
                    attrs["phase"] = parts[7]

            parts[8] = format_attr(attrs)
            out.write("\t".join(parts) + "\n")

    print("Done.")
    print(f"Input kept transcript IDs: {len(keep_transcript_ids)}")
    print(f"Renamed transcripts: {len(transcript_old_to_new)}")
    print(f"Output GFF3: {args.out_gff}")
    print(f"ID map: {args.map}")

if __name__ == "__main__":
    main()
