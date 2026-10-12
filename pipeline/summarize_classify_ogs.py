#!/usr/bin/env python3
import os
import re
import csv
import argparse
from collections import Counter, defaultdict

BAD_NAME_PAT = re.compile(
    r'^(na|hypothetical protein|uncharacterized protein|predicted protein)$',
    re.I
)
LOW_QUAL_PAT = re.compile(
    r'low quality|hypothetical protein|uncharacterized protein|predicted protein|partial',
    re.I
)

FAMILY_RULES = [
    ("Collagens",    [r'collagen', r'type iv collagen', r'\bcol4\b']),
    ("Laminin",      [r'laminin', r'\blama\b', r'\blamb\b', r'\blamc\b', r'nidogen']),
    ("Integrin",     [r'integrin', r'\bitga\b', r'\bitgb\b']),
    ("Cadherin",     [r'cadherin', r'protocadherin']),
    ("Fibronectin",  [r'fibronectin', r'\bfn1\b']),
    ("Elastin",      [r'elastin', r'tropoelastin', r'fibrillin', r'fibulin', r'emilin']),
    ("Proteoglycan", [r'proteoglycan', r'glypican', r'syndecan', r'perlecan', r'aggrecan', r'versican', r'decorin', r'biglycan']),
    ("MMP",          [r'matrix metalloproteinase', r'\bmmp\b', r'metalloproteinase', r'adam', r'adamts']),
    ("Paxillin",     [r'paxillin', r'\bpxn\b']),
    ("Actin",        [r'\bactin\b', r'actinin', r'actin-related', r'\barp\b']),
    ("Myosin",       [r'\bmyosin\b', r'tropomyosin', r'troponin', r'titin', r'nebulin', r'desmin', r'dystrophin', r'myomesin']),
]

def read_tsv(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)

def simplify_name(name):
    if name is None:
        return "NA"
    x = name.strip()
    if not x:
        return "NA"
    x = re.sub(r'^LOW QUALITY PROTEIN:\s*', '', x, flags=re.I)
    x = re.sub(r'\s+isoform\s+X\d+', '', x, flags=re.I)
    x = re.sub(r'\s+putative$', '', x, flags=re.I)
    x = re.sub(r'\s+homologue$', '', x, flags=re.I)
    x = re.sub(r'\s+like$', '', x, flags=re.I)
    x = re.sub(r'\s+', ' ', x).strip()
    return x

def is_good_name(name):
    if name is None:
        return False
    x = simplify_name(name)
    if BAD_NAME_PAT.match(x):
        return False
    if LOW_QUAL_PAT.search(x):
        return False
    return True

def assign_family(text_blob, seed_category=""):
    t = (text_blob + " " + seed_category).lower()
    for fam, pats in FAMILY_RULES:
        for p in pats:
            if re.search(p, t):
                return fam
    if "ecm" in seed_category.lower():
        return "Other_ECM"
    if "muscle" in seed_category.lower():
        return "Other_Muscle"
    return "Ambiguous"

def count_fasta_headers(faa):
    n = 0
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root_dir", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    root_dir = args.root_dir
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    summary_tsv = os.path.join(outdir, "OG_family_summary.tsv")

    rows_out = []

    og_dirs = [os.path.join(root_dir, x) for x in os.listdir(root_dir)
               if os.path.isdir(os.path.join(root_dir, x))]

    for og_dir in sorted(og_dirs):
        og = os.path.basename(og_dir)

        ann_file = os.path.join(og_dir, f"{og}.annotation.tsv")
        domain_file = os.path.join(og_dir, f"{og}.pfam.domain.tsv")
        info_file = os.path.join(og_dir, f"{og}.info.tsv")
        faa_file = os.path.join(og_dir, f"{og}.faa")
        spmap_file = os.path.join(og_dir, f"{og}.species_map.tsv")

        if not os.path.exists(ann_file) or not os.path.exists(faa_file):
            continue

        ann = read_tsv(ann_file)
        n_tip = count_fasta_headers(faa_file)

        species = set()
        dj_copy = 0
        good_names = []
        all_names = []
        lowq_count = 0

        for r in ann:
            sp = r.get("species", "Unknown")
            gid = r.get("new_id", r.get("gene_id", "NA"))
            pname = r.get("protein_name", "NA")

            species.add(sp)
            if sp == "Dicyema_japonicum":
                dj_copy += 1

            pname2 = simplify_name(pname)
            all_names.append(pname2)

            if LOW_QUAL_PAT.search(pname2) or pname2.upper() == "NA":
                lowq_count += 1

            if is_good_name(pname2):
                good_names.append(pname2)

        seed_category = ""
        if os.path.exists(info_file):
            info = read_tsv(info_file)
            if len(info) > 0:
                seed_category = info[0].get("category", "")

        domain_counter = Counter()
        if os.path.exists(domain_file):
            dom = read_tsv(domain_file)
            seen = set()
            for r in dom:
                gid = r.get("gene_id", "")
                dname = r.get("domain_name", "").strip()
                dacc = r.get("domain_acc", "").strip()
                label = dname if dname else dacc
                key = (gid, label)
                if gid and label and key not in seen:
                    seen.add(key)
                    domain_counter[label] += 1

        top_good_name, top_good_name_n = ("NA", 0)
        if good_names:
            c = Counter(good_names)
            top_good_name, top_good_name_n = c.most_common(1)[0]

        top_domain, top_domain_n = ("NA", 0)
        if domain_counter:
            top_domain, top_domain_n = domain_counter.most_common(1)[0]

        good_name_ratio = round(len(good_names) / n_tip, 4) if n_tip else 0
        lowq_ratio = round(lowq_count / n_tip, 4) if n_tip else 0
        dom_ratio = round(top_domain_n / n_tip, 4) if n_tip else 0
        name_ratio = round(top_good_name_n / n_tip, 4) if n_tip else 0

        consensus_name = top_good_name
        if consensus_name == "NA" and top_domain != "NA":
            consensus_name = top_domain
        if consensus_name == "NA" and seed_category:
            consensus_name = seed_category

        text_blob = " ".join([top_good_name, top_domain, consensus_name] + good_names[:20])
        family = assign_family(text_blob, seed_category)

        if n_tip >= 5 and (good_name_ratio >= 0.40 or dom_ratio >= 0.50) and (name_ratio >= 0.25 or dom_ratio >= 0.50):
            quality = "HQ"
        elif n_tip >= 5 and (good_name_ratio >= 0.20 or dom_ratio >= 0.30):
            quality = "MQ"
        else:
            quality = "LQ"

        rows_out.append({
            "Orthogroup": og,
            "n_tip": n_tip,
            "species_n": len(species),
            "D_japonicum_copy": dj_copy,
            "seed_category": seed_category,
            "good_name_ratio": good_name_ratio,
            "low_quality_ratio": lowq_ratio,
            "top_good_name": top_good_name,
            "top_good_name_count": top_good_name_n,
            "top_domain": top_domain,
            "top_domain_count": top_domain_n,
            "dominant_domain_ratio": dom_ratio,
            "dominant_name_ratio": name_ratio,
            "consensus_name": consensus_name,
            "assigned_family": family,
            "quality_class": quality
        })

    # 写 summary
    fieldnames = list(rows_out[0].keys()) if rows_out else []
    with open(summary_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    # 建目录并软链接
    for r in rows_out:
        fam = r["assigned_family"]
        q = r["quality_class"]
        og = r["Orthogroup"]
        src = os.path.abspath(os.path.join(root_dir, og))
        dst_dir = os.path.join(outdir, q, fam)
        os.makedirs(dst_dir, exist_ok=True)
        dst = os.path.join(dst_dir, og)
        if not os.path.exists(dst):
            os.symlink(src, dst)

    print("Done.")
    print("Summary:", summary_tsv)

if __name__ == "__main__":
    main()
