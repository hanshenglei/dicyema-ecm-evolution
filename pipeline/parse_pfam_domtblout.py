#!/usr/bin/env python3
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as f, open(outfile, "w") as out:
    out.write(
        "gene_id\tdomain_acc\tdomain_name\tevalue\tscore\t"
        "ali_start\tali_end\tenv_start\tenv_end\tdescription\n"
    )

    for line in f:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split(maxsplit=22)
        if len(parts) < 23:
            continue

        domain_name = parts[0]
        domain_acc = parts[1]
        gene_id = parts[3]
        evalue = parts[6]
        score = parts[7]
        ali_start = parts[17]
        ali_end = parts[18]
        env_start = parts[19]
        env_end = parts[20]
        description = parts[22]

        out.write(
            f"{gene_id}\t{domain_acc}\t{domain_name}\t{evalue}\t{score}\t"
            f"{ali_start}\t{ali_end}\t{env_start}\t{env_end}\t{description}\n"
        )
