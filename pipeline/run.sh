#!/usr/bin/env bash
set -euo pipefail
mkdir -p 02.longest_cds 02.longest_map

for cds in 01.main_cds/*.cds.fa
do
    sp=$(basename "$cds" .cds.fa)

    python keep_longest_cds_by_gene.py \
      -i "$cds" \
      -o 02.longest_cds/${sp}.longest.cds.fa \
      --map 02.longest_map/${sp}.longest.map.tsv \
      --prefer_complete
done
seqkit stats 01.main_cds/*.cds.fa 02.longest_cds/*.longest.cds.fa \
  > 02.longest_cds.stats.txt
mkdir -p 03.strict_filter

for cds in 02.longest_cds/*.longest.cds.fa
do
    sp=$(basename "$cds" .longest.cds.fa)

    python filter_cds_strict_and_translate.py \
      --cds "$cds" \
      --out_prefix 03.strict_filter/${sp} \
      --min_aa 30
done
seqkit stats 03.strict_filter/*.strict.complete.cds.fa \
  > 03.strict_filter.cds.stats.txt

seqkit stats 03.strict_filter/*.strict.complete.faa \
  > 03.strict_filter.faa.stats.txt
mkdir -p 03.strict_filter_reason

for report in 03.strict_filter/*.strict_cds_filter.report.tsv
do
    sp=$(basename "$report" .strict_cds_filter.report.tsv)

    awk -F '\t' '$2=="FAIL"{print $3}' "$report" \
      | tr ';' '\n' \
      | sort \
      | uniq -c \
      | sort -k1,1nr \
      > 03.strict_filter_reason/${sp}.reason.summary.txt
done
mkdir -p 04.renamed_cds 04.renamed_faa 04.annotation_map

while read sp prefix
do
    echo "Processing $sp"

    python rename_strict_cds_faa_with_product.py \
      --cds 03.strict_filter/${sp}.strict.complete.cds.fa \
      --faa 03.strict_filter/${sp}.strict.complete.faa \
      --longest_map 02.longest_map/${sp}.longest.map.tsv \
      --species ${sp} \
      --prefix ${prefix} \
      --out_cds 04.renamed_cds/${sp}.cds.fa \
      --out_faa 04.renamed_faa/${sp}.faa \
      --out_map 04.annotation_map/${sp}.annotation.map.tsv

done < species_prefix.tsv
head -n 1 04.annotation_map/Dicyema_japonicum.annotation.map.tsv \
  > 04.annotation_map/all_species.annotation.map.tsv

for f in 04.annotation_map/*.annotation.map.tsv
do
    tail -n +2 "$f"
done >> 04.annotation_map/all_species.annotation.map.tsv
