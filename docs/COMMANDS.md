# Commands and analysis scope

Commands below preserve recorded options while replacing machine-specific paths
with environment variables. Thread counts do not define biological replicates.

## CDS processing

```bash
python pipeline/keep_longest_cds_by_gene.py \
  -i "$CDS_INPUT" -o "$LONGEST_CDS" --map "$LONGEST_MAP" --prefer_complete
python pipeline/filter_cds_strict_and_translate.py \
  --cds "$LONGEST_CDS" --out_prefix "$STRICT_PREFIX" --min_aa 30
```

`pipeline/run.sh` is the recovered preprocessing batch script. Run it in the
numbered input directory with its helper scripts and `species_prefix.tsv`
present. It does not run OrthoFinder or create a complete ECM annotation.

## Orthogroups, gene trees and motifs

OrthoFinder result tables, resolved gene trees and duplication events are inputs
to this repository. A fresh inference can be run as follows; this is a rerun
example, not a claim that the original invocation used every default identically.

```bash
orthofinder -f "$ECM_INPUT_ROOT/05.OrthoFinder_input" -t 16 -a 16
bash pipeline/run_one_OG_tree_motif.sh "$OG_DIR" 10 10 20000000 yes
```

The recovered per-OG helper executes MAFFT `--auto --thread`, trimAl
`-automated1`, IQ-TREE 2 `-m MFP -B 1000 -alrt 1000`, and MEME
`-protein -mod zoops -nmotifs 10 -minw 6 -maxw 50 -maxsize 20000000`.
If trimAl yields an empty output it records an untrimmed-alignment fallback.
The separate batch helper skips MEME above its configured sequence limit.
Keep those records: absent MEME output means not analyzed, not zero motifs.

## Domains and genomic rescue

The independent worker executes two acceptance schemes on twelve datasets:
HMMER `--cut_ga`, and joint sequence/independent-domain E <= 1e-5. It also runs
`--max` for fourteen core profiles in both full and strict focal proteomes.
Total: 28 searches. It records normalized scan sequences, raw outputs and hashes.
HMMER coordinates are 1-based inclusive. The optional TBLASTN command is:

```bash
tblastn -query queries.faa -subject "$GENOME" -num_threads 4 \
  -seg yes -evalue 1e-5 -max_target_seqs 100 -max_hsps 20 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp sframe' \
  -out hits.outfmt6.tsv
```

The worker selects four C4/Laminin_N domain queries from two positive comparison
species and two within-Dicyema controls; it applies a 180-second bound and records
timeouts without interpreting them as negative results.

## Coordinate sensitivity and evidence

Coordinate analysis merges overlapping/adjacent CDS intervals, then builds
same-OG, same-strand connected components using shared CDS / shorter CDS length
at 0.25, 0.50 and 0.80. Evidence extraction expands and deduplicates explicitly
exported orthology pairs. Terminal duplication descendant sets are compared
against overlap components; internal-node events are reported separately.
Seed inventories are conditional keyword-derived sets, not a genome-wide loss
denominator. Final figure scripts use the existing trees; no species phylogeny
or ancestor-loss reconstruction is added by plotting.
