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
to this repository. The archived run log records **OrthoFinder 3.0.1b1** and the
following options. Install that historical executable separately if rerunning
the upstream inference; the downstream environment recipe does not silently
substitute version 2.5.5 or a current release.

```bash
orthofinder -f "$ECM_INPUT_ROOT/05.OrthoFinder_input" -t 32 -a 32 -M msa
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

## Matched annotation benchmarks and motif checks

The revision workers accept explicit paths and never require remote-access
credentials. All three protein files are passed directly to BUSCO without
normalization. The worker obtains BUSCO 5.8.2 and the fixed 2024-01-08
eukaryota_odb10 dataset from their official sources, records hashes, and runs
protein mode with eight CPUs, offline operation and run-statistics opt-out.

```bash
python scripts/publication_busco_matched.py \
  --project-root "$ANNOTATION_ROOT" --output "$NEW_BUSCO_OUTPUT" \
  --hmmer-bin "$HMMER_BIN_DIR" --cpu 8
python scripts/publication_original_profiles.py \
  --input "$PUBLISHED_PROTEINS" --profiles "$TARGET_PROFILES" \
  --core-profiles "$CORE_PROFILES" --hmmer "$HMMSEARCH" \
  --output "$NEW_PUBLISHED_PROFILE_OUTPUT"
python scripts/publication_revision_analysis.py \
  --package-root "$EVIDENCE_PACKAGE_ROOT" --revision-root "$REVISION_ROOT"
```

The profile worker adds the published 5,012-protein resource using the same
GA/E-value rules and core exhaustive searches as the current prediction sets.
The summarizer joins per-marker BUSCO statuses and maps all archived focal MEME
occurrences to original GA-accepted Pfam intervals. Predefined metzincin/astacin
and NPxY string checks are separate from learned MEME motifs; they are not
activity, localization or receptor-signaling assays. Input contracts and
annotation-source interpretation are in `docs/ANNOTATION.md`.
