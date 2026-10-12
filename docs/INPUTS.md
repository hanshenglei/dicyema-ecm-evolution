# Input contracts

`ECM_INPUT_ROOT` is the original numbered ECM input tree, staged outside outputs.
The domain validation worker requires:

```text
inputs/ecm/
  Dicyema_japonicum/
    Dicyema_japonicum.faa
    Dicyema_japonicum.strict.complete.faa
    Dicyema_japonicum.strict_cds_filter.report.tsv
    Dicyema_japonicum.genome.fa
  05.OrthoFinder_input/
    C_gigas.fa C_teleta.fa D_japonicum.fa H_robusta.fa I_linei.fa
    L_gigantea.fa O_sinensis.fa O_vulgaris.fa P_canaliculata.fa S_mansoni.fa
    OrthoFinder/Results_Jun11/Orthogroups/Orthogroups.tsv
  08.ECM_Muscle_annotation/pfam/ECM_Muscle.with_Dj.pfam.domtblout
  10.family_analysis/selected_orthogroups.info.tsv
  12.family_groups/OG_family_summary.tsv
```

`ECM_PFAM_HMM` points to a locally obtained Pfam-A HMM database with the exact
40 accession-version profiles in `config/profile_metadata.tsv`. The worker
records profile hashes, gathering scores, source hashes and exact commands.
The original local database had no explicit release label; profile versions
and hashes, not a guessed Pfam release, define the model provenance.

For coordinate and comparative evidence scripts, retain the following mirror
under `ECM_PROJECT_ROOT`:

```text
remote_results_cache/                 # original small table/XML result mirror
  05.OrthoFinder_input/OrthoFinder/Results_Jun11/
    Orthogroups/Orthogroups.GeneCount.tsv
    Orthologues/D_japonicum.tsv
    Gene_Duplication_Events/Duplications.tsv
  07.ECM_Muscle_candidates/D_japonicum.ECM.related.og_map.tsv
  08.ECM_Muscle_annotation/pfam/ECM_Muscle.with_Dj.pfam.domain.tsv
  09.presence_absence/ECM.orthogroup.presence_absence.tsv
  09.presence_absence/Muscle.orthogroup.presence_absence.tsv
  10.family_analysis/selected_OGs/OG*/OG*.motif.tsv
  12.family_groups/OG_family_summary.tsv
publication_package/remote_sources/   # original sequence/GFF/tree mirror
  05.OrthoFinder_input/D_japonicum.fa
  Dicyema_japonicum/Dicyema_japonicum.strict.complete.renamed.faa
  Dicyema_japonicum/Dicyema_japonicum.strict.complete.renamed.gff3
  08.ECM_Muscle_annotation/pfam/ECM_Muscle.with_Dj.pfam.domtblout
  10.family_analysis/selected_OGs/OG*/  # original FAA, trimmed alignment, tree
```

Keep the OrthoFinder snapshot name `Results_Jun11` for this archived run, or
adapt the `ORTH`/source path constants when analyzing another run. Never join
renamed models by an unverified identifier prefix: coordinate analysis checks
exact sequences before linking the GFF. The scripts create their result folders
under `publication_package`. Genome-check outputs deliberately refuse overwrite.

The supplied comparative sampling has ten species, including two octopus
comparators, and no human proteome. The workflow validates a fixed archived
study; assertions about dataset counts are intentional, not general API limits.
