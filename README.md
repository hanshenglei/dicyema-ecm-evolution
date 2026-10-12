# Dicyema ECM evolution

Code for comparative orthogroup, Pfam-domain, MEME-motif and genomic-coordinate
analysis in *Dicyema japonicum*. The public repository contains analysis code,
commands and configuration only. It does not contain unpublished manuscripts,
sequence inputs, Pfam databases or server-access credentials.

## Organization

- `pipeline/`: recovered preprocessing and original family-analysis helpers.
- `scripts/`: independent domain validation, model-overlap sensitivity,
  orthology/duplication evidence extraction and manuscript-figure generation.
- `config/`: species prefixes, environment example and profile metadata.
- `docs/`: input contracts, recorded parameters and code provenance.
- `tests/`: small synthetic unit tests without biological sequence inputs.

## Run

Use a Linux environment for the external bioinformatics tools. The environment
recipe is a reconstruction recipe, not a byte-identical historical lockfile.
HMMER 3.4 and BLAST+ 2.17.0 were used in the independent validation. Preserve the
versions and hashes recorded with each original analysis when reproducing it.

```bash
conda env create -f environment.yml
conda activate dicyema-ecm
source config/project.env.example
# Stage the inputs described in docs/INPUTS.md, then:
python scripts/publication_domain_validation.py
python scripts/publication_domain_supplement.py
python scripts/publication_domain_genome_check.py
python scripts/publication_domain_verify.py
python scripts/publication_coordinates.py
python scripts/publication_analyze.py
python scripts/publication_figures.py
python -m unittest discover -s tests -v
```

Set `ECM_PROJECT_ROOT` to the staged workspace. The domain worker will not
overwrite an existing run manifest and will not write into its input tree.
Genome validation is a bounded six-query analysis, not genome-wide reannotation.
For upstream commands and their scope, see `docs/COMMANDS.md`.

## Interpretation

An orthogroup can contain both orthologs and paralogs. Use exported pair records
and species-tree event labels rather than OG membership alone. Protein-model
counts and overlap components are different units. HQ/MQ/LQ measure annotation
and domain agreement; they are not phylogenetic branch confidence. A missing
profile hit is detection evidence under the recorded thresholds, not a proof of
genomic deletion. Motif IDs are local to individual MEME runs.

Scientific input and result tables accompany the manuscript, separately from
this code repository. Reproduction needs those files; a code-only clone is not
a complete dataset. Public source assembly: Lu et al. (2019),
doi:10.1093/gbe/evz157, BioProject PRJNA496385, WGS RDEW00000000.
The study's integrated annotation is distinct from the original published gene
models. Do not silently replace one with the other.

## License

Project code is released under the MIT license. External programs, databases
and independently obtained inputs retain their own licenses.
