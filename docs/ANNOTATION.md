# Annotation input provenance

The assembly is associated with BioProject PRJNA496385 and WGS RDEW01
(GenBank GCA_011109175.1), as documented by the NCBI assembly report and Lu et al.
(2019), doi:10.1093/gbe/evz157. The integrated GFF supplied with the study is
distinct from the original published gene-model inventory.

Retained GFF source labels include AUGUSTUS, GeneMark.hmm3, gmst, transdecoder,
miniprot and EVM. Source labels establish retained-model provenance, not software
version numbers or the complete execution history. Final source composition and
coordinate-phase signature matches are provided with Supplementary Data S2.
An AGAT v1.4.1 standardization log is present in the archived source records.
GffCompare v0.12.6 records a standardized AUGUSTUS reference and four input
tracks (BRAKER, EVM, miniprot and TransDecoder). A CD-HIT-format intermediate
clustering file is not evidence that the final resource is locus-deduplicated.
Helixer v0.3.5 and EGAPx explorations are documented separately from final
retained-model provenance. Annotation-stage DIAMOND logs record v2.0.11.149,
whereas currently installed executables must not be assigned retrospectively.

`pipeline/filter_rename_gff_by_cds_ids.py` is the recovered focal GFF filtering
and renaming helper. It connects passing CDS identifiers and mapping records to
the coordinate models used in the downstream overlap analysis. Keep the
exact-input checksums and identifier maps supplied with the study.

This release reproduces the downstream ECM analysis from the supplied integrated
annotation; it does not claim to regenerate the annotation from scratch. A
historical weighted-integration command was not recovered and is not fabricated
as a rerun command. Archived BUSCO outputs refer to their explicitly recorded
protein inputs, which must not silently be substituted for another export.

The matched revision benchmark uses BUSCO v5.8.2, HMMER v3.4 and the same
eukaryota_odb10 dataset (2024-01-08, 255 markers) for the published 5,012 proteins,
the actual full 17,260 predictions and the strict 12,402 predictions. Raw input
hashes, complete/fragmented/missing marker transitions and duplicated-marker
counts accompany the manuscript. Greater marker recovery is not automatically
greater independent-locus precision or evidence for biological duplication.
