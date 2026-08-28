# Data availability

No sequencing reads, alignments, genome references, indexes, or large genomic
tables are redistributed in this repository.

## Public accessions

The Mature versus Memory comparison used:

- Mature NK: SRR7650764, SRR7650846, SRR7650883, SRR7650911
- Memory NK: SRR7650766, SRR7650808, SRR7650848, SRR7650885, SRR7650913,
  SRR7650922

The separate Week 1 training sample was SRR7650763 (Immature NK).

These accessions are associated with Calderon et al. (2019),
https://doi.org/10.1038/s41588-019-0505-9. Retrieve public data from the
appropriate NCBI archive and verify current archive metadata before reuse.

## Inputs not included

Reproduction requires independently obtained and validated inputs:

- paired FASTQ files for the Week 1 teaching exercise;
- the ten comparison BAM/BAI pairs and sample-level broadPeak files, or a
  documented method to regenerate equivalent inputs;
- GRCh38 reference resources;
- an hg38 blacklist with compatible chromosome naming;
- GENCODE Release 50 GRCh38 primary-assembly GTF and FASTA resources;
- locally installed analysis software listed in `environment/software_versions.tsv`.

Local filesystem paths from the internship environment were deliberately
removed. `config/input_manifest.example.tsv` shows the intended public format.
