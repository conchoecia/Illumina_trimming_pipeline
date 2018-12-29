# Illumina trimming pipeline

This trims libraries from Illumina sequencing.

- Features
  - Uses snakemake for parallelization and resource management.
  - Library names are automatically determined just by specifying directories.
  - Options to specify trimming with TruSeq2 or Nextera Adapters.
  - Options to trim Dovetail Chicago and HiC data. The junction sequence is user-specified.
  - Generates a report of trimming efficiency and linker content sequence.

# Requirements

- This pipeline requires
  - Snakemake (for python 3)
  - bbmerge
  - trimmomatic
  - gcc

# Instructions

There is an example called `directories.yaml(example)`. Rename your
actual file to `directories.yaml`. Run `snakemake --cores <num_cores>
-p -r` in the directory with the `directories.yaml`, `bin/`, and
`Snakemake` files.
