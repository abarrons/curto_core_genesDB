# Curto Core Genes DB

Pipeline to build a custom, taxonomically-labeled **nucleotide BLAST database of
bacterial core-gene marker sequences**, intended for use as a reference database
to classify metagenomic reads (e.g. via BLAST against the resulting database).

> **Note:** This repository documents a research pipeline originally built for
> execution on a SLURM HPC cluster (UCI HPC3, account `JMARTINY_LAB`). Scripts
> contain hardcoded absolute paths (e.g. `/dfs5/bio/abarrons/...`) that will
> need to be updated for any other environment.

## Pipeline overview

The scripts in [`process_executables/`](process_executables) are numbered in the
order they should be run:

```
genome assemblies
      │
      ▼
[1-2] copy assemblies, filter scaffolds, rename to isolate IDs
      │
      ▼
[3] Prodigal (ORF prediction) → hmmsearch vs. core-gene HMM profiles
    → extract best-hit nucleotide sequence per gene per genome
      │   (automatically submits step 4)
      ▼
[4] CD-HIT dereplication per gene + prevalence summary table
      │
      ▼
[5] add NCBI taxid-tagged headers, concatenate, makeblastdb
      │
      ▼
 referenceDB/total_coregenes.fna  (BLAST nucleotide database)
```

| Step | Script | Purpose |
|------|--------|---------|
| 1 | [`1_move_og_genomes_blobplots.sh`](process_executables/1_move_og_genomes_blobplots.sh) | Copy mapped assemblies + BlobToolKit outputs into `og_genomes/` / `blobtools_out/` using the batch mapping files (`og_isolates_b1_md.txt`, `og_isolates_b2_md.txt`). |
| 2 | [`2_filterbyGCandLenght.sh`](process_executables/2_filterbyGCandLenght.sh) | Filter scaffolds using BlobTools hit-table columns (length/coverage), then rename filtered assemblies from sequencing ID to isolate ID. |
| 3 | [`3_prodigal2hmmscan.sh`](process_executables/3_prodigal2hmmscan.sh) | Run Prodigal, `hmmsearch` each protein against the core-gene HMM profiles, and extract the best-hit nucleotide sequence per gene/genome into `all_genomes4DB/filtered_seqs/`. The last array task concatenates per-gene files and submits step 4. |
| 4 | [`4_derep_seqs.sh`](process_executables/4_derep_seqs.sh) | Dereplicate each gene's sequences with CD-HIT (100% identity, including sequences contained in longer ones) into `clean_genes/`, and write `final_summary_table.txt` with per-gene prevalence/uniqueness stats. |
| 5 | [`5_build_blastDB.sh`](process_executables/5_build_blastDB.sh) | Map genomes to NCBI taxids, rename sequence headers (`gnl\|BACT\|...`), concatenate all genes, and run `makeblastdb`. |
| 6 | [`6_count_genes_by_genome.sh`](process_executables/6_count_genes_by_genome.sh) | Optional summary: count how many core genes were recovered per genome. |
| 7 | [`7_calculate_gene_lengths.py`](process_executables/7_calculate_gene_lengths.py) | Optional summary: build a gene-length matrix (genes x genomes). |

Supporting files (not numbered steps):

- [`hmmpress_optional.sh`](process_executables/hmmpress_optional.sh): `hmmpress` the HMM profiles. Not needed by step 3, which uses `hmmsearch` on the profiles directly.
- [`get_taxIDs.py`](process_executables/get_taxIDs.py): one-off NCBI Entrez query used to create `ncbi_taxids.txt`.
- [`curto_core_genes_analysis.R`](curto_core_genes_analysis.R): reads step 4's `final_summary_table.txt`, selects genes passing a uniqueness threshold (`DBgenes.txt`), and plots prevalence/uniqueness.

### Known gaps / things to verify before re-running end-to-end

- **Gene-list selection isn't wired through**: `curto_core_genes_analysis.R`
  writes a filtered `DBgenes.txt`, but steps 4 and 5 iterate over *all* `.hmm`
  profiles rather than reading it. Decide whether all profiled genes or only
  the curated subset should end up in the database.
- Step 3 is named for `hmmscan` but actually calls `hmmsearch`.

## Legacy scripts

[`process_executables/legacy/`](process_executables/legacy) holds scripts that are no longer part of the pipeline:

- `analyze_core_genes_duplicates.sh` + `find_duplicates.py`: former step 4. Only
  flagged exact-duplicate sequences (reports never used downstream) and moved
  `filtered_seqs` into `analysis_out_rd/`. Superseded by CD-HIT in step 4.
- `array_test.sh`: debug copy of the former step 4.
- `blast_hmmbuild.sh`: older SGE/`qsub` workflow that builds a *protein* BLAST database from a different directory layout.
- `move_excluded_genomes.sh`: superseded by the exclusion-list handling in step 3.
- `remove_noisy_genomes.sh`: targets the old `analysis_out/` path and hardcodes two genome IDs.
- `reftrees.sh`: legacy RAxML tree-building script from an earlier project.

## Requirements

- SLURM (scripts contain `#SBATCH` directives; adapt or strip these to run
  outside a SLURM cluster)
- [Prodigal](https://github.com/hyattpd/Prodigal) 2.6.3
- [HMMER](http://hmmer.org/) 3.3 (`hmmsearch`, `hmmpress`)
- [BBMap](https://sourceforge.net/projects/bbmap/) (`filterbyname.sh`)
- [CD-HIT](https://github.com/weizhongli/cdhit)
- [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) 2.13.0 (`makeblastdb`, `blastdbcmd`)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- Python 3 with `biopython`, `pandas`
- R with `tidyverse` (for `curto_core_genes_analysis.R`)

## Repository layout

```
.
├── process_executables/    # numbered pipeline scripts (+ legacy/)
├── og_refseq_genomes/      # reference genome downloads (gitignored, large)
├── *.txt / *.xlsx          # metadata/mapping tables used by the pipeline
├── final_summary_table.txt # per-gene prevalence/uniqueness summary (step 4 output)
├── gene_prevalence.png,
│   gene_uniqueness.png     # plots produced by curto_core_genes_analysis.R
└── curto_core_genes_analysis.R
```

Large/generated artifacts (raw reference genomes, compiled FASTA/BLAST inputs,
intermediate pipeline directories, SLURM logs) are excluded via `.gitignore`
since they are reproducible outputs rather than source files.
