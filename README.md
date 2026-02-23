# nf-frip

Compute FRiP (Fraction of Reads in Peaks) for ChIP-seq samples.

FRiP = reads overlapping peaks / total mapped reads.

## Input Modes (Priority Order)

1. `--frip_samplesheet` (explicit; if file exists)
2. `--samples_master` auto mode

## Mode 1: Explicit FRiP samplesheet

Provide `--frip_samplesheet` CSV with columns:

```text
sample_id,bam,peaks
```

## Mode 2: Auto from `samples_master`

Required columns:
```text
sample_id,condition
```
Optional columns used:
```text
is_control,enabled
```

Auto behavior:
- resolves BAM from `${chipfilter_output}/${sample_id}*.clean.bam`
- resolves peak file from IDR output:
  - default (`frip_peak_mode=condition`): `${condition}_idr.sorted.chr.narrowPeak`
  - optional (`frip_peak_mode=sample`): `${sample_id}_idr.sorted.chr.narrowPeak`
- by default excludes control samples (`frip_include_controls=false`)

## Output

Output directory: `${project_folder}/${frip_output}`

Per sample:
- `<sample>.frip.tsv`

Columns:
- `sample`, `bam`, `peaks`, `in_peaks`, `total_mapped`, `FRiP`

## Run

Explicit sheet:
```bash
nextflow run main.nf -profile hpc --frip_samplesheet frip_samplesheet.csv
```

Auto from `samples_master`:
```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/nf-chipfilter/chipfilter_output \
  --idr_output /path/to/nf-idr/idr_output
```

Resume:
```bash
nextflow run main.nf -profile hpc -resume
```
