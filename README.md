# nf-frip

Compute FRiP (Fraction of Reads in Peaks) for ChIP-seq samples.

FRiP = reads overlapping peaks / total mapped reads.

## Input

Provide `--frip_samplesheet` CSV with columns:

```text
sample_id,bam,peaks
```

Example:

```text
sample_id,bam,peaks
GAR0968,/path/to/GAR0968.clean.bam,/path/to/GAR0968_idr.sorted.chr.narrowPeak
GAR0979,/path/to/GAR0979.clean.bam,/path/to/GAR0979_idr.sorted.chr.narrowPeak
```

## Output

Output directory: `${project_folder}/${frip_output}`

Per sample:

- `<sample>.frip.tsv`

Columns:

- `sample`, `bam`, `peaks`, `in_peaks`, `total_mapped`, `FRiP`

## Run

```bash
nextflow run main.nf -profile hpc --frip_samplesheet frip_samplesheet.csv
```

Resume:

```bash
nextflow run main.nf -profile hpc --frip_samplesheet frip_samplesheet.csv -resume
```
