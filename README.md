# nf-frip

`nf-frip` calculates FRiP: the fraction of mapped reads that overlap peak regions.

This module is optional. It can use IDR peaks, consensus peaks, or a custom FRiP samplesheet.

## Inputs

Auto mode:

```bash
--samples_master /path/to/samples_master.csv
--chipfilter_output /path/to/chipfilter_output
```

Select peak sources:

```bash
--frip_peak_sources idr,consensus_q0.01,consensus_q0.05
```

Provide only the directories required by the selected sources:

- `idr` needs `--idr_output`
- consensus sources need `--peak_consensus_output`

Explicit mode:

```csv
sample_id,bam,peaks,peak_set
WT_rep1,/path/to/sample.bam,/path/to/peaks.bed,custom
```

## Output

```text
${project_folder}/${frip_output}/sample.peak_set.frip.tsv
```

## Run

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --chipfilter_output /path/to/chipfilter_output \
  --idr_output /path/to/idr_output \
  --peak_consensus_output /path/to/peak_consensus_output \
  --project_folder /path/to/output_project
```

Actual execution should be tested where Nextflow is installed.
