#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def frip_output = params.frip_output ?: 'frip_output'

process compute_frip {
  tag "${sample_id}"
  stageInMode 'symlink'
  stageOutMode 'move'

  publishDir "${params.project_folder}/${frip_output}", mode: 'copy'

  input:
    tuple val(sample_id), path(bam), path(peaks)

  output:
    path("${sample_id}.frip.tsv")

  script:
  """
  set -euo pipefail

  bam_ref=\$(samtools idxstats ${bam} | awk '\$1 != "*" {print \$1; exit}')
  peak_ref=\$(awk 'NR==1{print \$1; exit}' ${peaks})
  peaks_norm="peaks.norm.bed"

  # Harmonize contig naming between BAM and peak file (chr1 vs 1)
  if [[ "\$bam_ref" == chr* && "\$peak_ref" != chr* ]]; then
    awk 'BEGIN{OFS="\\t"} {if(\$1 !~ /^chr/) \$1="chr"\$1; print}' ${peaks} > "\$peaks_norm"
  elif [[ "\$bam_ref" != chr* && "\$peak_ref" == chr* ]]; then
    awk 'BEGIN{OFS="\\t"} {sub(/^chr/,"",\$1); print}' ${peaks} > "\$peaks_norm"
  else
    cp ${peaks} "\$peaks_norm"
  fi

  total_mapped=\$(samtools view -c -F 260 ${bam})
  in_peaks=\$(bedtools intersect -u -abam ${bam} -b "\$peaks_norm" | samtools view -c -)

  frip=0
  if [[ "\$total_mapped" -gt 0 ]]; then
    frip=\$(awk -v a="\$in_peaks" -v b="\$total_mapped" 'BEGIN{printf "%.6f", a/b}')
  fi

  printf "sample\tbam\tpeaks\tin_peaks\ttotal_mapped\tFRiP\n%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample_id}" "${bam}" "${peaks}" "\$in_peaks" "\$total_mapped" "\$frip" > ${sample_id}.frip.tsv
  """
}

workflow {
  if (!params.frip_samplesheet) {
    exit 1, "ERROR: Please provide --frip_samplesheet with columns: sample_id,bam,peaks"
  }

  def input_ch = Channel
    .fromPath(params.frip_samplesheet, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
      assert row.sample_id && row.bam && row.peaks : "frip_samplesheet must contain: sample_id,bam,peaks"
      def sid = row.sample_id.toString().trim()
      def b = file(row.bam.toString())
      def p = file(row.peaks.toString())
      assert b.exists() : "BAM not found for ${sid}: ${b}"
      assert p.exists() : "Peak file not found for ${sid}: ${p}"
      tuple(sid, b, p)
    }

  compute_frip(input_ch)
}
