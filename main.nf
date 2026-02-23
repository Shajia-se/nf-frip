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
  def input_ch

  if (params.frip_samplesheet && file(params.frip_samplesheet).exists()) {
    input_ch = Channel
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
  } else if (params.samples_master) {
    def master = file(params.samples_master)
    assert master.exists() : "samples_master not found: ${params.samples_master}"

    def header = null
    def records = []
    master.eachLine { line, n ->
      if (!line?.trim()) return
      def cols = line.split(',', -1)*.trim()
      if (n == 1) {
        header = cols
      } else {
        def rec = [:]
        header.eachWithIndex { h, i -> rec[h] = i < cols.size() ? cols[i] : '' }
        records << rec
      }
    }

    assert header : "samples_master header not found: ${params.samples_master}"
    assert header.contains('sample_id') : "samples_master missing required column: sample_id"
    assert header.contains('condition') : "samples_master missing required column: condition"

    def isEnabled = { rec ->
      def v = rec.enabled?.toString()?.trim()?.toLowerCase()
      (v == null || v == '' || v == 'true')
    }
    def isControl = { rec ->
      rec.is_control?.toString()?.trim()?.toLowerCase() == 'true'
    }

    def includeControls = (params.frip_include_controls == null) ? false : params.frip_include_controls
    def peakMode = (params.frip_peak_mode ?: 'condition').toString()
    def peakSuffix = (params.frip_peak_suffix ?: '_idr.sorted.chr.narrowPeak').toString()
    def bamDir = file(params.chipfilter_output)
    def idrDir = file(params.idr_output)
    assert bamDir.exists() : "chipfilter_output not found: ${params.chipfilter_output}"
    assert idrDir.exists() : "idr_output not found: ${params.idr_output}"

    def resolveBam = { sid ->
      def hits = bamDir.listFiles()?.findAll { f ->
        f.isFile() && f.name.endsWith('.clean.bam') && (f.name == "${sid}.clean.bam" || f.name.startsWith("${sid}_"))
      } ?: []
      if (hits.isEmpty()) throw new IllegalArgumentException("No clean BAM found for sample_id '${sid}' under: ${params.chipfilter_output}")
      if (hits.size() > 1) throw new IllegalArgumentException("Multiple clean BAM files matched sample_id '${sid}': ${hits*.name.join(', ')}")
      file(hits[0].toString())
    }

    def rows = records
      .findAll { rec -> isEnabled(rec) }
      .findAll { rec -> includeControls ? true : !isControl(rec) }
      .collect { rec ->
        def sid = rec.sample_id?.toString()?.trim()
        def cond = rec.condition?.toString()?.trim()
        if (!sid || !cond) return null
        def bam = resolveBam(sid)
        def peakName = (peakMode == 'sample') ? "${sid}${peakSuffix}" : "${cond}${peakSuffix}"
        def peaks = file("${params.idr_output}/${peakName}")
        assert peaks.exists() : "Peak file not found for sample '${sid}': ${peaks}"
        [sid: sid, bam: bam, peaks: peaks]
      }
      .findAll { it != null }

    input_ch = Channel
      .fromList(rows)
      .ifEmpty { exit 1, "ERROR: No FRiP input rows generated from samples_master: ${params.samples_master}" }
      .map { r -> tuple(r.sid, r.bam, r.peaks) }
  } else {
    exit 1, "ERROR: Provide --frip_samplesheet (existing file) or --samples_master for auto FRiP input."
  }

  compute_frip(input_ch)
}
