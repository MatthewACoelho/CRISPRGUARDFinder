#!/usr/bin/env nextflow

params.pam = "NGG"
params.genome = "hg38"
params.id = "unknown"
params.guide_min_pvalue = 0
params.guide_length = 20
params.guide_mismatches = 5
params.guard_length = 14
params.guard_mismatches = 3
params.out_path = "."
params.chr = ""
params.start = ""
params.end = ""
params.strand = ""
params.max_guard_distance = 0

process find_guide_off_targets {

	output:
		file "guide.txt.stderr" into off_target_detail1, off_target_detail2
		file "guide.txt.stdout" into off_target_summary1, off_target_summary2

	script:
		"""
		\$guard_root/bin/ot query -quiet -n 4 \
			-seq ${params.guide} \
			-mismatch ${params.guide_length} ${params.guide_mismatches} \
			-info \$guard_root/data/${params.genome}/info \
			-format annotated \
			\$guard_root/data/${params.genome}/${params.pam}/*.inx \
			> guide.txt.stdout 2> guide.txt.stderr
		"""
}

process find_guards {

	input:
    	file off_targets from off_target_detail1
    	file off_target_summary from off_target_summary1

    output:
		file "*_tmp.txt" into guard_tmp
		file "*_guides.txt" into guards

    script:
    	"""
		Rscript \$guard_root/bin/find_guards.R -n \
			--id ${params.id} \
			--genome ${params.genome}\
			--guard_len ${params.guard_length} \
			--guard_mismatches ${params.guard_mismatches} \
			--max_guard_distance ${params.max_guard_distance} \
			--crrna ${params.guide} \
			--chr "${params.chr}" \
			--start "${params.start}" \
			--end "${params.end}" \
			--strand "${params.strand}" \
			--crrna_min_pvalue ${params.guide_min_pvalue}
    	"""
}

process find_guard_off_targets {
	tag "$guards"

	input:
		each file(guards) from guards

	output:
		file "*.stderr" into guard_off_target_detail
		file "*.stdout" into guard_off_target_summary

	script:
		"""
		\$guard_root/bin/ot query -quiet -n 4 \
			-list $guards \
			-mismatch ${params.guard_length} ${params.guard_mismatches} \
			-info \$guard_root/data/${params.genome}/info \
			-format exact_annotated \
			\$guard_root/data/${params.genome}/${params.pam}/*.inx \
			> ${guards}.stdout 2> ${guards}.stderr
		"""
}

process score_guards {
	tag "score_guards"

    publishDir "${params.out_path}", mode: 'copy'

	input:
    	file off_targets from off_target_detail2
    	file off_target_summary from off_target_summary2
		file tmp from guard_tmp.collect()
    	file off_targets from guard_off_target_detail.collect()
    	file off_target_summary from guard_off_target_summary.collect()

    output:
		file "${params.id}_final.txt"

    script:
    	"""
		Rscript \$guard_root/bin/find_guards.R \
			--id ${params.id} \
			--genome ${params.genome} \
			--guard_len ${params.guard_length} \
			--guard_mismatches ${params.guard_mismatches} \
			--max_guard_distance ${params.max_guard_distance} \
			--crrna ${params.guide} \
			--chr "${params.chr}" \
			--start "${params.start}" \
			--end "${params.end}" \
			--strand "${params.strand}" \
			--crrna_min_pvalue ${params.guide_min_pvalue}
    	"""
}

workflow.onComplete
{
    log.info "--------------------------"
    log.info "Pipeline execution summary"
    log.info "--------------------------"
    log.info "Started at  : ${workflow.start}"
    log.info "Completed at: ${workflow.complete}"
    log.info "Duration    : ${workflow.duration}"
    log.info "Success     : ${workflow.success}"
    log.info "Work dir    : ${workflow.workDir}"
    log.info "Exit status : ${workflow.exitStatus}"
    log.info "Error report: ${workflow.errorReport ?: '-'}"
}
