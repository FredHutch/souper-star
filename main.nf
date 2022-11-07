// Enable DSL 2 syntax
nextflow.enable.dsl = 2

log.info """\
S O U P E R - S T A R
================================
input_folder    : $params.input_folder

"""

// Import processes
include { sam_to_bam } from './processes/samtools.nf'

workflow {

    // Get the input BAM/SAM files
    bam_ch = Channel.fromPath("${params.input_folder}/*.bam")
    sam_ch = Channel.fromPath("${params.input_folder}/*.sam")

    // SAM -> BAM
    sam_to_bam(sam_ch)

    input_ch.view()
}