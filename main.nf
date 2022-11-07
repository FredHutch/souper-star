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

    // Define channel of BAM inputs
    sam_to_bam
        .out
        .mix(bam_ch)
        .set { input_ch }

    // Add a different numeric index to each file
    ix_ch = input_ch
        .reduce( [] ) { list_of_bams, bam -> println "list_of_bams: $list_of_bams bam: $bam"; return list_of_bams.add([bam, list_of_bams.size()]) }
        .flatten()

    input_ch.view()
    ix_ch.view()
}