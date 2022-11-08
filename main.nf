// Enable DSL 2 syntax
nextflow.enable.dsl = 2

log.info """\
S O U P E R - S T A R
================================
input_folder    : $params.input_folder

"""

// Import processes
include { 
    sam_to_bam;
    add_tags
} from './processes.nf'

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

    // Add unique tags for each input file
    add_tags(input_ch)

}