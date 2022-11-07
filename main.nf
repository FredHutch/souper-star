// Enable DSL 2 syntax
nextflow.enable.dsl = 2

log.info """\
S O U P E R - S T A R
================================
input_folder    : $params.input_folder

"""

workflow {

    // Get the input BAM/SAM files
    Channel
        .fromPath(
            "${params.input_folder}/*.bam",
            "${params.input_folder}/*.sam"
        )
        .set { input_ch }

    input_ch.view()
}