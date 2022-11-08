// Enable DSL 2 syntax
nextflow.enable.dsl = 2

log.info """\
S O U P E R - S T A R
==========================================
samplesheet      : $params.samplesheet
samplesheet_sep  : $params.samplesheet_sep
umi_len          : $params.umi_len

CONTAINERS
=========================================
samtools         : $params.container__samtools
simplesam        : $params.container__simplesam
soupercell       : $params.container__soupercell

"""

// Import processes
include { 
    sam_to_bam;
    add_tags;
    merge
} from './processes.nf'

workflow {

    // Get the input BAM/SAM files from the samplesheet
    Channel
        .fromPath(
            "${params.samplesheet}",
            checkIfExists: true
        )
        .splitCsv(
            header: true,
            sep: "${params.samplesheet_sep}",
            strip: true
        )
        .map {
            it -> [
                file(
                    "${it.path}",
                    checkIfExists: true
                ),
                "${it.sample}"
            ]
        }
        .branch {
            bam: it[0].name.endsWith(".bam")
            sam: true
        }
        .set { input }

    // SAM -> BAM
    sam_to_bam(input.sam)

    // Define channel of BAM inputs
    sam_to_bam
        .out
        .mix(input.bam)
        .set { bam_ch }

    // Add unique tags for each input file
    add_tags(bam_ch)

    // Merge
    merge(add_tags.out.groupTuple())

    // Remove duplicates
    dedup(merge.out)

    // Index the BAM
    index(dedup.out)

}