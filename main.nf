// Enable DSL 2 syntax
nextflow.enable.dsl = 2

log.info """\
S O U P E R - S T A R
==========================================
samplesheet      : $params.samplesheet
samplesheet_sep  : $params.samplesheet_sep
umi_len          : $params.umi_len
genome_fasta     : $params.genome_fasta
k                : $params.k
flags            : $params.flags
results          : $params.results

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
    merge_sample;
    dedup;
    index;
    make_bed;
    merge_all;
    get_barcodes;
    soupercell;
} from './processes.nf'

workflow {

    if ( "${params.results}" == "false" ){error "Must provide parameter: results"}
    if ( "${params.samplesheet}" == "false" ){error "Must provide parameter: samplesheet"}
    if ( "${params.genome_fasta}" == "false" ){error "Must provide parameter: genome_fasta"}

    // Get the input BAM/SAM files from the samplesheet
    Channel
        .fromPath(
            "${params.samplesheet}",
            checkIfExists: true
        )
        .splitCsv(
            header: true,
            sep: params.samplesheet_sep,
            strip: true
        )
        .map {
            it -> [
                "${it[params.sample_col]}",
                file(
                    "${it[params.path_col]}",
                    checkIfExists: true
                )
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

    // Remove duplicates
    dedup(bam_ch)

    // Add unique tags for each input file
    add_tags(dedup.out)

    // Merge
    merge_sample(
        add_tags.out.groupTuple(
            sort: true
        )
    )

    // Index the BAM
    index(merge_sample.out)

    // Make a channel containing:
    //    tuple val(sample), path(bam), path(bai)
    merge_sample.out.join(index.out).set { indexed_bam }

    // Make a BED file from the BAM with its index
    make_bed(indexed_bam)

    // Merge the sample-level BAMs together
    merge_all(indexed_bam.toSortedList())

    // Get the barcodes which were actually used
    get_barcodes(merge_all.out)

    // Make sure that the genome FASTA exists
    genome = file(
        "${params.genome_fasta}",
        checkIfExists: true
    )

    // Make sure that the genome index exists
    genome_index = file(
        "${params.genome_fasta}.fai",
        checkIfExists: true
    )

    // Run soupercell
    soupercell(
        merge_all.out,
        get_barcodes.out,
        genome,
        genome_index
    )

}