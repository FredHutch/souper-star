process sam_to_bam {
    container "${params.container__samtools}"
    label "cpu_large"
    tag "${sam}"
    
    input:
        tuple path(sam), val(sample)

    output:
        tuple path("${sam}.bam"), val(sample)

    script:
    """
    samtools sort -@ ${task.cpus} "${sam}" -o "${sam}.bam"
    """
}

process add_tags {
    container "${params.container__simplesam}"
    label "io_limited"
    tag "${bam}"

    input:
        tuple path(bam), val(sample)

    output:
        tuple val(sample), path("${bam.name.replaceAll(/.bam$/, '')}.tagged.bam")

    script:
    """
    add_tags.py -u ${params.umi_len} -i ${task.index} "${bam}" > "${bam.name.replaceAll(/.bam$/, '')}.tagged.bam"
    """
}

process merge {
    container "${params.container__samtools}"
    label "cpu_large"

    input:
        tuple val(sample), path("inputs/")

    output:
        tuple val(sample), path("${sample}.bam")

    script:
    """
    samtools merge "${sample}.bam" inputs/* -@ ${task.cpus}
    """
}

process dedup {
    container "${params.container__samtools}"
    label "cpu_large"

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path("${bam.name.replaceAll(/.bam$/, '')}.dedup.bam")

    script:
    """
samtools sort -n -m 2G -@ ${task.cpus} "${bam}" \
    | samtools fixmate -m -@ ${task.cpus} - - \
    | samtools sort -m 2G -@ ${task.cpus} - \
    | samtools markdup -r -s -f dup.out --barcode-tag CB -@ ${task.cpus} \
        - "${bam.name.replaceAll(/.bam$/, '')}.dedup.bam"
    """
}

process index {
    container "${params.container__samtools}"
    label "cpu_large"

    input:
        tuple val(sample), path(bam)

    output:
        tuple val(sample), path("${bam}.bai")

    script:
    """
samtools index "${bam}"
    """
}

process make_bed {
    publishDir "${params.results}/${sample}/", mode: 'copy', overwrite: true
    container "${params.container__samtools}"
    label "io_limited"

    input:
        tuple val(sample), path(bam), path(bai)

    output:
        path "${sample}.bed"

    script:
    """
make_bed.sh "${bam}" > "${sample}.bed"
    """
}