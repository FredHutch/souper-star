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
        tuple path("${bam.name.replaceAll(/.bam$/, '')}.tagged.bam"), val(sample)

    script:
    """
    add_tags.py -u ${params.umi_len} -i ${task.index} "${bam}" > "${bam.name.replaceAll(/.bam$/, '')}.tagged.bam"
    """
}

process merge {
    container "${params.container__samtools}"
    label "cpu_large"

    input:
        tuple val(mark), path("inputs/")

    output:
        path "${mark}.bam"

    script:
    """
    samtools merge "${mark}.bam" inputs/* -@ ${task.cpus}
    """
}

process dedup {
    container "${params.container__samtools}"
    label "cpu_large"

    input:
        path bam

    output:
        path "${bam.name.replaceAll(/.bam$/, '')}.dedup.bam"

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
        path bam

    output:
        path "${bam}.bai"

    script:
    """
samtools index "${bam}"
    """
}