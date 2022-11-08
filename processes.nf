process sam_to_bam {
    container "${params.container__samtools}"
    label "cpu_large"
    tag "${sam}"
    
    input:
        path sam

    output:
        path "${sam}.bam"

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
        path bam

    output:
        path "${bam.replaceAll(/.bam$/, '')}.tagged.bam"

    script:
    """
    add_tags.py -u ${params.umi_len} -i ${task.index} "${bam}" > "${bam.replaceAll(/.bam$/, '')}.tagged.bam"
    """
}