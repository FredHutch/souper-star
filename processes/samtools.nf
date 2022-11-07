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