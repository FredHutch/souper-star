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

process merge_sample {
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
        tuple val(sample), path("${sample}.bed")

    script:
    """
make_bed.sh "${bam}" > "${sample}.bed"
    """
}

process merge_all {
    container "${params.container__samtools}"
    label "cpu_large"

    input:
        path "inputs/"

    output:
        path "merged.*"

    script:
    """#!/bin/bash
set -e
samtools merge merged.bam inputs/*.bam -@ ${task.cpus}
samtools index merged.bam
    """
}

process get_barcodes {
    container "${params.container__samtools}"
    label "io_limited"

    input:
        path "*"

    output:
        path "barcodes.tsv"

    script:
    """
get_barcodes.sh merged.bam > barcodes.tsv
    """
}

process soupercell {
    publishDir "${params.results}", mode: 'copy', overwrite: true
    container "${params.container__soupercell}"
    label "cpu_large"

    input:
        path "*"
        path "barcodes.tsv"
        path "genome.fa"
        path "genome.fa.fai"

    output:
        path "*"

    script:
    """
souporcell_pipeline.py \
    -i merged.bam \
    -b barcodex.tsv \
    -f genome.fa \
    -t ${task.cpus} \
    -o ./ \
    -k ${params.k} \
    ${params.flags}
    """
}