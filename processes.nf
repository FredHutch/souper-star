process sam_to_bam {
    container "${params.container__samtools}"
    label "cpu_large"
    tag "${sam}"
    
    input:
        tuple val(sample), path(sam), val(ix)

    output:
        tuple val(sample), path("${sam.name.replaceAll(/.sam$/, '')}.bam"), val(ix)

    script:
    """
    samtools sort -@ ${task.cpus} "${sam}" -o "${sam.name.replaceAll(/.sam$/, '')}.bam"
    """
}

process add_tags {
    container "${params.container__misc}"
    label "io_limited"
    tag "${bam}"

    input:
        tuple val(sample), path(bam), val(ix)

    output:
        tuple val(sample), path("${bam.name.replaceAll(/.bam$/, '')}.tagged.bam")

    script:
    """
    add_tags.py -u ${params.umi_len} -i ${ix} "${bam}"\
    | samtools view -b - \
    > "${bam.name.replaceAll(/.bam$/, '')}.tagged.bam"
    """
}

process merge_sample {
    container "${params.container__samtools}"
    label "cpu_large"
    tag "${sample}"

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
    tag "${bam}"

    input:
        tuple val(sample), path(bam), val(ix)

    output:
        tuple val(sample), path("${bam.name.replaceAll(/.bam$/, '')}.dedup.bam"), val(ix)

    script:
    """
samtools sort -n -m 2G -@ ${task.cpus} "${bam}" \
    | samtools fixmate -m -@ ${task.cpus} - - \
    | samtools sort -m 2G -@ ${task.cpus} - \
    | samtools markdup -r -s -f dup.out --barcode-name -@ ${task.cpus} \
        - "${bam.name.replaceAll(/.bam$/, '')}.dedup.bam"
    """
}

process index {
    container "${params.container__samtools}"
    label "cpu_large"
    tag "${sample}"

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
    container "${params.container__misc}"
    label "io_limited"
    tag "${sample}"

    input:
        tuple val(sample), path(bam), path(bai)

    output:
        tuple val(sample), path("${sample}.bed.gz")

    script:
    """
make_bed.sh "${bam}" > "${sample}.bed.gz"
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
    container "${params.container__misc}"
    label "io_limited"
    tag "${bam}"

    input:
        tuple val(sample), path(bam)

    output:
        path "${bam.name.replaceAll(/.bam$/, '')}.barcodes.tsv.gz"

    script:
    """
get_barcodes.sh "${bam}" > "${bam.name.replaceAll(/.bam$/, '')}.barcodes.tsv.gz"
    """
}

process join_barcodes {
    container "${params.container__misc}"
    label "io_limited"

    input:
        path "input/"

    output:
        path "barcodes.tsv.gz"

    script:
    """
cat input/*.barcodes.tsv.gz > barcodes.tsv.gz
    """
}

process souporcell {
    publishDir "${params.results}", mode: 'copy', overwrite: true
    container "${params.container__souporcell}"
    label "cpu_large"

    input:
        path "*"
        path "barcodes.tsv.gz"
        path "genome.fa"
        path "genome.fa.fai"

    output:
        path "souporcell/*"

    script:
    """
#!/bin/bash
set -e

# Set up a local temporary directory
mkdir tmp
export TMPDIR=$PWD/tmp/

gunzip -c barcodes.tsv.gz > barcodes.tsv

souporcell_pipeline.py \
    -i merged.bam \
    -b barcodes.tsv \
    -f genome.fa \
    -t ${task.cpus} \
    -o souporcell \
    -k ${params.k} \
    --no_umi \
    ${params.flags}
    """
}