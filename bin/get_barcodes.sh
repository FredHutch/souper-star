#!/bin/bash

set -euo pipefail

mkdir tmp

export TMPDIR=$PWD/tmp/

samtools view merged.bam \
    | perl -nle '@reads=split(/\t/,$_); if (m/CB:Z:([^\t\n]+)\t/) { print $1; }' \
    | sort -u \
    | gzip -c \
    > barcodes.tsv.gz
