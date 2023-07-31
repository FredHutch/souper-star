#!/bin/bash

set -e

samtools view "${1}" \
    | perl -nle '@reads = split(/\t/, $_); if ((m/CB:Z:([^\t\n]+)\t/) && ($reads[1] & 64)) {$start = $reads[3] - 1; $reads_length = length($reads[9]); $TLEN = $reads[8]; if ($TLEN >= 0) {$end = $reads[3] + abs($TLEN) - 1;} else {$end = $reads[3] + $reads_length - 1;} print "$reads[2]\t$start\t$end\t$1";}' \
    | bgzip -c
