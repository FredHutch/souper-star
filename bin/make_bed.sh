#!/bin/bash

set -e

samtools view "${1}" \
    | perl -nle '@reads=split(/\t/,$_); if (m/CB:Z:([^\t\n]+)\t/) {$start=$reads[3]-1; $end=$reads[3]+length($reads[9])-1; print "$reads[2]\t","$start\t","$end\t",$1; }' \
    | bgzip -c
