#!/bin/bash

set -e

samtools view "${1}" \
    | perl -nle '@reads=split(/\t/,$_); if (m/CB:Z:([^\t\n]+)\t/) { print "$reads[2]\t","$reads[3]\t","$reads[7]\t",$1; }' \
    | gzip -c
