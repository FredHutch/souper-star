#!/bin/bash
set -e

samtools view -B -h -f 0x40 "${1}" \
    | perl -nle '@reads=split(/\t/,$_); if ($_ =~ /CB:Z:([^\t\n]+)/ && $reads[8] > 0) { $start=$reads[3]-1; $end=$start + $reads[8]; print "$reads[2]\t$start\t$end\t$1"; }' \
    | LC_COLLATE=C sort -k1,1 -k2n,2n -k3n,3n | bgzip -c
