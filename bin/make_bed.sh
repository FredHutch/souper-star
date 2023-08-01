#!/bin/bash

set -e

samtools view "${1}" \
    | perl -nle '
        @reads = split(/\t/, $_); 
        if ((m/CB:Z:([^\t\n]+)\t/)) {
            $reads_length = length($reads[9]);
            $TLEN = $reads[8];
            # Check if the read is aligned
            if (($reads[1] & 4) == 0 && $TLEN >= 0) {
                $start = $reads[3] - 1;
                $end = $reads[3] + abs($TLEN) - 1;
                print "$reads[2]\t$start\t$end\t$1";
            } 
        }' \
    | bgzip -c
