#!/bin/bash

set -e

samtools sort -n "$1" | \
samtools view -h | \
awk -F'\t' '{
    if ($0 ~ /^@/) { print; next }
    if (match($0, /CB:Z:([^\t\n]+)\t/, arr)) {
        split($1, a, "_");
        $1 = a[1] "_" arr[1];
    }
    print
}' | \
samtools view -b | \
# Convert to BEDPE format
bedtools bamtobed -bedpe -i stdin | \
cut -f1,2,6,7 | \
sort -k1,1 -k2n,2n -k3n,3n | \
bgzip -c










