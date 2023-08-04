#!/bin/bash

set -e

samtools sort -n "$1" | \
samtools view -h | \
awk -F'\t' '{
    if ($0 ~ /^@/) { print; next }
    match($0, /CB:Z:([^\t\n]+)\t/, arr)
    if (arr[1] != "") {
        sub(/_[A-Za-z0-9_]+$/, "_" arr[1], $1)
    }
    print $0
}' | \
bedtools bamtobed -bedpe -i stdin | \
cut -f1,2,6,7 | \
sort -k1,1 -k2n,2n -k3n,3n | \
bgzip -c





