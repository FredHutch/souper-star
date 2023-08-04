#!/bin/bash

set -e

# Sort the SAM file by read name
samtools sort -n "$1" | \
# Convert to BAM and add modified read name with barcode from CB tag
samtools view -h | \
awk -F'\t' '{
    if ($0 ~ /^@/) { print; next }
    if (match($0, /CB:Z:([^\t\n]+)\t/, arr)) {
        split($1, a, "_");
        $1 = a[1] "_" arr[1];
    }
    print $0
}' | \
samtools view -b | \
# Convert to BEDPE format
bedtools bamtobed -bedpe -i stdin | \
cut -f1,2,6,7 | \
sort -k1,1 -k2n,2n -k3n,3n | \
bgzip -c






