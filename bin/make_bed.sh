#!/bin/bash

set -e

# Debug file
debug_file="/home/djanssen/henikoff/for_dj/DJ_Hs_CD34_Exp123/souper-star_CSCB_D1-6_500reads_080323/debug_output.txt"

# Sort the SAM file by read name
samtools sort -n "$1" | \
# Convert to BAM and add modified read name with barcode from CB tag
samtools view -h | \
awk -F'\t' -v debug_file="$debug_file" '{
    if ($0 ~ /^@/) { print; next }
    if (match($0, /CB:Z:([^\t\n]+)\t/, arr)) {
        print "Original: " $1 >> debug_file;  # Debug: print original read name
        split($1, a, "_");
        print "First part: " a[1] >> debug_file;  # Debug: print first part of read name
        print "Barcode: " arr[1] >> debug_file;  # Debug: print barcode
        $1 = a[1] "_" arr[1];
        print "Modified: " $1 >> debug_file;  # Debug: print modified read name
    }
    print $0
}' | \
samtools view -b | \
# Convert to BEDPE format
bedtools bamtobed -bedpe -i stdin | \
cut -f1,2,6,7 | \
sort -k1,1 -k2n,2n -k3n,3n | \
bgzip -c







