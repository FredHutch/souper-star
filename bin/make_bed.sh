#!/bin/bash

set -e

# Extract filename without extension to use in the intermediate output file path
filename=$(basename -- "$1")
filename_no_ext="${filename%.*}"

samtools sort -n "$1" | \
samtools view -h | \
awk -F'\t' '{
    if ($0 ~ /^@/) { print; next }
    if (match($0, /CB:Z:([^\t\n]+)\t/, arr)) {
        split($1, a, "_");
        modified_read_name = a[1] "_" arr[1];
        $1 = modified_read_name;
    }
    print $0
}' > "/home/djanssen/henikoff/for_dj/DJ_Hs_CD34_Exp123/souper-star_CSCB_D1-6_500reads_080323/intermediate_output_${filename_no_ext}.txt"











