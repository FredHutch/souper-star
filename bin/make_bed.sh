#!/bin/bash

set -e

# Extract filename without extension to use in the intermediate output file path
filename=$(basename -- "$1")
filename_no_ext="${filename%.*}"

debug_file="/home/djanssen/henikoff/for_dj/DJ_Hs_CD34_Exp123/souper-star_CSCB_D1-6_500reads_080323/debug.txt"

echo "Debugging $filename_no_ext" >> "$debug_file"

samtools sort -n "$1" | \
samtools view -h | \
awk -F'\t' -v debug_file="$debug_file" '{
    if ($0 ~ /^@/) { print; next }
    if (match($0, /CB:Z:([^\t\n]+)\t/, arr)) {
        print "Match found: " $0 >> debug_file;
        split($1, a, "_");
        $1 = a[1] "_" arr[1];
        # Rebuild $0
        $0 = $1;
        for(i=2; i<=NF; i++) {
            $0 = $0 FS $i;
        }
        print "Modified line: " $0 >> debug_file;
    } else {
        print "No match: " $0 >> debug_file;
    }
    print $0
}' > "/home/djanssen/henikoff/for_dj/DJ_Hs_CD34_Exp123/souper-star_CSCB_D1-6_500reads_080323/intermediate_output_${filename_no_ext}.txt"









