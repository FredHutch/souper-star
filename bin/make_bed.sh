#!/bin/bash

# Modify the QNAME to include the barcode as a comment, then convert to BED
# Extract chromosome, start, end, and barcode and write to compressed output BED file
samtools view -h "$1" \
    | awk -v OFS='\t' '{
        if ($0 ~ /^@/) {
            print $0
            next
        }
        split($1, arr, "_")
        barcode = arr[2]
        # Add the barcode as a comment in the QNAME
        print $1 "#CB:" barcode, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15
      }' \
    | samtools view -b -h - \
    | bedtools bamtobed -bedpe -i stdin \
    | awk -v OFS='\t' '{
        # Extract barcode from QNAME
        split($4, arr, "#CB:")
        split(arr[2], barcode_arr, "_")
        barcode = barcode_arr[1]
        # Print chromosome, start, end, and barcode
        print $1, $2, $6, barcode
      }' \
    | sort -k1,1 -k2n,2n -k3n,3n \
    | bgzip -c
