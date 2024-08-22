#!/bin/bash
CORES=${1}
BAM_IN="${2}"
set -e

samtools sort -n -m 2G -@ ${CORES} "${BAM_IN}" \
	| samtools view -B -h - \
	| bedtools bamtobed -bedpe -i stdin | cut -f1,2,6,7 | sort -k1,1 -k2n,2n -k3n,3n | awk -v OFS='\t' '{len = $3 - $2 ; print $0, len }' \
	| bgzip -c
