#!/bin/bash

set -e

#module load BEDTools/2.30.0-GCC-12.2.0 
#module load SAMtools/1.17-GCC-12.2.0

samtools view -h "$1"  | \
bedtools bamtobed -bedpe -i stdin | \
perl -nle 'BEGIN { $header=1 }
if ($header and /^@/) { print; next }
if ($header and !/^@/) { $header=0 }
@fields = split(/\t/);
if ($fields[6] =~ m/CB:Z:([^\t\n]+)\t/) {
  $barcode = $1;
  print join("\t", $fields[0], $fields[1], $fields[5], $barcode);
}' | \
sort -k1,1 -k2n,2n -k3n,3n | \
bgzip -c
