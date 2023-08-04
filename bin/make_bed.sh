#!/bin/bash

set -e

# Perform sorting by read name and pipe into Perl for processing
samtools sort -n "$1" | \
samtools view -h | \
perl -nle '
    @reads = split(/\t/, $_);

    # If the line starts with '@', it's a header line, so just print it
    if ($_ =~ /^@/) {
        print;
        next;
    }

    # Check for the specific pattern in the line
    if (m/CB:Z:([^\t\n]+)\t/) {
        # Store the value found in the pattern (e.g., "GGTCATTTACTTGTTGATTACTCGCCTATCCT-18")
        $new_barcode = $1;

        # Replace the barcode part of the read name with the new barcode
        $reads[0] =~ s/_[A-Za-z0-9_]+$/_$new_barcode/;

        # Join the modified fields back together and print
        print join("\t", @reads);
    } else {
        # If pattern not found, print line unchanged
        print;
    }
' |  bedtools bamtobed -bedpe -i stdin | \
  cut -f1,2,6,7 | \
sort -k1,1 -k2n,2n -k3n,3n | \

  bgzip -c



