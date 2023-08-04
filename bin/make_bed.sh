#!/bin/bash

# Check if the input file is provided
if [ -z "$1" ]; then
    echo "Please provide a BAM file as input."
    exit 1
fi

samtools sort -n "$1" | \
samtools view -h | \
perl -nle '
    @reads = split(/\t/, $_);

    # If the line starts with "@", it's a header line, so just print it
    if ($_ =~ /^@/) {
        print;
        next;
    }

    # Check for the specific pattern in the line
    if (m/CB:Z:([^\t\n]+)\t/) {
        # Store the value found in the pattern
        $new_barcode = $1;

        # Replace the barcode part of the read name with the new barcode
        $reads[0] =~ s/_[A-Za-z0-9_]+$/_$new_barcode/;

        # Join the modified fields back together and print
        print join("\t", @reads);
    } else {
        # If the pattern is not found, print the line unchanged
        print;
    }
'




