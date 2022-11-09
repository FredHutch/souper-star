#!/bin/bash

set -e

samtools view "${1}" \
    | perl -nle '@reads=split(/\t/,$_); if (m/CB:Z:([^\t\n]+)\t/) { print $1; }' \
    | sort -u \
    | gzip -c
