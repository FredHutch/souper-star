#!/bin/bash

set -e

#!/bin/bash

samtools view -h "$1"  | bedtools bamtobed -bedpe -i stdin | cut -f1,2,6,7 | sort -k1,1 -k2n,2n -k3n,3n | bgzip -c
