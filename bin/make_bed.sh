#!/bin/bash

set -e

#!/bin/bash

samtools sort -n "$1" | \
  samtools view -h | \
  bedtools bamtobed -bedpe -i stdin | \
  cut -f1,2,6,7 | \
  bgzip -c


