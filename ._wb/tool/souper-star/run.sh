#!/bin/bash

set -euo pipefail

date
echo
echo "Running workflow from ${PWD}"
echo

# Run the workflow
echo Starting workflow
nextflow \
    run \
    "${TOOL_REPO}" \
    --results "${PWD}" \
    -params-file ._wb/tool/params.json \
    -resume

echo
date
echo Done
