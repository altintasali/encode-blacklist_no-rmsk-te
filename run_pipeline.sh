#!/usr/bin/env bash

set -eo pipefail

# Create log directory if it doesn't exist yet
mkdir -p logs

# Redirect master output to file while printing live to terminal
exec > >(tee -i logs/pipeline_execution.log) 2>&1

echo "======================================================================="
echo " STARTING PIPELINE: ENCODE Blacklist Repetitive Element Filter Engine"
echo " Start Time: $(date)"
echo "======================================================================="

echo -e "\n[STAGE 1/2] Running genomic sequence filtering (Bash module)..."
if [ -f "scripts/filter_blacklist.sh" ]; then
    chmod +x scripts/filter_blacklist.sh
    ./scripts/filter_blacklist.sh
else
    echo "ERROR: scripts/filter_blacklist.sh not found!" >&2
    exit 1
fi

echo -e "\n[STAGE 2/2] Running downstream statistical evaluation (R module)..."
if [ -f "scripts/analyze_repeats.R" ]; then
    Rscript scripts/analyze_repeats.R
else
    echo "ERROR: scripts/analyze_repeats.R not found!" >&2
    exit 1
fi

echo -e "\n======================================================================="
echo " PIPELINE COMPLETE!"
echo " End Time: $(date)"
echo " Master Log Saved: logs/pipeline_execution.log"
echo "======================================================================="
