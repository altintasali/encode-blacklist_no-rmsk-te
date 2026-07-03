#!/usr/bin/env bash

set -eo pipefail

# Define folder targets relative to repo root
RAW_DIR="data/raw"
TEMP_DIR="data/process_temp"
FILTERED_DIR="results/cleaned_blacklist"

# Ensure all subdirectories exist
mkdir -p "$RAW_DIR" "$TEMP_DIR" "$FILTERED_DIR"

# Define specific file paths
ORIGINAL_BLACKLIST="${RAW_DIR}/mm10-blacklist.v2.bed"
RMSK_DATA="${RAW_DIR}/rmsk.txt.gz"

ALL_RMSK_BED="${TEMP_DIR}/mm10_all-rmsk.bed"
PURE_TES="${TEMP_DIR}/mm10_only-TEs.bed"
TRAPPED_TES="${TEMP_DIR}/mm10_TEs_trapped_in_blacklist.bed"
TRAPPED_ALL_RMSK="${TEMP_DIR}/mm10_all-rmsk_trapped_in_blacklist.bed"

OUTPUT_BLACKLIST_NO_TE="${FILTERED_DIR}/mm10-blacklist.v2_no-TEs.bed"
OUTPUT_BLACKLIST_NO_RMSK="${FILTERED_DIR}/mm10-blacklist.v2_no-all-rmsk.bed"

echo "-> Step 1: Checking genomic resources..."
if [ ! -f "$ORIGINAL_BLACKLIST" ]; then
    echo "Downloading ENCODE mm10 v2 blacklist..."
    wget -q --show-progress https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz -O ${ORIGINAL_BLACKLIST}.gz
    gunzip ${ORIGINAL_BLACKLIST}.gz
else
    echo "   Original blacklist found."
fi

if [ ! -f "$RMSK_DATA" ]; then
    echo "Downloading mm10 RepeatMasker track from UCSC..."
    wget -q --show-progress http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz -O "$RMSK_DATA"
else
    echo "   RepeatMasker dataset found."
fi

echo "-> Step 2: Formatting entire RepeatMasker dataset to BED..."
zcat "$RMSK_DATA" | awk 'BEGIN {OFS="\t"} {print $6, $7, $8, $11, "0", $10, $12, $13}' > "$ALL_RMSK_BED"

echo "-> Step 3: Extracting pure TEs..."
awk 'BEGIN {FS=OFS="\t"} $7 ~ /^(LINE|SINE|LTR|DNA|Retroposon|RC)\??$/' "$ALL_RMSK_BED" > "$PURE_TES"

echo "-> Step 4: Intersecting unfiltered RepeatMasker with blacklist..."
bedtools intersect -u -wa -a "$ALL_RMSK_BED" -b "$ORIGINAL_BLACKLIST" > "$TRAPPED_ALL_RMSK"

echo "-> Step 5: Finding TEs trapped in the blacklist..."
bedtools intersect -u -wa -a "$PURE_TES" -b "$ORIGINAL_BLACKLIST" > "$TRAPPED_TES"

echo "-> Step 6: Subtracting TEs only..."
bedtools subtract -a "$ORIGINAL_BLACKLIST" -b "$PURE_TES" > "$OUTPUT_BLACKLIST_NO_TE"

echo "-> Step 7: Subtracting ALL RepeatMasker elements..."
bedtools subtract -a "$ORIGINAL_BLACKLIST" -b "$ALL_RMSK_BED" > "$OUTPUT_BLACKLIST_NO_RMSK"

echo "-> Step 8: Verifying calculated bounds..."
LINES_ORIG=$(wc -l < "$ORIGINAL_BLACKLIST" | xargs)
LINES_NO_TE=$(wc -l < "$OUTPUT_BLACKLIST_NO_TE" | xargs)
LINES_NO_RMSK=$(wc -l < "$OUTPUT_BLACKLIST_NO_RMSK" | xargs)

BP_ORIG=$(awk '{sum+=$3-$2} END {print sum}' "$ORIGINAL_BLACKLIST")
BP_NO_TE=$(awk '{sum+=$3-$2} END {print sum}' "$OUTPUT_BLACKLIST_NO_TE")
BP_NO_RMSK=$(awk '{sum+=$3-$2} END {print sum}' "$OUTPUT_BLACKLIST_NO_RMSK")

echo "-----------------------------------------------------------------------"
echo "Metrics                     | Original      | No TEs        | No rmsk"
echo "-----------------------------------------------------------------------"
printf "Total Genomic Blocks (Lines) | %-13s | %-13s | %-13s\n" "$LINES_ORIG" "$LINES_NO_TE" "$LINES_NO_RMSK"
printf "Total Blacklisted Size (bp)  | %-13s | %-13s | %-13s\n" "$BP_ORIG" "$BP_NO_TE" "$BP_NO_RMSK"
echo "-----------------------------------------------------------------------"
