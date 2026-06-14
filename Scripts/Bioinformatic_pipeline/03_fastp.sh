#!/bin/bash

## Define variables for paths
DATA="/Volumes/MARTA_11/PP_Jorgi/fastp_controls"
TRIM_DIR="/Volumes/MARTA_11/PP_Jorgi/results/fastp2"
QC_TRIM_DIR="/Volumes/MARTA_11/PP_Jorgi/results/fastp/trim_qc_reports"
FAILED_DIR="/Volumes/MARTA_11/PP_Jorgi/results/fastp/failed_reads"

## Create output directories if they don't exist
mkdir -p "$TRIM_DIR" "$QC_TRIM_DIR" "$FAILED_DIR"

echo "Starting fastp (trimming + QC only, dedup will run post-host removal)..."

for R1 in "$DATA"/*_R1.fastq.gz; do

    base=$(basename "$R1" _R1.fastq.gz)
    R2="${DATA}/${base}_R2.fastq.gz"

    # Check R2 exists before proceeding
    if [[ ! -f "$R2" ]]; then
        echo "WARNING: R2 not found for $base, skipping..."
        continue
    fi

    echo "Processing sample: $base"

    fastp \
      -i "$R1" \
      -I "$R2" \
      -o "$TRIM_DIR/${base}_fastp_R1.fastq.gz" \
      -O "$TRIM_DIR/${base}_fastp_R2.fastq.gz" \
      --detect_adapter_for_pe \
      --correction \
      --length_required 50 \
      --low_complexity_filter \
      --complexity_threshold 30 \
      --html "$QC_TRIM_DIR/${base}_fastp_report.html" \
      --json "$QC_TRIM_DIR/${base}_fastp_report.json" \
      --failed_out "$FAILED_DIR/${base}_failed.fastq.gz" \
      --thread 7

    if [[ $? -ne 0 ]]; then
        echo "ERROR: fastp failed for sample $base"
        exit 1
    fi

    # Print read counts for tracking
    echo "Finished: $base"
    echo "  Input:  $(zcat "$R1" | awk 'NR%4==1' | wc -l) read pairs"
    echo "  Output: $(zcat "$TRIM_DIR/${base}_fastp_R1.fastq.gz" | awk 'NR%4==1' | wc -l) read pairs remaining"
    echo "---"

done

echo ""
echo "All samples processed."
echo "Run MultiQC with:"
echo "  multiqc $QC_TRIM_DIR -o $QC_TRIM_DIR/multiqc_report"
