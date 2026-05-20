#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Initialize Conda 
eval "$(conda shell.bash hook)"

# Output directories
CLEAN_DIR="../clean"
OUTPUT_DIR="../output"
ABUND_DIR="../abundances"

mkdir -p "$CLEAN_DIR" "$OUTPUT_DIR" "$ABUND_DIR"

# Define sample array
mapfile -t SAMPLES < <(
    for f in ../reads/*_R1.fastq; do
        basename "$f" _R1.fastq
    done
)

# Warn if reads are not where they should be
if [ ${#SAMPLES[@]} -eq 0 ]; then
    echo "No samples found in ../reads/"
    exit 1
fi

# Print sample names
echo "Detected samples:"
printf ' - %s\n' "${SAMPLES[@]}"

# Check to make sure every sample has the reverse read
for sample in "${SAMPLES[@]}"; do
    if [[ ! -f ../reads/"${sample}"_R2.fastq ]]; then
        echo "Missing paired read for sample: $sample"
        exit 1
    fi
done

# QC and polyG tail removal
echo "Running fastp for QC and trimming..."

conda activate fastp 

for sample in "${SAMPLES[@]}"; do 
    echo "Processing QC for: $sample"
    fastp \
        -i ../reads/"${sample}"_R1.fastq \
        -I ../reads/"${sample}"_R2.fastq \
        -o "$CLEAN_DIR"/"${sample}"_c_R1.fastq.gz \
        -O "$CLEAN_DIR"/"${sample}"_c_R2.fastq.gz \
        -h "${sample}".html \
        -j "${sample}".json \
        -g -w 9
done

conda deactivate

# Kallisto quantification
conda activate kallisto

if [[ ! -f dual.index.idx ]]; then
    echo "Running Kallisto indexing"
    kallisto index -i dual.index.idx peri.micro.ffn
fi

echo "Running Kallisto quantification..."
for sample in "${SAMPLES[@]}"; do
    echo "Quantifying: $sample"
    kallisto quant \
        -i dual.index.idx \
        -o "$OUTPUT_DIR"/"$sample" \
        -t 14 -b 100 \
        "$CLEAN_DIR"/"${sample}"_c_R1.fastq.gz \
        "$CLEAN_DIR"/"${sample}"_c_R2.fastq.gz
done

conda deactivate

# Aggregating data
echo "Organizing abundance files..."

# Copy and rename abundance files per sample name
# Warn if run fails for any samples

for f in "$OUTPUT_DIR"/*; do
    base=$(basename "$f")

    if [[ -f "$f/abundance.tsv" ]]; then
        cp "$f/abundance.tsv" "$ABUND_DIR/$base.abund.tsv"
    else
        echo "Warning: missing abundance.tsv for $base"
    fi
done

# Group samples by condition
# target_id=f1, length = f2, eff_length = f3, est_counts=f4, tpm = f5
# Change to your sample conditions

echo "Aggregating output..."
cd "$ABUND_DIR"

CONDITIONS=("C_06h" "C_24h" "C_72h" "E_06h" "E_24h" "E_72h")

for cond in "${CONDITIONS[@]}"; do
    {
        echo -e "target_id\t${cond}_f1\t${cond}_f2\t${cond}_f3\t${cond}_f4"

        paste \
            <(tail -n +2 "${cond}_f1.abund.tsv" | cut -f1,5) \
            <(tail -n +2 "${cond}_f2.abund.tsv" | cut -f5) \
            <(tail -n +2 "${cond}_f3.abund.tsv" | cut -f5) \
            <(tail -n +2 "${cond}_f4.abund.tsv" | cut -f5)

    } > "${cond}.tpm.tsv"
done

echo "Enjoy your quantified data!"