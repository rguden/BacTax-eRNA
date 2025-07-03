#!/bin/bash
set -e  # Stop on error

echo "Setting up Conda environment"

# Parameters
THREADS=8
EMAIL="rodgeemae.guden@ugent.be"
START_TIME=$(date +"%Y-%m-%d %H:%M:%S")

# Input and output directories
INPUT_DIR="/home/rguden/Documents/ORCHESTRA/sequences/RNA_raw_data"
OUTPUT_DIR="/home/rguden/Documents/ORCHESTRA/RNA_output"
KRAKEN_DB="/home/rguden/Documents/databases/kraken2_standard"
LOG_DIR="$OUTPUT_DIR/logs"
LOG_FILE="$LOG_DIR/rRNA_pipeline.log"

# SortMeRNA database directory
SORTMERNA_DB_DIR="/home/rguden/Documents/databases/sortmerna"

# Ensure directories exist
mkdir -p "$OUTPUT_DIR/quality_control"
mkdir -p "$OUTPUT_DIR/trimming/fastp"
mkdir -p "$OUTPUT_DIR/trimming/fastqc/fastp"
mkdir -p "$OUTPUT_DIR/sortmerna"
mkdir -p "$LOG_DIR"
mkdir -p "$SORTMERNA_DB_DIR"


# Function to send email notifications
send_email() {
    echo "[INFO] Email disabled"
}

# Function to log messages
echo_log() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] $1" | tee -a "$LOG_FILE"
}

# Error Handling
trap 'send_email "FAILED" "Metatranscriptomics Pipeline failed at: $(date)"; exit 1' ERR
set -e  # Exit on error

# Start logging
echo_log "Pipeline started at $START_TIME"

# Parse command-line arguments (allow multiple steps)
while getopts "s:" opt; do
    case ${opt} in
        s)
            IFS=',' read -r -a STEPS < <(echo "$OPTARG")  # Split steps by commas
            ;;
        *)
            echo_log "Usage: $0 -s <step1,step2,...>" | tee -a $LOG_FILE
            exit 1
            ;;
    esac
done


# Loop through each specified step
for STEP in "${STEPS[@]}"; do
    case "$STEP" in

    qc)
        echo_log "Step 2: Quality Control (FastQC)" | tee -a $LOG_FILE
        fastqc --noextract -t 8 -o $OUTPUT_DIR/quality_control $INPUT_DIR/*.fastq.gz &>> $LOG_FILE
        multiqc $OUTPUT_DIR/quality_control -o $OUTPUT_DIR/quality_control/multiqc_report &>> $LOG_FILE
        echo_log "Step 2: Quality Control (FastQC) completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step 2: Quality Control (FastQC) completed successfully!"
        ;;

    trim_fastp)
        echo_log "Step 3: Trimming adpaters with fastp..." | tee -a $LOG_FILE
        for R1 in "$INPUT_DIR"/*_R1.fastq.gz; do
            if [[ -f "$R1" ]]; then
                R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"
                SAMPLE_NAME=$(basename "$R1" _R1.fastq.gz)

                FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastp"

                TRIMMED_R1="$FASTP_OUTPUT/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
                TRIMMED_R2="$FASTP_OUTPUT/${SAMPLE_NAME}_trimmed_R2.fastq.gz"
                JSON_REPORT="$FASTP_OUTPUT/${SAMPLE_NAME}_fastp.json"
                HTML_REPORT="$FASTP_OUTPUT/${SAMPLE_NAME}_fastp.html"

                echo_log "Running fastp for $SAMPLE_NAME..." | tee -a $LOG_FILE
                fastp --in1 "$R1" --in2 "$R2" \
                      --out1 "$TRIMMED_R1" --out2 "$TRIMMED_R2" \
                      --detect_adapter_for_pe \
                      --length_required 30 \
                      --thread $THREADS \
                      --json "$JSON_REPORT" \
                      --html "$HTML_REPORT" &>> "$LOG_FILE"
                echo_log "fastp trimming completed for $SAMPLE_NAME!" | tee -a $LOG_FILE
            fi
        done
        echo_log "Step 3: Trimming with fastp completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step 3: Trimming with fastp completed successfully!"
        ;;

    postqc_trim_fastp)
        echo_log "Step 4: Running FastQC and MultiQC after fastp trimming..." | tee -a $LOG_FILE

        FASTQC_FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastqc/fastp"

        fastqc --noextract -t $THREADS -o "$FASTQC_FASTP_OUTPUT" "$FASTP_OUTPUT"/*.fastq.gz &>> "$LOG_FILE"
        multiqc "$FASTQC_FASTP_OUTPUT" -o "$FASTQC_FASTP_OUTPUT/multiqc_report" &>> "$LOG_FILE"

        echo_log "Step 4: Post-fastp Quality Control completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step 4: Post-fastp Quality Control completed successfully!"
        ;;

    download_db)
        echo_log "Step 7: Downloading SortMeRNA database..." | tee -a $LOG_FILE
        wget -O "$SORTMERNA_DB_DIR/database.tar.gz" "https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz"
        tar -xzf "$SORTMERNA_DB_DIR/database.tar.gz" -C "$SORTMERNA_DB_DIR"
        rm "$SORTMERNA_DB_DIR/database.tar.gz"
        echo_log "SortMeRNA database downloaded!" | tee -a $LOG_FILE
        ;;

    index_db)
        echo_log "Step 8: Indexing SortMeRNA database..." | tee -a $LOG_FILE
        indexdb_rna --ref "$SORTMERNA_DB_DIR/smr_v4.3_default_db.fasta","$SORTMERNA_DB_DIR/smr_v4.3_default_db.idx"
        echo_log "Indexing completed!" | tee -a $LOG_FILE
        ;;

    remove_rrna_unmerged)
        echo_log "Step 10: Removing rRNA for unmerged fast-p reads..." | tee -a $LOG_FILE


        FASTP_OUTPUT="$OUTPUT_DIR/trimming/fastp"
        SORTMERNA_OUTPUT="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp"

        IDX_DIR="$SORTMERNA_OUTPUT/idx"
        KVDB_DIR="$SORTMERNA_OUTPUT/kvdb"
        READB_DIR="$SORTMERNA_OUTPUT/readb"
        OUT_ALIGNED="$SORTMERNA_OUTPUT/out/aligned"
        OUT_MRNA="$SORTMERNA_OUTPUT/out/mrna"

        mkdir -p "$SORTMERNA_OUTPUT"
        mkdir -p "$IDX_DIR"
        mkdir -p "$KVDB_DIR"
        mkdir -p "$READB_DIR"
        mkdir -p "$OUT_ALIGNED"
        mkdir -p "$OUT_MRNA"

        # Loop through paired-end files
        for forward in "$FASTP_OUTPUT/"*_R1.fastq.gz; do
            if [[ -f "$forward" ]]; then
                SAMPLE_NAME=$(basename "$forward" _R1.fastq.gz)
                reverse="$FASTP_OUTPUT/${SAMPLE_NAME}_R2.fastq.gz"

                if [[ ! -f "$reverse" ]]; then
                    echo_log "Skipping $SAMPLE_NAME: Missing reverse read ($reverse)!" | tee -a "$LOG_FILE"
                    continue
                fi

                # Create per-sample directories
                SAMPLE_KVDB="$KVDB_DIR/$SAMPLE_NAME"
                SAMPLE_READB="$READB_DIR/$SAMPLE_NAME"
                mkdir -p "$SAMPLE_KVDB" "$SAMPLE_READB"


        echo_log "Processing $SAMPLE_NAME for rRNA removal..." | tee -a "$LOG_FILE"

        # Run SortMeRNA separately on forward and reverse reads
        sortmerna \
            --ref "$SORTMERNA_DB_DIR/silva-euk-28s-id98.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-euk-18s-id95.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-bac-23s-id98.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-bac-16s-id90.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-arc-23s-id98.fasta" \
            --ref "$SORTMERNA_DB_DIR/silva-arc-16s-id95.fasta" \
            --ref "$SORTMERNA_DB_DIR/rfam-5s-database-id98.fasta" \
            --ref "$SORTMERNA_DB_DIR/rfam-5.8s-database-id98.fasta" \
            --reads "$forward" --reads "$reverse" \
            --out2 \
            --paired_out \
            --fastx \
            --blast "1 cigar qcov" \
            --aligned "$OUT_ALIGNED/${SAMPLE_NAME}_aligned" \
            --other "$OUT_MRNA/${SAMPLE_NAME}_mrna" \
            --kvdb "$SAMPLE_KVDB" \
            --readb "$SAMPLE_READB" \
            --v \
            --threads "$THREADS" --index 2 \
            --workdir "$SORTMERNA_OUTPUT" \
            --idx-dir "$IDX_DIR" &>> "$LOG_FILE"

        # Capture exit status
        SORTMERNA_EXIT=$?
        if [[ $SORTMERNA_EXIT -ne 0 ]]; then
            echo_log "Error: SortMeRNA failed for $SAMPLE_NAME with exit code $SORTMERNA_EXIT" | tee -a "$LOG_FILE"
            send_email "FAILED" "Step 10: SortMeRNA failed for $SAMPLE_NAME (Exit Code: $SORTMERNA_EXIT)"
            exit 1
        fi

            echo_log "rRNA removal completed for $SAMPLE_NAME!" | tee -a $LOG_FILE
            else
            echo_log "Skipping $SAMPLE_NAME: File does not exist!" | tee -a $LOG_FILE
            fi
        done
        ;;

    postqc_sortmerna)
        echo_log "Step: Running FastQC and MultiQC after removing rRNA..." | tee -a $LOG_FILE

        SORTMERNA_OUTPUT="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp"
        OUT_MRNA="$SORTMERNA_OUTPUT/out/mrna"
        SORTMERNA_QC="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp/qc"

        fastqc --noextract -t $THREADS -o "$SORTMERNA_QC" "$OUT_MRNA"/*.fq.gz &>> "$LOG_FILE"
        multiqc "$SORTMERNA_QC" -o "$SORTMERNA_QC/multiqc_report" &>> "$LOG_FILE"

        echo_log "Step: Post-sortmerna Quality Control completed successfully!" | tee -a $LOG_FILE
        send_email "COMPLETED" "Step: Post-sorting Quality Control completed successfully!"
        ;;

    summary_sortmerna)
        echo_log "Step: Summarizing sortmerna results..." | tee -a $LOG_FILE

        # Directories
        FASTP_OUTPUT="/home/rguden/Documents/ORCHESTRA/RNA_output/trimming/fastp"
        SORTMERNA_OUTPUT="/home/rguden/Documents/ORCHESTRA/RNA_output/sortmerna/sortmerna_unmerged_fastp"
        OUT_MRNA="$SORTMERNA_OUTPUT/out/mrna"
        OUT_RRNA="$SORTMERNA_OUTPUT/out/aligned"

        # Sample metadata file (make sure this is tab-separated!)
        METADATA_FILE="/home/rguden/Documents/ORCHESTRA/metadata/RNA/sample_condition.tsv"

        # Output summary files
        SUMMARY_FILE="$SORTMERNA_OUTPUT/sortmerna_summary.tsv"
        SUMMARY_COMBINED="$SORTMERNA_OUTPUT/sortmerna_summary_combined.tsv"

        # Create TSV headers (now include condition)
        echo -e "Sample\tCondition\tTotal_Reads_Trimmed_R1\tmRNA_Reads_R1\trRNA_Reads_R1\tmRNA_Percentage_R1\trRNA_Percentage_R1\tTotal_Reads_Trimmed_R2\tmRNA_Reads_R2\trRNA_Reads_R2\tmRNA_Percentage_R2\trRNA_Percentage_R2" > "$SUMMARY_FILE"
        echo -e "Sample\tCondition\tTotal_Reads_Trimmed\tmRNA_Reads\trRNA_Reads\tmRNA_Percentage\trRNA_Percentage" > "$SUMMARY_COMBINED"

        # Read metadata into associative arrays (requires Bash >= 4)
        declare -A SAMPLE_TO_CONDITION

        while IFS=$'\t' read -r sample condition replicate; do
            [[ $sample == "sample" ]] && continue  # skip header
            SAMPLE_TO_CONDITION["$sample"]="$condition"
        done < "$METADATA_FILE"

        # Loop through forward reads (_R1)
        for forward in "$FASTP_OUTPUT/"*_R1.fastq.gz; do
            if [[ -f "$forward" ]]; then
                # Get the true file base (with _trimmed) for file lookup
                SAMPLE_FILEBASE=$(basename "$forward" _R1.fastq.gz)
                # Get the sample name (strip _trimmed) for metadata lookup
                SAMPLE_NAME=${SAMPLE_FILEBASE%_trimmed}

                # Now find the correct reverse filename (needs _trimmed)
                reverse="$FASTP_OUTPUT/${SAMPLE_FILEBASE}_R2.fastq.gz"

                # Metadata: use SAMPLE_NAME
                CONDITION="${SAMPLE_TO_CONDITION[$SAMPLE_NAME]}"
                [[ -z "$CONDITION" ]] && CONDITION="UNKNOWN"

                # Ensure reverse file exists
                if [[ ! -f "$reverse" ]]; then
                    echo "Skipping $SAMPLE_NAME: Missing reverse read!" | tee -a "$LOG_FILE"
                    continue
                fi

                # Read counts
                TOTAL_R1=$(zcat "$forward" | awk '{s++} END {print s/4}')
                TOTAL_R2=$(zcat "$reverse" | awk '{s++} END {print s/4}')

                # mRNA
                if [[ -f "$OUT_MRNA/${SAMPLE_FILEBASE}_mrna_fwd.fq.gz" ]]; then
                    MRNA_R1=$(zcat "$OUT_MRNA/${SAMPLE_FILEBASE}_mrna_fwd.fq.gz" | awk '{s++} END {print s/4}')
                else
                    MRNA_R1=0
                fi
                if [[ -f "$OUT_MRNA/${SAMPLE_FILEBASE}_mrna_rev.fq.gz" ]]; then
                    MRNA_R2=$(zcat "$OUT_MRNA/${SAMPLE_FILEBASE}_mrna_rev.fq.gz" | awk '{s++} END {print s/4}')
                else
                    MRNA_R2=0
                fi

                # rRNA
                if [[ -f "$OUT_RRNA/${SAMPLE_FILEBASE}_aligned_fwd.fq.gz" ]]; then
                    RRNA_R1=$(zcat "$OUT_RRNA/${SAMPLE_FILEBASE}_aligned_fwd.fq.gz" | awk '{s++} END {print s/4}')
                else
                    RRNA_R1=0
                fi
                if [[ -f "$OUT_RRNA/${SAMPLE_FILEBASE}_aligned_rev.fq.gz" ]]; then
                    RRNA_R2=$(zcat "$OUT_RRNA/${SAMPLE_FILEBASE}_aligned_rev.fq.gz" | awk '{s++} END {print s/4}')
                else
                    RRNA_R2=0
                fi

                # Percentages
                if [[ $TOTAL_R1 -gt 0 ]]; then
                    MRNA_PERC_R1=$(awk "BEGIN {print ($MRNA_R1/$TOTAL_R1)*100}")
                    RRNA_PERC_R1=$(awk "BEGIN {print ($RRNA_R1/$TOTAL_R1)*100}")
                else
                    MRNA_PERC_R1=0
                    RRNA_PERC_R1=0
                fi
                if [[ $TOTAL_R2 -gt 0 ]]; then
                    MRNA_PERC_R2=$(awk "BEGIN {print ($MRNA_R2/$TOTAL_R2)*100}")
                    RRNA_PERC_R2=$(awk "BEGIN {print ($RRNA_R2/$TOTAL_R2)*100}")
                else
                    MRNA_PERC_R2=0
                    RRNA_PERC_R2=0
                fi

                # Output detailed per-direction summary
                echo -e "$SAMPLE_NAME\t$CONDITION\t$TOTAL_R1\t$MRNA_R1\t$RRNA_R1\t$MRNA_PERC_R1\t$RRNA_PERC_R1\t$TOTAL_R2\t$MRNA_R2\t$RRNA_R2\t$MRNA_PERC_R2\t$RRNA_PERC_R2" >> "$SUMMARY_FILE"

                # Combined counts
                TOTAL_COMBINED=$((TOTAL_R1 + TOTAL_R2))
                MRNA_COMBINED=$((MRNA_R1 + MRNA_R2))
                RRNA_COMBINED=$((RRNA_R1 + RRNA_R2))

                if [[ $TOTAL_COMBINED -gt 0 ]]; then
                    MRNA_PERC_COMBINED=$(awk "BEGIN {print ($MRNA_COMBINED/$TOTAL_COMBINED)*100}")
                    RRNA_PERC_COMBINED=$(awk "BEGIN {print ($RRNA_COMBINED/$TOTAL_COMBINED)*100}")
                else
                    MRNA_PERC_COMBINED=0
                    RRNA_PERC_COMBINED=0
                fi

                echo -e "$SAMPLE_NAME\t$CONDITION\t$TOTAL_COMBINED\t$MRNA_COMBINED\t$RRNA_COMBINED\t$MRNA_PERC_COMBINED\t$RRNA_PERC_COMBINED" >> "$SUMMARY_COMBINED"
            fi
        done

        echo "SortMeRNA summary table generated: $SUMMARY_FILE"
        echo_log "Summarizing sortmerna results completed successfully" | tee -a $LOG_FILE
        ;;

    classify_rrna)
        echo_log "Taxonomic Classification with Kraken2 (rRNA reads)"
        KRAKEN_OUTPUT="$OUTPUT_DIR/community_analyses/kraken_rrna"
        rRNA_DIR="$OUTPUT_DIR/sortmerna/sortmerna_unmerged_fastp/out/aligned"
        mkdir -p "$KRAKEN_OUTPUT"

        for R1 in "$rRNA_DIR"/*_aligned_fwd.fq.gz; do
            if [[ -f "$R1" ]]; then
                SAMPLE_NAME=$(basename "$R1" _aligned_fwd.fq.gz)
                R2="$rRNA_DIR/${SAMPLE_NAME}_aligned_rev.fq.gz"
                if [[ ! -f "$R2" ]]; then
                    echo_log "Skipping $SAMPLE_NAME: Missing aligned reverse read" | tee -a "$LOG_FILE"
                    continue
                fi

                echo_log "Running Kraken2 on $SAMPLE_NAME (paired-end rRNA)..."

                kraken2 --db "$KRAKEN_DB" \
                    --threads "$THREADS" \
                    --paired \
                    --report "$KRAKEN_OUTPUT/${SAMPLE_NAME}_report.txt" \
                    --output "$KRAKEN_OUTPUT/${SAMPLE_NAME}_kraken.txt" \
                    "$R1" "$R2" &>> "$LOG_FILE"

                echo_log "Kraken2 classification completed for $SAMPLE_NAME!"
            fi
        done

        echo_log "Kraken2 classification completed for all rRNA samples!" | tee -a "$LOG_FILE"
        send_email "Kraken2 rRNA Step Completed" "Kraken2 classification for rRNA reads has completed at $(date '+%Y-%m-%d %H:%M:%S')."
        ;;


    krona_kraken)
        echo_log "Generating Krona visualization"
        mkdir -p "$OUTPUT_DIR/community_analyses/kraken_rrna/krona"

        # Load krona taxonomy database
        # ktUpdateTaxonomy.sh

            # Check if there are any formatted Krona input files
        if ls "$OUTPUT_DIR/community_analyses/kraken_rrna"/*_kraken.txt 1> /dev/null 2>&1; then
            # Generate a single Krona chart combining multiple samples
            echo_log "Generating combined Krona chart from multiple samples"

            # Run ktImportTaxonomy with multiple input files
            ktImportTaxonomy -t 5 -m 3 -o "$OUTPUT_DIR/community_analyses/kraken_rrna/krona/multi-krona.html" "$OUTPUT_DIR/community_analyses/kraken_rrna"/*_kraken.txt &>> "$LOG_FILE"

            # Check if the multi-sample Krona file was successfully created
            if [[ -f "$OUTPUT_DIR/community_analyses/kraken_rrna/krona/multi-krona.html" ]]; then
                echo_log "Generated combined Krona chart: $OUTPUT_DIR/community_analyses/kraken_rrna/krona/multi-krona.html"
            else
                echo_log "Error: Failed to generate combined Krona chart"
            fi
        else
            echo_log "Warning: No Krona input files found in $OUTPUT_DIR/community_analyses/kraken_rrna/krona"
        fi
        ;;

    filter_bacteria_rrna)
        echo_log "Step: Filtering Kraken2 rRNA reports to retain only Bacteria"

        KRAKEN_DIR="$OUTPUT_DIR/community_analyses/kraken_rrna"
        FILTERED_DIR="$KRAKEN_DIR/filtered_bacteria"
        mkdir -p "$FILTERED_DIR"

        for report in "$KRAKEN_DIR"/*_report.txt; do
            [[ -f "$report" ]] || continue

            SAMPLE_NAME=$(basename "$report" _report.txt)
            OUTFILE="$FILTERED_DIR/${SAMPLE_NAME}_report_bacteria_only.txt"

            # Find line number of 'Bacteria'
            BACT_LINE=$(awk -F'\t' '$4 ~ /^R2$/ && $6 ~ /Bacteria/ {print NR; exit}' "$report")
            if [[ -z "$BACT_LINE" ]]; then
                echo_log "Skipping $SAMPLE_NAME: No Bacteria entry found"
                continue
            fi

            # Get indentation of Bacteria line
            BACT_INDENT=$(awk -v line="$BACT_LINE" 'NR == line {
                match($6, /[^ ]/, m); print RSTART - 1
            }' "$report")

            # Output Bacteria subtree (greater indent and after Bacteria line)
            awk -v start="$BACT_LINE" -v indent="$BACT_INDENT" -F'\t' '
                NR == start { print; next }
                NR > start {
                    match($6, /[^ ]/, m);
                    cur_indent = RSTART - 1;
                    if (cur_indent > indent) print;
                    else exit
                }' "$report" > "$OUTFILE"

            echo_log "Filtered Bacteria-only report written: $OUTFILE"
        done

        echo_log "Filtering complete for all Kraken2 reports." | tee -a "$LOG_FILE"
        ;;

    generate_bacteria_otu_R)
        echo_log "Step: Generating species- and genus-level OTU tables from Kraken2 Bacteria reports"

        R_SCRIPT="/home/rguden/Documents/ORCHESTRA/RNA_output/community_analyses/kraken_rrna/filtered_bacteria/combine_bacteria_reports.R"
        FILTERED_DIR="$OUTPUT_DIR/community_analyses/kraken_rrna/filtered_bacteria"
        OTU_DIR="$OUTPUT_DIR/community_analyses/kraken_rrna/filtered_bacteria/otu_table"
        mkdir -p "$OTU_DIR"

        OTU_SPECIES="$OTU_DIR/otu_table_species.tsv"
        OTU_GENUS="$OTU_DIR/otu_table_genus.tsv"

        THREADS=${THREADS:-6}  # Default to 4 threads if not set
        
        export TMPDIR="$HOME/tmp_orchestra"
        mkdir -p "$TMPDIR"

        if [[ ! -f "$R_SCRIPT" ]]; then
            echo_log "ERROR: R script not found at $R_SCRIPT"
            exit 1
        fi

        echo_log "Running R script to create species and genus OTU tables..."
        Rscript "$R_SCRIPT" "$FILTERED_DIR" "$OTU_SPECIES" "$OTU_GENUS" "$THREADS" \
            && echo_log "✅ OTU tables written to $OTU_DIR" \
            || echo_log "❌ R script failed for OTU generation"

        send_email "Genus/Species OTU Table Complete" "OTU tables (genus, species) have been generated at $OTU_DIR"
        ;;

    autofix_otu_header)
        echo_log "Step: Auto-fixing OTU table sample names for genus and species" | tee -a "$LOG_FILE"

        OTU_DIR="$OUTPUT_DIR/community_analyses/kraken_rrna/filtered_bacteria/otu_table"
        METADATA="/home/rguden/Documents/ORCHESTRA/metadata/RNA/sample_condition.tsv"

        # Get metadata sample IDs (excluding header)
        meta_samples=($(awk 'NR>1{print $1}' "$METADATA"))

        for LEVEL in genus species; do
            OTU_TABLE="$OTU_DIR/otu_table_${LEVEL}.tsv"
            OUT_OTU="${OTU_TABLE%.tsv}_matchmeta.tsv"

            if [[ ! -f "$OTU_TABLE" ]]; then
                echo_log "Skipping $LEVEL: File not found: $OTU_TABLE"
                continue
            fi

            echo_log "Auto-fixing sample names in $OTU_TABLE"

            # Read header and clean sample names
            IFS=$'\t' read -r -a header < <(head -n1 "$OTU_TABLE")
            new_header=("${header[0]}")

            for ((i=1; i<${#header[@]}; i++)); do
                clean=${header[i]}
                clean=${clean%_trimmed}
                clean=${clean%_R1}
                clean=${clean%_R2}
                clean=${clean%_1}
                clean=${clean%_2}
                # You can add more suffix patterns here if needed

                new_header+=("$clean")
            done

            # Write new OTU table with fixed header
            {
                printf "%s" "${new_header[0]}"
                for ((i=1; i<${#new_header[@]}; i++)); do
                    printf "\t%s" "${new_header[i]}"
                done
                printf "\n"
                tail -n +2 "$OTU_TABLE"
            } > "$OUT_OTU"

            echo_log "✅ Cleaned OTU table written: $OUT_OTU"
        done

        send_email "OTU Table Header Fix Complete" "Genus and species OTU tables have been matched to metadata format."
        ;;

    run_diversity_analysis)
        echo_log "Step: Diversity Analysis (Alpha & Beta)" | tee -a "$LOG_FILE"

        OTU_TABLE="/home/rguden/Documents/ORCHESTRA/RNA_output/community_analyses/kraken_rrna/filtered_bacteria/otu_table/otu_table_genus_matchmeta.tsv"
        METADATA="/home/rguden/Documents/ORCHESTRA/metadata/RNA/sample_condition.tsv"
        R_SCRIPT="/home/rguden/Documents/ORCHESTRA/RNA_output/community_analyses/kraken_rrna/filtered_bacteria/otu_diversity_analyses.R"
        DIVERSITY_OUTDIR="/home/rguden/Documents/ORCHESTRA/RNA_output/community_analyses/kraken_rrna/filtered_bacteria/genus_diversity_results_double_filtered"
        R_LOG="$DIVERSITY_OUTDIR/species_diversity_analysis_R_error.log"

        mkdir -p "$DIVERSITY_OUTDIR"
        
        export TMPDIR="$HOME/tmp_orchestra"
        mkdir -p "$TMPDIR"

        if [[ ! -f "$R_SCRIPT" ]]; then
            echo_log "ERROR: R script not found at $R_SCRIPT" | tee -a "$LOG_FILE"
            exit 1
        fi

        echo_log "Running diversity analysis with Rscript..."
        if ! Rscript "$R_SCRIPT" "$OTU_TABLE" "$METADATA" "$DIVERSITY_OUTDIR" &> "$R_LOG"; then
            echo_log "❌ ERROR: Diversity analysis failed. Check R log: $R_LOG" | tee -a "$LOG_FILE"
            echo_log "Here is the last few lines of the R error log:" | tee -a "$LOG_FILE"
            tail -n 10 "$R_LOG" | tee -a "$LOG_FILE"
            exit 1
        fi

        echo_log "✅ Diversity analysis completed successfully. Results in: $DIVERSITY_OUTDIR"
        ;;

    run_deseq2)
        echo_log "Step: Differential Abundance Analysis with DESeq2" | tee -a "$LOG_FILE"

        OTU_TABLE="/home/rguden/Documents/ORCHESTRA/RNA_output/community_analyses/kraken_rrna/filtered_bacteria/otu_table/otu_table_genus_matchmeta.tsv"
        METADATA="/home/rguden/Documents/ORCHESTRA/metadata/RNA/sample_condition.tsv"
        R_SCRIPT="/home/rguden/Documents/ORCHESTRA/RNA_output/community_analyses/kraken_rrna/filtered_bacteria/rrna_deseq2.R"
        DESEQ2_OUTDIR="/home/rguden/Documents/ORCHESTRA/RNA_output/community_analyses/kraken_rrna/filtered_bacteria/genus_deseq2_results_double_filtered"
        R_LOG="$DESEQ2_OUTDIR/deseq2_differential_abundance_R_error.log"

        mkdir -p "$DESEQ2_OUTDIR"
        
        export TMPDIR="$HOME/tmp_orchestra"
        mkdir -p "$TMPDIR"

        if [[ ! -f "$R_SCRIPT" ]]; then
            echo_log "ERROR: DESeq2 R script not found at $R_SCRIPT" | tee -a "$LOG_FILE"
            exit 1
        fi

        echo_log "Running DESeq2 differential abundance analysis with Rscript..."
        if ! Rscript "$R_SCRIPT" "$OTU_TABLE" "$METADATA" "$DESEQ2_OUTDIR" &> "$R_LOG"; then
            echo_log "❌ ERROR: DESeq2 analysis failed. Check R log: $R_LOG" | tee -a "$LOG_FILE"
            echo_log "Here are the last few lines of the R error log:" | tee -a "$LOG_FILE"
            tail -n 10 "$R_LOG" | tee -a "$LOG_FILE"
            exit 1
        fi

        echo_log "✅ DESeq2 analysis completed successfully. Results in: $DESEQ2_OUTDIR"
        ;;


    *)
        echo_log "Invalid step: $STEP" | tee -a "$LOG_FILE"
        exit 1
        ;;

    esac
done

END_TIME=$(date +"%Y-%m-%d %H:%M:%S")
echo_log "Pipeline completed at $END_TIME! Logs stored at $LOG_FILE" | tee -a $LOG_FILE
send_email "COMPLETED" "Metatranscriptomics Pipeline finished successfully at $END_TIME"

