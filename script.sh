# # Directory that stores the Fastq file downloaded
# DATA_DIR="/home/eaderogba279/Bionformatics_Projects/TNBC_Lespedeza_RNASeq_Analysis/data"
genome_index="grch38_genome/grch38/genome"

# # Directory that stores the FastQC reports
# REPORT_DIR="/home/eaderogba279/Bionformatics_Projects/TNBC_Lespedeza_RNASeq_Analysis/fastp_output"

# # List of accession numbers
acc_no=("SRR28420795" "SRR28420796" "SRR28420797" "SRR28420798")

# ### FASTQC ###

# # Looping through the array of accession numbers
# for SRA in "${acc_no[@]}"
# do
#     echo "Downloading $acc_no..."
#     # Use fastq-dump with --split-files to download paired-end reads
#     fastq-dump --split-files "$acc_no" -O $DATA_DIR
#     echo "$acc_no download complete!"
# done

# echo "All downloads are completed"

# # ### GENERATING MULTIQC REPORT ###

# cd $REPORT_DIR
# echo "Generating MultiQC Report"

# multiqc .

# echo "MultiQC Report Generated"

# ### FASTP ###

# for FILE in $DATA_DIR/*_1.fastq
# do
#     # Extract the basename without the _1.fastq.gz or _2.fastq.gz extension
#     BASE_NAME=$(basename $FILE _1.fastq)

#     # Define the names of the forward and the reverse read files
#     FORWARD_READS="${DATA_DIR}/${BASE_NAME}_1.fastq"
#     REVERSE_READS="${DATA_DIR}/${BASE_NAME}_2.fastq"

#     # Check if the corresponding reverse read file exists
#     if [ -f $REVERSE_READS ]; then
#         echo "Running Fastp on paired-end reads: $FORWARD_READS and $REVERSE_READS..."
#         fastp -i $FORWARD_READS -I $REVERSE_READS -o "${REPORT_DIR}/${BASE_NAME}_1.fastq" -O "${REPORT_DIR}/${BASE_NAME}_2.fastq"

#         echo "Fastp processing completed for $BASE_NAME"
#     else
#         echo "Reverse read file for $BASE_NAME does not exist. Skipping..."
#     fi
# done

# echo "FastQC analysis completed for all paired-end files!"

#### MAPPING TO REFERENCE SEQUENCE ####
#hisat2 -x [path_to_index] -1 [forward_reads.fq] -2 [reverse_reads.fq] -S [output.sam]

# for sample in "${acc_no[@]}"; do
#     # Input files for the paired-end reads
#     fastq1="data/${sample}_1.fastq"
#     fastq2="data/${sample}_2.fastq"

#     # Output SAM file path
#     output_sam="alignment/${sample}.sam"

#     # Run HISAT2 and redirect output to SAM file
#     hisat2 -x "$genome_index" -p 7 -1 "$fastq1" -2 "$fastq2" -S "$output_sam"

#     # Sort the SAM file and convert to BAM
#     samtools sort -o "${output_sam%.sam}.sorted.bam" "${output_sam}"
# done

# run feature counts
# featureCounts -T 7 -s 2 -p -a data/Homo_sapiens.GRCh38.111.gtf -o counts/featureCounts.txt alignment/SRR28420795.sorted.bam alignment/SRR28420796.sorted.bam alignment/SRR28420797.sorted.bam alignment/SRR28420798.sorted.bam

# Cleaning up the featureCounts matrix

#cut -f1,7,8,9,10,11,12 counts/featureCounts.txt > counts/featureCounts.Rmatrix.txt
cat counts/featureCounts.txt | sed "2s/alignment\///g" | sed "2s/\.sorted.bam//g" > counts/featureCounts_results.mod.txt
