#!/bin/bash
#SBATCH --job-name=featurecounts-pair-ended-combined   # Job name (same as the script name without the extension)
#SBATCH --output=featurecount-final.out                # Output file for standard output
#SBATCH --error=featurecount-final.err                 # Output file for standard error
#SBATCH --cpus-per-task=4                              # Number of CPU cores
#SBATCH --mem=25000                                    # Memory requirement (MB)
#SBATCH --time=13:00:00                                # Maximum runtime
#SBATCH --partition=pcoursea                           # Partition to submit the job to

#SBATCH --mail-user=nihui.shao@students.unibe.ch      
#SBATCH --mail-type=BEGIN,END,FAIL                     # Send email notifications at job start, end, and failure

# Set variables
ANNOTATION_GTF="/data/users/nshao/rnaseq/reference-genome/Homo_sapiens.GRCh37.87.gtf.gz"
THREADS=5
#INPUT_DIR=/data/users/nshao/rnaseq/mapping/hisat2_mapping/
#OUTPUT_DIR=/data/users/nshao/rnaseq/featurecounts/output_counts_paired2/
OUTPUT_COUNTS_PAIRED=/data/users/nshao/rnaseq/featurecounts/counts.txt
INPUT_BAM_SORTED=/data/users/nshao/rnaseq/mapping/hisat2_mapping/*sorted.bam

# Run featureCounts to generate a combined count matrix for all samples
apptainer exec --bind /data/ /containers/apptainer/subread_2.0.6.sif \
    featureCounts -p --countReadPairs -t exon -g gene_id -a $ANNOTATION_GTF -o $OUTPUT_COUNTS_PAIRED $INPUT_BAM_SORTED