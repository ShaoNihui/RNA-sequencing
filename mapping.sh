#!/bin/bash
#SBATCH --job-name=mapping-job          
#SBATCH --cpus-per-task=4              
#SBATCH --mem=25000                    
#SBATCH --time=48:00:00            
#SBATCH --partition=pcmpg_el8        
#SBATCH --error=/data/users/nshao/rnaseq/outputs/histat2-mapping.err
#SBATCH --output=/data/users/nshao/rnaseq/outputs/hisat2-mapping.out

#SBATCH --mail-user=nihui.shao@students.unibe.ch     
#SBATCH --mail-type=BEGIN,END,FAIL



# Set variables
REF_GENOME="/data/users/nshao/rnaseq/2.index/reference_genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"  # the reference genome
GTF_FILE="/data/users/nshao/rnaseq/2.index/reference_genome/Homo_sapiens.GRCh38.113.gtf"  # the gif file of reference genome
INDEX_BASE="/data/users/nshao/rnaseq/2.index/hisat_build_out/hisat_build_out"  # index base 

#use for loop to loop over every sample 
for FASTQ_1 in /data/courses/rnaseq_course/breastcancer_de/reads/*_R1.fastq.gz 
do
    # remove tail to get fastq_2
    SAMPLE=$(basename $FASTQ_1 1.fastq.gz)
    
    # set the file name for mate2
    FASTQ_2="/data/courses/rnaseq_course/breastcancer_de/reads/${SAMPLE}2.fastq.gz"
    
    # set the output directory
    OUTPUT_SAM="/data/users/nshao/rnaseq/mapping/hisat2_mapping/${SAMPLE}_output.sam"
    OUTPUT_BAM="/data/users/nshao/rnaseq/mapping/hisat2_mapping/${SAMPLE}_output.bam"
    SORTED_BAM="/data/users/nshao/rnaseq/mapping/hisat2_mapping/${SAMPLE}_sorted.bam"
    OUTPUT_LOG="/data/users/nshao/rnaseq/mapping/hisat2_mapping_log/${SAMPLE}_hisat2.log"
    THREADS=4  

    # Step 1: Run Hisat2 to map the reads to the reference genome with RF strandedness
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
      hisat2 -x $INDEX_BASE -1 $FASTQ_1 -2 $FASTQ_2 -S $OUTPUT_SAM -p $THREADS  >> $OUTPUT_LOG
      #the samples are unstranded
      #save output log

    # Step 2: Convert SAM to BAM format using Samtools
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
      samtools view -hbS $OUTPUT_SAM > $OUTPUT_BAM

    # Step 3: Sort the BAM file by genomic coordinates using Samtools
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
      samtools sort -m 2G -@ $THREADS -o $SORTED_BAM -T temp $OUTPUT_BAM

    # Step 4: Index the sorted BAM file using Samtools
    apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
      samtools index $SORTED_BAM
done


#!/bin/bash
#SBATCH --job-name=multiqc-job
#SBATCH --cpus-per-task=4
#SBATCH --mem=8000
#SBATCH --time=04:00:00
#SBATCH --partition=pibu_el8
#SBATCH --error=/data/users/nshao/rnaseq/3.mapping/multiqc.err
#SBATCH --output=/data/users/nshao/rnaseq/3.mapping/multiqc.out

#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

# Set input and output directories
INPUT_DIR="/data/users/nshao/rnaseq/3.mapping/hisat2_mapping_log"  # Directory containing analysis results (e.g., logs, BAM files)
OUTPUT_DIR="/data/users/nshao/rnaseq/3.mapping/multiqc_report"        # Directory to store the MultiQC report
CONTAINER_IMAGE="/containers/apptainer/multiqc-1.19.sif"   # MultiQC container image

# Load MultiQC using the container
apptainer exec --bind /data/ $CONTAINER_IMAGE \
  multiqc $INPUT_DIR -o $OUTPUT_DIR

