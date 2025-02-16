#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=1:30:00
#SBATCH --output=fastqc.out
#SBATCH --output=/home/nshao/data/courses/rnaseq/breastcancer/fastqc/output_%j.o
#SBATCH --error=/home/nshao/data/courses/rnaseq/breastcancer/fastqc/error_%j.e
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL


apptainer exec --bind /data/ /containers/apptainer/fastqc-0.12.1.sif fastqc \
  -o /home/nshao/data/courses/rnaseq/breastcancer/fastqc \
  /data/courses/rnaseq_course/breastcancer_de/reads/*.fastq.gz

#!/bin/bash
#SBATCH --job-name=rnaseq-qc
#SBATCH --output=job_output2.log
#SBATCH --error=job_error2.log
#SBATCH --partition=pibu_el8 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#singularity exec /containers/apptainer/fastqc-0.12.1.sif fastqc -o /data/users/kxia/project/fastqc input_file.fastq
#singularity pull docker://quay.io/biocontainers/hisat2:2.2.1--h87f3376_4
#SBATCH --mail-user=nihui.shao@students.unibe.ch     # set my mail 
#SBATCH --mail-type=BEGIN,END,FAIL             

# Step 2: Combine FastQC results using MultiQCdat
apptainer exec --bind /data/ /containers/apptainer/multiqc-1.19.sif multiqc -o /data/users/lyang/courses/rnaseq/breastcancer/fastqc /data/users/lyang/courses/rnaseq/breastcancer/fastqc

#mount the directory: /data/ to the container
#ues container image: /containers/apptainer/multiqc-1.19.sif
#the tool name from container: multiqc
#version: multiqc-1.19
#the output direvtor is: /data/users/nshao/rnaseq/1.fastqc/fastqc_summary
#the input directory is: /data/users/nshao/rnaseq/1.fastqc
#the input file is all the file  produced by fastqc