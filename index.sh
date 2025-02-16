#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --partition=pshort_el8
#SBATCH --mem=8000
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=hisat2_indexing
#SBATCH --error=histat2build.err
#SBATCH --output=hisat2build.out
#SBATCH --mail-user=nihui.shao@students.unibe.ch     
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL,LANG=en_US.UTF-8,LC_ALL=en_US.UTF-8

# Define path variables and mkdir
OUTPUT_DIR="/data/users/nshao/rnaseq/2.index/hisat_build_out" 
IMAGE="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"
ASS="Homo_sapiens.GRCh38.dna.primary_assembly.fa" # the reference genome sequence downloaded from the Ensembl ftp site.


apptainer exec --bind /data:/data "$IMAGE" hisat2-build -p 4 "/data/users/nshao/rnaseq/2.index/reference_genome/$ASS" "/data/users/nshao/rnaseq/2.index/hisat_build_out/hisat_build_out"


