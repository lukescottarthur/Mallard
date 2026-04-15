#!/bin/bash
#SBATCH --job-name=splitvcfmallard                      # Job name 
#SBATCH --partition=batch                           # Partition name 
#SBATCH --ntasks=1                                  # 
#SBATCH --cpus-per-task=8                           # CPU core count per task
#SBATCH --mem=32G                                    # Memory per node
#SBATCH --time=4:00:00                              # Time limit hrs:mins:secs
#SBATCH --output=/home/las80898/Mallard/%x_%j.out  ####change to your path
#SBATCH --error=/home/las80898/Mallard/%x_%j.error ####change to your path
#SBATCH --mail-user=las80898@uga.edu                #####change to your path
#SBATCH --mail-type=END,FAIL                        # Mail events 

# Load conda environment 
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate test_env

# Navigate to directory (change to your directory)
cd /home/las80898/mallard_wholegenome_data


## split into three files
# GFM x WM
bcftools view -Oz -S /home/las80898/Mallard/GFMxWM_samplelist.txt Wholegenomes_Allpopulations_reheader.vcf.gz -o GFMxWM.vcf.gz

# GFM x KC
bcftools view -Oz -S /home/las80898/Mallard/GFMxKC_samplelist.txt Wholegenomes_Allpopulations_reheader.vcf.gz -o GFMxKC.vcf.gz

# WM x KC
bcftools view -Oz -S /home/las80898/Mallard/WMxKC_samplelist.txt Wholegenomes_Allpopulations_reheader.vcf.gz -o WMxKC.vcf.gz

# make index files
bcftools index GFMxWM.vcf.gz
bcftools index GFMxKC.vcf.gz
bcftools index WMxKC.vcf.gz