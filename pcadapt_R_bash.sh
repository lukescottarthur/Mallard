#!/bin/bash
#SBATCH --job-name=Revised_GFMxWMA_1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=4:00:00
#SBATCH --error=/home/las80898/Mallard/Revised_GFMxWMA_1.%j.err
#SBATCH --output=/home/las80898/Mallard/Revised_GFMxWMA_1.%j.out
#SBATCH --mail-user=las80898@uga.edu
#SBATCH --mail-type=START,END,FAIL

# activate R environment
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate test_env

#set output directory variable
OUTDIR="/scratch/las80898/pcadapt_output_5"                 

#if output directory doesn't exist, create it
if [ ! -d "$OUTDIR" ]
then
    mkdir -p "$OUTDIR"
fi

# run R script
R --no-save < /home/las80898/Mallard/ZZZ_pcadapt_GFMxWMA_revised.R