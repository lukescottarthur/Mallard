#!/bin/bash
#SBATCH --job-name=pcadapt_bash_GFMxWM
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=12:00:00
#SBATCH --error=/home/las80898/Mallard/pcadapt_bash_GFMxWM.%j.err
#SBATCH --output=/home/las80898/Mallard/pcadapt_bash_GFMxWM.%j.out
#SBATCH --mail-user=las80898@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# activate R environment
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate test_env

#set output directory variable
OUTDIR="/scratch/las80898/pcadapt_output_3"                 

#if output directory doesn't exist, create it
if [ ! -d "$OUTDIR" ]
then
    mkdir -p "$OUTDIR"
fi

# run R script
R --no-save < /home/las80898/Mallard/pcadapt_bash_wholegenome_2.R