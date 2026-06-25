#!/bin/bash
#SBATCH --job-name=cmplot_figures
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8gb
#SBATCH --time=1:00:00
#SBATCH --error=/home/las80898/Mallard/cmplot_figures.%j.err
#SBATCH --output=/home/las80898/Mallard/cmplot_figures.%j.out
#SBATCH --mail-user=las80898@uga.edu
#SBATCH --mail-type=START,END,FAIL

# activate R environment
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate test_env

#set output directory variable
OUTDIR="/scratch/las80898/mallard_cmplot"

#if output directory doesn't exist, create it
if [ ! -d "$OUTDIR" ]
then
    mkdir -p "$OUTDIR"
fi

cd $OUTDIR

# run R script
R --no-save < /home/las80898/Mallard/ZZZ_cmplot_figures.R