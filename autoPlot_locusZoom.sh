#!/bin/bash

#SBATCH --job-name=plotExprGeno
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=40G
#SBATCH --time=23:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR


module load R




# Plot Expr vs Geno

Rscript brbSeq_scripts/autoPlot_locusZoom.R -c resFiles/ncsCands.txt
Rscript brbSeq_scripts/autoPlot_locusZoom.R -c resFiles/trsCands.txt




