#!/bin/bash

#SBATCH --job-name=runMatrixEQTL
#SBATCH --account=vasccell
#SBATCH --nodes=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cs806@leicester.ac.uk
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --export=NONE

cd $SLURM_SUBMIT_DIR

module load R

# tnfa_pilot_BRBseq
Rscript brbSeq_scripts/autoPrepCovariates.R tnfa_pilot_BRBseq 2 3
Rscript brbSeq_scripts/autoRunMatrixEQTL_eigenMT.R tnfa_pilot_BRBseq 2 3
#Rscript brbSeq_scripts/autoRunMT.R tnfa_pilot_BRBseq 2 3

# tnfa_pilot_BRBseq
Rscript brbSeq_scripts/autoPrepCovariates.R untreated_pilot_BRBseq 2 3
Rscript brbSeq_scripts/autoRunMatrixEQTL_eigenMT.R untreated_pilot_BRBseq 2 3

# huvec50
Rscript brbSeq_scripts/autoPrepCovariates.R huvec50 2 3
Rscript brbSeq_scripts/autoRunMatrixEQTL_eigenMT.R huvec50 2 3

