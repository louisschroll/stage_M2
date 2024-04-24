#!/bin/sh
#SBATCH --job-name=nmix_1
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=cefe

### Email
#SBATCH --mail-user=louis.schroll@cefe.cnrs.fr         
#SBATCH --mail-type=ALL

### Output & Error
#SBATCH --output=/lustre/schrolll/3.results/nmix_1.out
#SBATCH --error=/lustre/schrolll/3.results/nmix_1.err              


echo « Running on: $SLURM_NODELIST »

##load R
module load singularity

singularity exec R-4.3.2-equipe-HAIR.img Rscript /lustre/schrolll/2.code/pt2_telemetry_and_count/nmix_all_models.R