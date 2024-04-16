#!/bin/sh
#SBATCH --job-name=sterne_goeland_model
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --mem-per-cpu=5G
#SBATCH --partition=cefe

### Email
#SBATCH --mail-user=louis.schroll@cefe.cnrs.fr         
#SBATCH --mail-type=ALL

### Output & Error
#SBATCH --output=/lustre/schrolll/3.results/sterne_goel_model.out
#SBATCH --error=/lustre/schrolll/3.results/sterne_goel_model.err              


echo « Running on: $SLURM_NODELIST »

##load R
module load singularity

singularity exec R-4.3.2-equipe-HAIR.img Rscript /lustre/schrolll/2.code/pt2_telemetry_and_count/all_models_telemetry_count.R