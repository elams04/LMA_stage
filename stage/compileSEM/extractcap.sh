#!/bin/sh
#SBATCH -J ExCap
#SBATCH -p skylake
#SBATCH -n 1
#SBATCH -A b401
#SBATCH -t 04:30:00
#SBATCH -o ./%N.%x.out
#SBATCH -e ./%N.%x.err
# #SBATCH --mail-type=BEGIN,END
# #SBATCH --mail-user=fall@lma.cnrs-mrs.fr
#SBATCH --mem=85G 

# chargement des modules
module purge
module load userspace/all
python3/3.12.0

# echo of commands
set -x

# To compute in the submission directory
cd ${SLURM_SUBMIT_DIR}

# run mesh
srun python3 GetCapteurs.py > out.log
