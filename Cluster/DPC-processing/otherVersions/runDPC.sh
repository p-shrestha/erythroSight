#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=00:20:00
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-user=pranavsh@mail.ubc.ca

source ~/ENV/bin/activate

module load python
#pip install --no-index -r requirementsDPC.txt

srun python export_DPC_XY_Multiple_CC.py

