#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=18:00:00
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-user=pranavsh@mail.ubc.ca

source ~/ENV/bin/activate

module load python
#pip install --no-index -r requirementsDPC.txt


# Compute node: 1) Extract ScriptsCCScratch directory with the script export_DPC_XY_Multiple_CC.py and Logs directory
# 2) Create directory layers in the local directory $SLURM_TMPDIR, as per the code
cd $SLURM_TMPDIR
tar -xf ~/scratch/ArchivedRAW/ScriptsCCScratch.tar
mkdir ImagesNepal
cd ImagesNepal
mkdir RAW
mkdir DPC

# ---------------------------------------------------------------------------------------------------------------------

#for i in SS-010 SS-011 SS-012 SS-013 SS-014 SS-015 SS-016 SS-017 SS-018
#for i in AA-002 AA-003 AA-004 AA-005 AA-006 AA-007 AA-008 AA-009 AA-010
#for i in AA-011 AA-012 AA-013 AA-014 AA-015 AA-016 AA-017 AA-018 AA-019 AA-020
#for i in SBeta-001 SBeta-002 SBeta-003 SBeta-004 SBeta-005 SBeta-006 SBeta-007 SBeta-008 SBeta-009 SBeta-010
for i in AS-021 AS-022 AS-023 AS-024 AS-025 AS-026 AS-027 AS-028 AS-029 AS-030 
do
  # Login node: Create archived folder of raw images for particular UID
  cd ~/scratch/ImagesNepal/RAW/
  tar -cf ~/scratch/ArchivedRAW/$i-RAW.tar $i
  echo "Archive created for RAW images for $i"

  # Compute node: Extract .tar file in local storage
  cd $SLURM_TMPDIR/ImagesNepal/RAW
  tar -xf ~/scratch/ArchivedRAW/$i-RAW.tar
  echo "Archive extracted for RAW images in local storage in compute node for $i"

  # Now do my computations here on the local disk using the contents of the extracted archive...
  # Compute node: Run the program using images extracted in the local storage
  cd $SLURM_TMPDIR/ScriptsCCScratch/
  srun python export_DPC_XY_Multiple_CC.py --UID $i
  echo "DPC processing completed for $i"

  # The computations are done, so clean up the data set...
  # Login node: Save the processed DPC images in the scratch space
  cd $SLURM_TMPDIR/ImagesNepal/DPC/
  tar -cf ~/scratch/ArchivedDPC/$i-DPC.tar $i
  echo "Archive of DPC images saved in ~/scratch/ArchivedDPC/ for $i"

  # Login node: Save the Log files in the scratch space
  cd $SLURM_TMPDIR/ScriptsCCScratch/Logs/NepalExportDPC/
  cp $i* ~/scratch/ScriptsCCScratch/Logs/NepalExportDPC/
  echo "Log file saved in ~/scratch/ScriptsCCScratch/Logs/NepalExportDPC/ for $i"

  # Login node: Extract the archived file in the scratch space
  cd ~/scratch/ImagesNepal/DPC/
  tar -xf ~/scratch/ArchivedDPC/$i-DPC.tar
  echo "DPC images extracted and saved in ~/scratch/ImagesNepal/DPC/ for $i"

  # Compute node: Once computations and archiving is completed, delete the files to free up space
  cd $SLURM_TMPDIR/ImagesNepal/RAW/
  rm -r $i
  echo "Copy of RAW images deleted from local storage in compute node for $i"
  cd $SLURM_TMPDIR/ImagesNepal/DPC/
  rm -r $i
  echo "Copy of DPC images deleted from local storage in compute node for $i"
done

# ---------------------------------------------------------------------------------------------------------------------
