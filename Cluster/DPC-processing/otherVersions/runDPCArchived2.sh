#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
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
#for i in AS-011 AS-012 AS-013 AS-014 AS-015 AS-016 AS-017 AS-018 AS-019 AS-020
#for i in ABeta-011 ABeta-012 ABeta-013 ABeta-014 ABeta-015 ABeta-016 ABeta-017 ABeta-018 ABeta-019 ABeta-020
#for i in AA-009 AA-010 ABeta-020 SBeta-009 SBeta010 AA-021 AA-022 AA-023 AA-024 ABeta-021 ABeta-022 ABeta-023 ABeta-024 AS-031 AS-032 AS-033 AS-034 AS-035 AS-036
for i in SBeta-010
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
