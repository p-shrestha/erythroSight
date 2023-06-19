#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-user=pranavsh@mail.ubc.ca

#source ~/ENV/bin/activate

#module load python
#pip install --no-index -r requirementsDPC.txt

module load gcc opencv python/3.10
virtualenv --no-download --clear $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate

# Compute node: 1) Extract ScriptsCCScratch directory with the script export_DPC_XY_Multiple_CC.py and Logs directory
# 2) Create directory layers in the local directory $SLURM_TMPDIR, as per the code
cd $SLURM_TMPDIR
tar -xf ~/scratch/ArchivedRAW/ScriptsCCScratch.tar
mkdir ImagesCanada
cd ImagesCanada
mkdir RAW
mkdir DPC

# ---------------------------------------------------------------------------------------------------------------------

for i in N-005
do
  # Login node: Create archived folder of raw images for particular UID
#  cd ~/scratch/ImagesCanada/RAW/RatioStudy/
#  tar -cf ~/scratch/ArchivedRAW/RS-$i-RAW.tar $i
#  echo "Archive created for RAW images for $i"

  # Compute node: Extract .tar file in local storage
  cd $SLURM_TMPDIR/ImagesCanada/RAW
  tar -xf ~/scratch/ArchivedRAW/RS-$i-RAW.tar
  echo "Archive extracted for RAW images in local storage in compute node for $i"

  # Now do my computations here on the local disk using the contents of the extracted archive...
  # Compute node: Run the program using images extracted in the local storage
  cd $SLURM_TMPDIR/ScriptsCCScratch/
  srun python export_DPC_XY_Multiple_CC_Canada.py --UID $i
  echo "DPC processing completed for $i"

  # The computations are done, so clean up the data set...
  # Login node: Save the processed DPC images in the scratch space
  cd $SLURM_TMPDIR/ImagesCanada/DPC/
  tar -cf ~/scratch/ArchivedDPC/RS-$i-DPC.tar $i
  echo "Archive of DPC images saved in ~/scratch/ArchivedDPC/ for $i"

  # Login node: Save the Log files in the scratch space
  cd $SLURM_TMPDIR/ScriptsCCScratch/Logs/CanadaExportDPC/
  cp $i* ~/scratch/ScriptsCCScratch/Logs/CanadaExportDPC/
  echo "Log file saved in ~/scratch/ScriptsCCScratch/Logs/CanadaExportDPC/ for $i"

  # Login node: Extract the archived file in the scratch space
  cd ~/scratch/ImagesCanada/DPC/RatioStudy/
  tar -xf ~/scratch/ArchivedDPC/RS-$i-DPC.tar
  echo "DPC images extracted and saved in ~/scratch/ImagesCanada/DPC/ for $i"

  # Compute node: Once computations and archiving is completed, delete the files to free up space
  cd $SLURM_TMPDIR/ImagesCanada/RAW/
  rm -r $i
  echo "Copy of RAW images deleted from local storage in compute node for $i"
  cd $SLURM_TMPDIR/ImagesCanada/DPC/
  rm -r $i
  echo "Copy of DPC images deleted from local storage in compute node for $i"
done

# ---------------------------------------------------------------------------------------------------------------------
