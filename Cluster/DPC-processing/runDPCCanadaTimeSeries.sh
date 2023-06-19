#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --array=1-14          # Run array of 14 simultaneous jobs in caseListTimeSeries
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-user=pranavsh@mail.ubc.ca

#source ~/ENV/bin/activate

echo "Starting task $SLURM_ARRAY_TASK_ID"
i=$(sed -n "${SLURM_ARRAY_TASK_ID}p" caseListTimeSeries)

module load gcc opencv python/3.10
virtualenv --no-download --clear $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate

#pip install --no-index -r requirementsDPC.txt


# Compute node: 1) Extract ScriptsCCScratch directory with the script export_DPC_XY_Multiple_CC.py and Logs directory
# 2) Create directory layers in the local directory $SLURM_TMPDIR, as per the code
cd $SLURM_TMPDIR
tar -xf /home/shrestha/scratch/ArchivedRAW/ScriptsCCScratch.tar
mkdir ImagesCanada
cd ImagesCanada
mkdir RAW
mkdir DPC
cd RAW
mkdir TimeSeries
cd ../DPC
mkdir TimeSeries

# ---------------------------------------------------------------------------------------------------------------------


# Login node: Create archived folder of raw images for particular UID
cd /home/shrestha/scratch/ImagesCanada/RAW/TimeSeries/
tar -cf /home/shrestha/scratch/ArchivedRAW/TimeSeries/$i-RAW.tar $i
echo "Archive created for RAW images for $i"

# Compute node: Extract .tar file in local storage
cd $SLURM_TMPDIR/ImagesCanada/RAW/TimeSeries/
tar -xf /home/shrestha/scratch/ArchivedRAW/TimeSeries/$i-RAW.tar
echo "Archive extracted for RAW images in local storage in compute node for $i"

# Now do my computations here on the local disk using the contents of the extracted archive...
# Compute node: Run the program using images extracted in the local storage
cd $SLURM_TMPDIR/ScriptsCCScratch/
srun python export_DPC_XY_Multiple_CC_Canada_TimeSeries.py --UID $i
echo "DPC processing completed for $i"

# The computations are done, so clean up the data set...
# Login node: Save the processed DPC images in the scratch space
cd $SLURM_TMPDIR/ImagesCanada/DPC/TimeSeries/
tar -cf /home/shrestha/scratch/ArchivedDPC/TimeSeries/$i-DPC.tar $i
echo "Archive of DPC images saved in ~/scratch/ArchivedDPC/ for $i"

# Login node: Save the Log files in the scratch space
cd $SLURM_TMPDIR/ScriptsCCScratch/Logs/CanadaExportDPC/
cp $i* /home/shrestha/scratch/ScriptsCCScratch/Logs/CanadaExportDPC/TimeSeries/
echo "Log file saved in ~/scratch/ScriptsCCScratch/Logs/CanadaExportDPC/ for $i"

# Login node: Extract the archived file in the scratch space
cd /home/shrestha/scratch/ImagesCanada/DPC/TimeSeries/
tar -xf /home/shrestha/scratch/ArchivedDPC/TimeSeries/$i-DPC.tar
echo "DPC images extracted and saved in ~/scratch/ImagesCanada/DPC/ for $i"

# Compute node: Once computations and archiving is completed, delete the files to free up space
cd $SLURM_TMPDIR/ImagesCanada/RAW/TimeSeries/
rm -r $i
echo "Copy of RAW images deleted from local storage in compute node for $i"
cd $SLURM_TMPDIR/ImagesCanada/DPC/TimeSeries/
rm -r $i
echo "Copy of DPC images deleted from local storage in compute node for $i"


# ---------------------------------------------------------------------------------------------------------------------
