#!/bin/bash
#SBATCH --account=def-stoeber
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --mail-type=begin         # send mail when job begins
#SBATCH --mail-type=end           # send mail when job ends
#SBATCH --mail-type=fail          # send mail if job fails
#SBATCH --mail-type=array_tasks   # send mail for all job array tasks
#SBATCH --mail-user=pranavsh@mail.ubc.ca

echo "Starting task $SLURM_ARRAY_TASK_ID"

# 1) Compute Node: Create virtual environment (node-local virtual environment) inside the job
module load gcc opencv python/3.10
virtualenv --no-download --clear $SLURM_TMPDIR/env && source $SLURM_TMPDIR/env/bin/activate

# 2) Compute Node (virtual): Activate matlab
module load matlab/2023a

# 3) Define array of classifier names
cNames=(
"1FineTree"
"2MediumTree"
"3CoarseTree"
"6EffLogReg"
"7EffLinSVM"
"9KernelNaiveBayes"
"10LinSVM"
"11QuadSVM"
"12CubicSVM"
"13FineGaussSVM"
"14MedGaussSVM"
"15CoarseGaussSVM"
"16FineKNN"
"17MedKNN"
"18CoarseKNN"
"19CosineKNN"
"20CubicKNN"
"21WeightedKNN"
"22BoostedTree"
"23BaggedTree"
"24SubspaceDis"
"25SubspaceKNN"
"26RUSBoostedTree"
"27NarrowNeuralNet"
"28MedNeuralNet"
"29WideNeuralNet"
"30BilayeredNeuralNet"
"31TrilayeredNeuralNet"
"32SVMKernel"
"33LogisticRegKernel"
)

# 4) Compute Node (virtual): Run MATLAB program
matlab -nodisplay -r "classifierName = \"${cNames[$SLURM_ARRAY_TASK_ID - 1]}\";classificationMultipleSplits_MisclassCost10_CS"
