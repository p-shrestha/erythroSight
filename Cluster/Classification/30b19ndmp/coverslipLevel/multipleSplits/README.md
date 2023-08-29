# Multiple training:testing splits

Classification is performed iteratively (e.g. 1000 times) by randomly splitting the data into 80:20 training:testing data split (participant-wise splits). 

The MATLAB program `classificationMultipleSplits_CS.m` is run in the cluster with the shell script `runClassificationMultipleSplits_CS.sh`, using the the following command for running 30 jobs (job array) to test the 30 different classifiers (saved as functions in separate MATLAB files): 
```
sbatch --time=30:00:00 --array=1-30 runClassificationMultipleSplits_CS.sh
```
