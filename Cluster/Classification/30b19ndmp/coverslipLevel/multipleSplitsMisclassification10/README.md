# Multiple splits with misclassification penalties

Classification is performed iteratively (e.g. 1000 times) by randomly splitting the data into 80:20 training:testing data split (participant-wise splits). The multiple iterations are performed to test the repeatability and determine confidence intervals for the classification results. 

Here, misclassification penalties (value of 10 instead of 1) for ABeta, AS or SCD misclassified as AA. 

The MATLAB program `classificationMultipleSplits_MisclassCost10_CS.m` is run in the cluster with the shell script `runClassificationMultipleSplits_MisClass10_CS.sh`, using the the following command for running 30 jobs (job array) to test the 30 different classifiers (saved as functions in separate MATLAB files): 
```
sbatch --time=30:00:00 --array=1-30 runClassificationMultipleSplits_MisClass10_CS.sh
```
