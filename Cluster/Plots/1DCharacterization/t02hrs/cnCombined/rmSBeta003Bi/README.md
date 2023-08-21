# Remove incorrect coverslip SBeta_003_Bi

This coverslip was wrongly mislabelled as SBeta, when it in fact was either AA or ABeta, likely due to mislabelling during data acquisition. 

The MATLAB program `rmSBeta003Bi.m` removes one coverslip from the `morphParamData.mat` file and fixes the `generalInfo.mat` file. The shell script `runRmSBeta003Bi.sh` is run using the following command (conservative time estimate provided): 
```
sbatch --time=02:00:00 runRmSBeta003Bi.sh
```
