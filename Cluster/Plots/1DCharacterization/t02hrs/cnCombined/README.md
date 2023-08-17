# Combined data for Canada and Nepal



## morphParamCombine.m
The MATLAB file `morphParamCombine.m` combines the morphParam variables (saved in `morphParamData.mat`) and general information (saved in `generalInfo.mat`) according to the following: 
1. AA (Nepal) + N (Canada) = AA
2. ABeta (Nepal) = ABeta (no combining required)
3. AS (Nepal) + SCT (Canada) = AS
4. Seta (Nepal) + SS (Nepal) + SCD (Canada) = SCD

The shell script `runMPcombine.sh` runs the MATLAB program in the cluster, with all groups running sequentially in the single run (no parallelization), using the following script: 
```
sbatch runMPcombine.sh
```
