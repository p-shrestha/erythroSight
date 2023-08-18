# Normalizing some basic parameters such as area and perimeter

This section includes a MATLAB program to read the `morphParamData.mat` and `generalInfo.mat` files from each of the groups (AA, ABeta, AS, and SCD) and normalize some basic parameters. The selected basic parameters are 1. Area, 2. Perimeter, 4. Major, 5. Minor, 6. Height, 7. Width, 13. Feret, 15. MinFeret (8 parameters selected). The total number of morphological and intensity parameters is increased from 32 to 40, and saved. 

The MATLAB program `normalizeBasicParam.m` is run using the shell script `runNormalizeBMP.sh` by using the command: 
```
sbatch --time=02:00:00 --array=1-4 runNormalizeBMP.sh
```
