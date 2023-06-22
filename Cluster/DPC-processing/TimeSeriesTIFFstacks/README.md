# Creating TIFF stacks of the time series data

This code runs ImageJ headless in the cluster to combine the DPC images for the time series into image stacks and saves it as `.tif` files. 

The shell script `runTimeSeriesTIFFStack.sh` runs ImageJ in headless mode and runs the python script `createTIFFstacks.py`. Sample sbatch commands are provided below: 
```
sbatch --time=00:45:00 --array=1 runTimeSeriesTIFFStack.sh 
sbatch --time=01:30:00 --array=2-9 runTimeSeriesTIFFStack.sh 
sbatch --time=02:30:00 --array=10-12 runTimeSeriesTIFFStack.sh
sbatch --time=03:30:00 --array=13-14 runTimeSeriesTIFFStack.sh
sbatch --time=04:30:00 --array=15 runTimeSeriesTIFFStack.sh
```

After the TIFF stacks are created, 8 bit versions of the TIFF stacks can also be saved by running the `runTimeSeriesTIFFStack8bit.sh` shell script, which runs ImageJ with the `createTIFFstacks8bit.py` python code. The shell scripts can be run as follows: 
```
sbatch --time=00:30:00 --array=1 runTimeSeriesTIFFStack8bit.sh 
sbatch --time=00:25:00 --array=2-9 runTimeSeriesTIFFStack8bit.sh 
sbatch --time=00:35:00 --array=10-12 runTimeSeriesTIFFStack8bit.sh
sbatch --time=00:50:00 --array=13-15 runTimeSeriesTIFFStack8bit.sh
```
