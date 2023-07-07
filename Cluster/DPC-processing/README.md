# Processing differential phase contrast (DPC) images

The raw images from the Octopi included image pairs illuminated using the left half and the right half of the LED matrix. For each location, two images were captured. Each image pair was used to create a differential phase contrast (DPC) image, which resulted in an increase in the contrast of the processed image (as shown in the figure below).

![DPC](https://github.com/p-shrestha/erythroSight/assets/62575038/426a8601-f03b-4385-a1e0-89cbefe35f19)

Three main versions are provided to process data from the following three groups (Note: the programs for processing the Nepal/Canada data are the same except for the file paths, while the program for processing the time series data is different because the file/directory structure is different): 
- Nepal data: Images collected during the research study in Nepal (December 2022)
- Canada data: Images collected during the research study in Canada (August 2022 - April 2023)
- Canada time series: Time series images (vidoes) collected in Canada to observe the sickling process

The python script, e.g. `export_DPC_XY_Multiple_CC.py`, processes the raw images to DPC files by traversing through the directory/file structure. The bash shell script, e.g. `runDPCArchived.sh`, is used to submit the job in the cluster. The job is submitted using commands such as `sbatch runrunDPCArchived.sh`. Note that `runDPCCanadaTimeSeries.sh` uses a job array that reads the name of the unique IDs (UIDs) from the file `caseListTimeSeries`. 
