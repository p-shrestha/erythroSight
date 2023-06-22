#@String UID
#@String dataPath
#@String resultsPath
from ij import IJ
from os import listdir
from os.path import join, isfile
import os
import time

t = time.time() # Record current time to calculate elapsed time

# Path and UID values, when arguments are not used
#UID = ["SCD-001"]
#dataPath = "/home/shrestha/scratch/ImagesCanada/DPC/TimeSeries/TIFFstacks/"
#resultsPath = "/home/shrestha/scratch/ImagesCanada/DPC/TimeSeries/TIFFstacks8bit/“

# Note: Here UID is entered as a string and not a list of strings, so do not need to loop through
dataDirUID = dataPath + UID + "/"
resultsDirUID = resultsPath + UID + "/"

# Make folder in results folder with UID name if not there
if not os.path.isdir(resultsDirUID):
	os.makedirs(resultsDirUID)
	#print('Created directory resultsDirUID =' + resultsDirUID)

stackList = os.listdir(dataDirUID) # Lists all the files in data folder

# Iterate through all the stacks in the UID
for fInd in range(len(stackList)):
	t1 = time.time()  # Record current time to calculate elapsed time
	filepath = dataDirUID
	file_name = stackList[fInd]
	if os.path.isfile(filepath + file_name):
		# print("File =" + file_name + " found in " + filepath)
		# Open TIFF stack corresponding to file_name
		macroImage = "open(\"" + filepath + file_name + "\");"
		IJ.runMacro(macroImage)  # Run macro to open DPC image

	else:
		print("TIFF stack file =" + file_name + " was not found in " + filepath)

	# Convert stack to 8 bit
	macro8bit = "run(\"8-bit\");"
	IJ.runMacro(macro8bit)  # Run macro to convert TIFF stack to 8 bit

	file_nameSave = file_name.rsplit(".", 1)[0] + "8bit.tif"
	macroSave = "saveAs(\"Tiff\", \"" + resultsDirUID + file_nameSave + "\");"
	IJ.runMacro(macroSave)  # Run macro to save file
	print('TIFF file saved')

	macroClose = "close();"
	IJ.runMacro(macroClose)  # Run macro to close stack

	elapsedTime1 = time.time() - t1
	strTime1 = "{:.3f}".format(elapsedTime1)
	strET1 = 'Elapsed time =' + strTime1 + 'seconds'
	print(strET1)

# Calculate and print elapsed time
elapsedTime = time.time() - t
strTime = "{:.3f}".format(elapsedTime)
strET = 'Total elapsed time for complete run =' + strTime + 'seconds'
print(strET)
