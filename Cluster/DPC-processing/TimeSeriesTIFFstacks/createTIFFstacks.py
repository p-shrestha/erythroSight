#@String UID
#@String dataPath
#@String resultsPath
from ij import IJ
from os import listdir
from os.path import join, isfile
import os
import time
import re

t = time.time() # Record current time to calculate elapsed time

# Path and UID values, when arguments are not used
#UID = ["SCD-001"]
#dataPath = "/home/shrestha/scratch/ImagesCanada/DPC/TimeSeries/"
#resultsPath = "/home/shrestha/scratch/ImagesCanada/DPC/TimeSeries/TIFFstacks/"


# Note: Here UID is entered as a string and not a list of strings, so do not need to loop through
dataDirUID = dataPath + UID + "/"
resultsDirUID = resultsPath + UID + "/"

# Make folder in results folder with UID name if not there
if not os.path.isdir(resultsDirUID):
	os.makedirs(resultsDirUID)
	#print('Created directory resultsDirUID =' + resultsDirUID)

#dirInDataUID = os.listdir(dataDirUID)  # Lists the directories or folders in the UID folder of the data folder

#tname = [indx for indx in dirInDataUID if indx[0].lower() == 't'] # Save all directory names starting with t in the UID data folder (Using lower(), just in case files were saved using 'T' instead of 't')

CSIDfull = os.listdir(dataDirUID)  # Lists the directories or folders in the tname folder
CSID = [x.split('_')[0] for x in CSIDfull]  # Create a list containing only the coverslip ID or CSID, such as Ai, Aii, Bi, Bii, and so on


# Iterate through the folders to open DPC images
for fInd in range(len(CSIDfull)):
	#if not os.path.isdir(resultsDirUID + CSIDfull[fInd]):
		#os.makedirs(resultsDirUID + CSIDfull[fInd])

	filepath = dataDirUID + CSIDfull[fInd] + '/'
	fileList = [f for f in listdir(filepath) if isfile(join(filepath, f))]

	file_name0 = fileList[0]
	# Find the left and right parts of the string around the image number
	match = re.split(r'i_(\d+)_dt', file_name0)

	if len(match) > 1:
		leftS = match[0]+ 'i_'
		rightS = '_dt' + match[2]
		print("Left string:" + leftS)
		print("Right string:" +  rightS)
	else:
		print("Splitting pattern not found.")

	for fLInd in range(len(fileList)):

		t1 = time.time() # Record current time to calculate elapsed time

		#file_name = fileList[fLInd]
		file_name = leftS + str(fLInd) + rightS
		print("UID =" + UID + " ; CSIDfull = " + CSIDfull[fInd])
		print("DPC .png file name = " + file_name)

		if os.path.isfile(filepath + file_name):
			# print("File =" + file_name + " found in " + filepath)
			# Open DPC image corresponding to the filename
			macroImage = "open(\"" + filepath + file_name + "\");"
			IJ.runMacro(macroImage)  # Run macro to open DPC image

		else:
			print("File =" + file_name + " was not found in " + filepath)

		# Calculate and print elapsed time
		elapsedTime1 = time.time() - t1
		strTime1 = "{:.3f}".format(elapsedTime1)
		strET1 = 'Elapsed time =' + strTime1 + 'seconds'
		print(strET1)
	# Combine all the opened images into a stack and save as TIFF
	macroStack = """
	run("Images to Stack", "use");
	run("Set Scale...", "distance=3000 known=900 unit=um");
	"""
	IJ.runMacro(macroStack)  # Run macro to combine all images into stack

	file_nameSave = leftS + rightS.rsplit(".", 1)[0] + ".tif"
	macroSave = "saveAs(\"Tiff\", \"" + resultsDirUID + file_nameSave + "\");"
	IJ.runMacro(macroSave)  # Run macro to save file

	macroClose = "close();"
	IJ.runMacro(macroClose)  # Run macro to close stack

# Calculate and print elapsed time
elapsedTime = time.time() - t
strTime = "{:.3f}".format(elapsedTime)
strET = 'Total elapsed time for complete run =' + strTime + 'seconds'
print(strET)
