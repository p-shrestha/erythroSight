#@String UID
#@String dataPath
#@String resultsPath
from ij import IJ
from ij.gui import PolygonRoi, Roi, Overlay
from ij.measure import ResultsTable, Measurements
from java.awt import FileDialog
from os import walk, listdir
from os.path import join, isfile
import os
import math
import argparse
import time

t = time.time() # Record current time to calculate elapsed time

# Open blank image with the dimensions of a processed DPC image; Note: the path needs to be updated based on where this program is running
macrotest1 = """
open("/home/shrestha/projects/def-stoeber/shrestha/MorphologicalParameters/Scripts/Blank_Image.tif");
selectWindow("Blank_Image.tif");
"""

IJ.runMacro(macrotest1)	# Run macro to open blank image

# Path and UID values, when arguments are not used
#UID = ["Test"]
#dataPath = "/home/shrestha/projects/def-stoeber/shrestha/Cellpose/OutlinesCanada/"
#resultsPath = "/home/shrestha/projects/def-stoeber/shrestha/MorphologicalParameters/Canada/"


# Note: Here UID is entered as a string and not a list of strings, so do not need to loop through
dataDirUID = dataPath + UID + "/"
resultsDirUID = resultsPath + UID + "/"

# Make folder in results folder with UID name if not there
if not os.path.isdir(resultsDirUID):
	os.makedirs(resultsDirUID)
	#print('Created directory resultsDirUID =' + resultsDirUID)

dirInDataUID = os.listdir(dataDirUID)  # Lists the directories or folders in the UID folder of the data folder

tname = [indx for indx in dirInDataUID if indx[0].lower() == 't'] # Save all directory names starting with t in the UID data folder (Using lower(), just in case files were saved using 'T' instead of 't')

for tnameInd in range(len(tname)): # Loop through all folders starting with t, such as t0, t2hrs, and so on
	CSIDfull = os.listdir(dataDirUID + tname[tnameInd] + '/')  # Lists the directories or folders in the tname folder
	CSID = [x.split('_')[0] for x in CSIDfull]  # Create a list containing only the coverslip ID or CSID, such as Ai, Aii, Bi, Bii, and so on

	# Make directory with tname in results
	if not os.path.isdir(resultsDirUID + tname[tnameInd]):
		os.makedirs(resultsDirUID + tname[tnameInd])
		#print('Created directory resultsDirUID + tname[tnameInd] =' + resultsDirUID + tname[tnameInd])

	# Iterate through the folders to read outline files and save morphological parameters
	for fInd in range(len(CSIDfull)):
		if not os.path.isdir(resultsDirUID + tname[tnameInd] + '/' + CSIDfull[fInd]):
			os.makedirs(resultsDirUID + tname[tnameInd] + '/' + CSIDfull[fInd])
			#print('Created directory resultsDirUID + tname[tnameInd] + CSIDfull[fInd] ='+ resultsDirUID + tname[tnameInd] + '/' + CSIDfull[fInd])

		filepath = dataDirUID + tname[tnameInd] + '/' + CSIDfull[fInd] + '/'
		fileList = [f for f in listdir(filepath) if isfile(join(filepath, f))]

		for fLInd in range(len(fileList)):

			t1 = time.time() # Record current time to calculate elapsed time

			file_name = fileList[fLInd]
			print("UID =" + UID + " ; tname = " + tname[tnameInd] + " ; CSIDfull = " + CSIDfull[fInd])
			print("Outlines .txt file name = " + file_name)

			tPosition = file_name.index('_t')
			IDNum = file_name[0:-20]

			# Open DPC image corresponding to the filename
			macroImage = "open(\"/home/shrestha/scratch/ImagesCanada/DPC/" + UID + "/" + tname[tnameInd] + "/" + CSIDfull[fInd] + "/" + IDNum +  "_DPC.png\");"
			IJ.runMacro(macroImage) # Run macro to open DPC image

			# Macro to select window
			macroSelect = "selectWindow(\""+ IDNum +  "_DPC.png\");"
			IJ.runMacro(macroSelect) # Run macro to select DPC image

			imp = IJ.getImage()

			overlay = imp.getOverlay()
			overlay = Overlay()

			textfile = open(filepath + file_name, "r")
			for line in textfile:
			    xy = map(int, line.rstrip().split(","))
			    X = xy[::2]
			    Y = xy[1::2]
			    imp.setRoi(PolygonRoi(X, Y, Roi.POLYGON))
			    roi = imp.getRoi()
			    imp.setRoi(roi)
			    overlay.add(roi)

			    # Note: The methods available like getStatistics() and getFeretsDiameter() can be found in link:
			    # https://imagej.nih.gov/ij/developer/api/ij/ij/gui/Roi.html

			    # NOTE: The list of parameters that work with roi.getStatistics() is in the following link:
			    # https://imagej.nih.gov/ij/developer/api/ij/ij/process/ImageStatistics.html

			    # Add measurements to the ResultsTable
			    rt = ResultsTable.getResultsTable()
			    rt.incrementCounter()
			    rt.addValue("Area", roi.getStatistics().area)
			    rt.addValue("Perimeter", roi.getLength())
			    rt.addValue("Angle", roi.getStatistics().angle)
			    rt.addValue("Major", roi.getStatistics().major)
			    rt.addValue("Minor", roi.getStatistics().minor)
			    rt.addValue("Height", roi.getStatistics().roiHeight)
			    rt.addValue("Width", roi.getStatistics().roiWidth)

			    rt.addValue("Mean",roi.getStatistics().mean)
			    rt.addValue("StdDev",roi.getStatistics().stdDev)
			    rt.addValue("Median",roi.getStatistics().median)
			    rt.addValue("Skewness",roi.getStatistics().skewness)
			    rt.addValue("Kurtosis",roi.getStatistics().kurtosis)

			    rt.addValue("Feret", roi.getFeretValues()[0])
			    rt.addValue("FeretAngle", roi.getFeretValues()[1])
			    rt.addValue("MinFeret", roi.getFeretValues()[2])
			    rt.addValue("FeretX", roi.getFeretValues()[3])
			    rt.addValue("FeretY", roi.getFeretValues()[4])

			    # Calculated values
			    circularity = (4 * math.pi * roi.getStatistics().area) / (roi.getLength() ** 2)
			    rt.addValue("Circularity", circularity)

			    AR = roi.getStatistics().major/roi.getStatistics().minor
			    rt.addValue("AR", AR)

			    Roundness = (4 * roi.getStatistics().area) / ( math.pi *(roi.getStatistics().major ** 2))
			    rt.addValue("Roundness", Roundness)

			    ConvexHullArea = PolygonRoi(roi.getConvexHull().xpoints, roi.getConvexHull().ypoints, roi.getConvexHull().npoints, PolygonRoi.POLYGON).getStatistics().area
			    rt.addValue("ConvexHullArea", ConvexHullArea)

			    Solidity = roi.getStatistics().area/ConvexHullArea
			    rt.addValue("Solidity", Solidity)

			    Eccentricity = math.sqrt(1-((0.5*(roi.getStatistics().minor ** 2))/(0.5*(roi.getStatistics().major ** 2))))
			    rt.addValue("Eccentricity", Eccentricity)

			    ESF = 1/AR
			    rt.addValue("ESF", ESF)

			    # Some of the following measurements are adapted from "Shape Analysis & Measurement", Michael A. Wirth, 2004
			    ElongationFW = (roi.getStatistics().roiWidth)/(roi.getFeretValues()[0]) # Calculated here as width/Feret's diameter
			    rt.addValue("ElongationFW", ElongationFW)

			    ElongationFmF = (roi.getFeretValues()[2])/(roi.getFeretValues()[0]) # Calculated here as Minimum Feret / Feret's diameter and minimum Feret (Similar to sphericity?)
			    rt.addValue("ElongationFmF", ElongationFmF)

			    ConvexPerimeter = PolygonRoi(roi.getConvexHull().xpoints, roi.getConvexHull().ypoints, roi.getConvexHull().npoints, PolygonRoi.POLYGON).getLength()
			    rt.addValue("ConvexPerimeter", ConvexPerimeter)

			    Convexity = ConvexPerimeter/roi.getLength()
			    rt.addValue("Convexity", Convexity)

			    FibreLength = (roi.getLength()-math.sqrt(abs((roi.getLength() **2)-(16*roi.getStatistics().area))))/4
			    rt.addValue("FibreLength", FibreLength)

			    FibreWidth = roi.getStatistics().area/FibreLength
			    rt.addValue("FibreWidth", FibreWidth)

			    Curl = roi.getFeretValues()[0]/FibreLength
			    rt.addValue("Curl",Curl)

			    CurlW = roi.getFeretValues()[0]/FibreWidth
			    rt.addValue("CurlW",CurlW)

			textfile.close()

			imp.setOverlay(overlay)

			if not os.path.isdir(resultsDirUID + tname[tnameInd] + '/' + CSIDfull[fInd]):
			    os.makedirs(resultsDirUID + tname[tnameInd] + '/' + CSIDfull[fInd])
			    print("Now here: Created directory"+ resultsDirUID + tname[tnameInd] + '/' + CSIDfull[fInd])

			#savefile = resultsPath + UID + "/ShapeParameters_" + UID + "-" + IDNum + ".csv"
			savefile = resultsDirUID + tname[tnameInd] + '/' + CSIDfull[fInd] + "/ShapeParameters_" + UID + "-" + IDNum + ".csv"

			rt.saveAs(savefile)

			IJ.run("Clear Results")

			# Macro to select window
                        macroClose = "close(\""+ IDNum +  "_DPC.png\");"
                        IJ.runMacro(macroClose) # Run macro to close DPC image

			# Calculate and print elapsed time
			elapsedTime1 = time.time() - t1
			strTime1 = "{:.3f}".format(elapsedTime1)
			strET1 = 'Elapsed time =' + strTime1 + 'seconds'
			print(strET1)

# Calculate and print elapsed time
elapsedTime = time.time() - t
strTime = "{:.3f}".format(elapsedTime)
strET = 'Total elapsed time for complete run =' + strTime + 'seconds'
print(strET)
