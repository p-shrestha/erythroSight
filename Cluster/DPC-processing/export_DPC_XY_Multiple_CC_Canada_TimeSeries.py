# Script to process DPC images based on Left and Right half LED illuminated images implemented in Compute Canada
# Script to analyze data from Canada
# Based on script from H. Li (2022)
# 2023 ps ubc

import time
import cv2
import os
import numpy as np
import json
from datetime import datetime
import pandas as pd
import argparse

# Script saved in folder ('ScriptsCC' - can be named anything), under the root directory with another directory "ImagesCanada" with all the images
# Folder structure in RAW image branch:
#		'RAW' > 'TimeSeries' > UID ('SS-001') > CIDfull (CSID_datestamp_timestamp) > [TSi [0,1,2,3,....n], 'acquisition parameters.json', 'configurations.xml'] > [0_0_0_BF_LED_matrix_left_half.bmp, 0_0_0_BF_LED_matrix_right_half.bmp, ..., coordinates.csv]
# Folder structure in processed DPC image branch:
#		'DPC' > 'TimeSeries' > UID ('SS-001') > CIDfull (CSID_datestamp_timestamp) > UID_CSID_TSi_dt_Nt_DPC.png
# Note the DPC image branch will be created in the /scratch space instead of the /project space due to space constraints

def generate_dpc(I1,I2,use_gpu=False):
	if use_gpu:
		# img_dpc = cp.divide(img_left_gpu - img_right_gpu, img_left_gpu + img_right_gpu)
		# to add
		I_dpc = 0
	else:
		I_dpc = np.divide(I1-I2,I1+I2)
		I_dpc = I_dpc + 0.5
	I_dpc[I_dpc<0] = 0
	I_dpc[I_dpc>1] = 1
	return I_dpc

parser = argparse.ArgumentParser(description='Program for processing DPC images')

parser.add_argument("--UID", "--UIDList", help = "Unique ID or UID refers to the unique participant ID specific to each participant, such as 'AA-001', 'SS-003', 'SCD-005', etc.; Input as a list of strings (with spaces): AA-001; or AA-001 AA-002", type = str, nargs='+',default=['Test1'])
args = parser.parse_args()

UID = args.UID
#UID = ['Test1', 'Test2'] # Unique ID or UID refers to the unique participant ID specific to each participant (e.g. 'SCT-004', 'SCD-003', and 'N-007' refer to the fourth, third and seventh participants for classes sickle cell trait, sickle cell disease and normal respectively)

t = time.time() # Record current time to calculate elapsed time

#strLog = '\n' # String with all the information to be saved in log file

for UIDInd in range(len(UID)):	# Run analysis for all UID listed
	strLog = '\n' # String with all the information from current UID to be saved in log file
	t3 = time.time()
	pathParent = os.path.dirname(os.getcwd())  # Path of the folder before current directory, where different UID folders are arranged
	RAWFolder = (pathParent + '/ImagesCanada/' + 'RAW/TimeSeries/' + UID[UIDInd] + '/')	# Path of RAW folder (with all raw images)

	# For the case when processed DPC images are saved in the same path of ImagesCanada
	DPCFolder = (pathParent + '/ImagesCanada/' + 'DPC/TimeSeries/'+ UID[UIDInd] + '/')	# Path of DPC folder (where processed images will be stored)
	# Make directory with tname in DPC
	try:
		os.mkdir(DPCFolder)
	except:
		pass
	# For the case when processed DPC images are saved in a different path, such as the /scratch space in compute Canada
	#DPCFolder = ('/home/shrestha/scratch' + '/ImagesCanada/' + '/DPC/' + UID[UIDInd] )	# Path of DPC folder (where processed images will be stored)

	LogsFolder = ('Logs/CanadaExportDPC/') # The 'Logs' folder is within the current folder (e.g. 'ScriptsCC' folder)

	CSIDfull = os.listdir(RAWFolder)  # Lists the directories or folders in the RAW folder
	CSID = [x.split('_')[0] for x in CSIDfull]  # Create a list containing only the coverslip ID or CSID, such as Ai, Aii, Bi, Bii, and so on

	# Iterate through the folders to read left and right half images and create/save DPC images
	for fInd in range(len(CSIDfull)):
		# Make directory with CSIDfull in DPC folder
		try:
			os.mkdir(DPCFolder + CSIDfull[fInd])
		except:
			pass

		image_format = 'png'
		#image_format = 'tiff'
		t2 = time.time()

		# Read json file
		json_file = (RAWFolder + CSIDfull[fInd] + '/acquisition parameters.json')

		with open(json_file, 'r') as f:
			acquisition_parameters = json.loads(f.read())

		dt = (acquisition_parameters['dt(s)'])
		Nt = (acquisition_parameters['Nt'])

		strCC = 'UID: ' +  UID[UIDInd] + '\t CSID: ' +  CSID[fInd] + '\t Number of image pairs = ' + str(Nt)

		print(strCC) # Print the current configuration
		strLog = strLog + strCC + '\n'	# Save the current configuration to log string

		# Iterate through different image pairs stored in folders 0 to Nt-1
		filesDirCSID = os.listdir(RAWFolder + CSIDfull[fInd] + '/')  # Lists the directories or folders AND files in the CSID folder
		TSimagepairs = [item for item in filesDirCSID if os.path.isdir(os.path.join((RAWFolder + CSIDfull[fInd] + '/'), item))] # Save only the directories (not files - acquisition parameters and configuration files) in TSimagepairs

		if (Nt == len(TSimagepairs)):  # The two are equal if the acquisition completes, and not if it stops before completion
			for TSi in range(Nt):
				i = 0
				j = 0
				k = 0

				file_id = str(i) + '_' + str(j) + '_' + str(k)  # Note: i,j,k here is always 0,0,0
				print(TSi)
				I_BF_left = cv2.imread(RAWFolder + CSIDfull[fInd] + '/'+ str(TSi) +'/' + file_id + '_' + 'BF_LED_matrix_left_half.bmp')
				I_BF_right = cv2.imread(RAWFolder + CSIDfull[fInd] + '/'+ str(TSi) +'/' + file_id + '_' + 'BF_LED_matrix_right_half.bmp')

				# make images float and 0-1
				I_BF_left = I_BF_left.astype('float') / 255
				I_BF_right = I_BF_right.astype('float') / 255
				# generate dpc
				I_DPC = generate_dpc(I_BF_left, I_BF_right)
				# I_DPC_u16 = (I_DPC*65535).astype(np.uint16)
				# File name: UID_CSID_TSi_dt_Nt_DPC.png
				filename = UID[UIDInd] + '_' + '_' + CSID[fInd] + '_' + str(TSi) + '_dt_' + str(dt) + 'seconds_Nt_' + str(Nt) + '_DPC'
				cv2.imwrite(DPCFolder + CSIDfull[fInd] + '/' + filename + '.' + image_format, I_DPC * 255)

				# Calculate and print elapsed time per iteration
				del_t = time.time() - t2
				str_time2 = "{:.3f}".format(del_t)
			strET2 = 'Elapsed Time =' + str_time2 + ' seconds'
			print(strET2)
			strLog = strLog + strET2 + '\n'

		else:
			# Run using i,j,k values from coordinates.csv file if acquisition did not complete
			strN = 'Note: Acquisition not completed for this case, and total number of images (' + str(len(TSimagepairs)) + ') is less than target value (' + str(Nt) + ')'
			print(strN)
			strLog = strLog + strN + '\n'
			for TSi in range(len(TSimagepairs)):
				i = 0
				j = 0
				k = 0

				file_id = str(i) + '_' + str(j) + '_' + str(k)  # Note: i,j,k here is always 0,0,0
				print(TSi)
				I_BF_left = cv2.imread(RAWFolder + CSIDfull[fInd] + '/' + str(TSi) + '/' + file_id + '_' + 'BF_LED_matrix_left_half.bmp')
				I_BF_right = cv2.imread(RAWFolder + CSIDfull[fInd] + '/' + str(TSi) + '/' + file_id + '_' + 'BF_LED_matrix_right_half.bmp')

				# make images float and 0-1
				I_BF_left = I_BF_left.astype('float') / 255
				I_BF_right = I_BF_right.astype('float') / 255
				# generate dpc
				I_DPC = generate_dpc(I_BF_left, I_BF_right)
				# I_DPC_u16 = (I_DPC*65535).astype(np.uint16)
				# File name: UID_CSID_TSi_dt_Nt_DPC.png
				filename = UID[UIDInd] + '_' + '_' + CSID[fInd] + '_' + str(TSi) + '_dt_' + str(dt) + '_Nt_' + str(Nt) + '_DPC'
				cv2.imwrite(DPCFolder + CSIDfull[fInd] + '/' + filename + '.' + image_format, I_DPC * 255)

				# Calculate and print elapsed time per iteration
				del_t = time.time() - t2
				str_time2 = "{:.3f}".format(del_t)
			strET2 = 'Elapsed Time =' + str_time2 + ' seconds'
			print(strET2)
			strLog = strLog + strET2 + '\n'
	# Calculate and print elapsed time for the current UID
	del_t_UID = time.time() - t3
	str_time = "{:.3f}".format(del_t_UID)
	strETUID = 'Total elasped time for ' + UID[UIDInd] + ' =' + str_time + ' seconds'
	print(strETUID)
	strLog = strLog + strETUID + '\n'

	# Write file into /Logs folder
	#logFileName = ' '.join([str(item) for item in UID])
	logFileName = str(UID[UIDInd])
	logFileFull = (LogsFolder + logFileName + '_' + datetime.now().strftime("%d-%m-%Y_%H-%M-%S") + '.txt')
	try:
		with open(logFileFull, 'w') as f:
			f.write(strLog)
	except FileNotFoundError:
		print("The 'Logs' directory does not exist")


# Calculate and print elapsed time
elapsed_time = time.time()-t
str_time = "{:.3f}".format(elapsed_time)
strET = 'Total elapsed time for complete run=' + str_time + ' seconds'
print(strET)
strLog = strLog + strET + '\n'


# Write file into /Logs folder
logFileName = '-'.join([str(item) for item in UID])
logFileFull = (LogsFolder + logFileName +'_' + datetime.now().strftime("%d-%m-%Y_%H-%M-%S") +'.txt')
try:
	with open(logFileFull, 'w') as f:
		f.write(strLog)
except FileNotFoundError:
	print("The 'Logs' directory does not exist")
