# Script to process DPC images based on Left and Right half LED illuminated images implemented in Compute Canada
# Based on script from H. Li (2022)
# 2023 ps ubc

import time
import cv2
import os
import numpy as np
import json
from datetime import datetime
import pandas as pd

# Script saved in folder ('ScriptsCC' - can be named anything), under the root directory with another directory "ImagesNepal" with all the images
# Folder structure in RAW image branch:
#		UID ('SS-001') > 'RAW' > tname ('t02hrs') > CIDfull (CSID_datestamp_timestamp) > [0, 'acquisition parameters.json', 'configurations.xml'] > [i_j_k_BF_LED_matrix_left_half.bmp, i_j_k_BF_LED_matrix_right_half.bmp, ..., coordinates.csv]
# Folder structure in processed DPC image branch:
#		UID ('SS-001') > 'DPC' > tname ('t02hrs') > CIDfull (CSID_datestamp_timestamp) > UID_tname_CSID_i_j_k_DPC.png
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


UID = ['AA-001', 'AA-002'] # Unique ID or UID refers to the unique participant ID specific to each participant (e.g. 'SCT-004', 'SCD-003', and 'N-007' refer to the fourth, third and seventh participants for classes sickle cell trait, sickle cell disease and normal respectively)

t = time.time() # Record current time to calculate elapsed time

#strLog = '\n' # String with all the information to be saved in log file

for UIDInd in range(len(UID)):	# Run analysis for all UID listed
	strLog = '\n' # String with all the information from current UID to be saved in log file
	t3 = time.time()
	pathParent = os.path.dirname(os.getcwd())  # Path of the folder before current directory, where different UID folders are arranged
	RAWFolder = (pathParent + '/ImagesNepal/' + UID[UIDInd] + '/RAW/')	# Path of RAW folder (with all raw images)

	# For the case when processed DPC images are saved in the same path of ImagesNepal
	DPCFolder = (pathParent + '/ImagesNepal/' + UID[UIDInd] + '/DPC/')	# Path of DPC folder (where processed images will be stored)
	# For the case when processed DPC images are saved in a different path, such as the /scratch space in compute Canada
	#DPCFolder = ('/home/shrestha/scratch' + '/ImagesNepal/' + UID[UIDInd] + '/DPC/')	# Path of DPC folder (where processed images will be stored)

	LogsFolder = ('Logs/NepalExportDPC/') # The 'Logs' folder is within the current folder (e.g. 'ScriptsCC' folder)

	dirInRAW = os.listdir(RAWFolder)  # Lists the directories or folders in the RAW folder

	tname = [indx for indx in dirInRAW if indx[0].lower() == 't'] # Save all directory names starting with t in the RAW folder (Using lower(), just in case files were saved using 'T' instead of 't')

	for tnameInd in range(len(tname)): # Loop through all folders starting with t, such as t0, t2hrs, and so on
		CSIDfull = os.listdir(RAWFolder + tname[tnameInd] + '/')  # Lists the directories or folders in the tname folder
		CSID = [x.split('_')[0] for x in CSIDfull]  # Create a list containing only the coverslip ID or CSID, such as Ai, Aii, Bi, Bii, and so on

		# Make directory with tname in DPC
		try:
			os.mkdir(DPCFolder + tname[tnameInd])
		except:
			pass

		# Iterate through the folders to read left and right half images and create/save DPC images
		for fInd in range(len(CSIDfull)):
			# Make directory with CSIDfull in current tname folder
			try:
				os.mkdir(DPCFolder + tname[tnameInd] + '/' + CSIDfull[fInd])
			except:
				pass

			image_format = 'png'
			#image_format = 'tiff'
			t2 = time.time()

			#Read csv file
			csv_file = (RAWFolder + tname[tnameInd] + '/' + CSIDfull[fInd] + '/0/coordinates.csv')
			df = pd.read_csv(csv_file)
			iList = df.i.astype(int)
			jList = df.j.astype(int)
			kList = df.k.astype(int)

			# Read json file
			json_file = (RAWFolder + tname[tnameInd] + '/' + CSIDfull[fInd] + '/acquisition parameters.json')

			with open(json_file, 'r') as f:
				acquisition_parameters = json.loads(f.read())

			Nx = (acquisition_parameters['Nx'])
			Ny = (acquisition_parameters['Ny'])

			strCC = 'UID: ' +  UID[UIDInd] + '\t tname: ' + tname[tnameInd] + '\t CSID: ' +  CSID[fInd] + '\t Number of image pairs = ' + str(len(iList))

			print(strCC) # Print the current configuration
			strLog = strLog + strCC + '\n'	# Save the current configuration to log string

			if ((Nx*Ny) == len(iList)):  # The two are equal if the acquisition completes, and not if it stops before completion

				# Run using Nx and Ny extracted from json file (if acquisition completed)
				for i in range(Ny):		# Note: i refers to Y and Ny
					for j in range(Nx): # Note: j refers to X and Nx
						k = 0  # No z scanning
						file_id = str(i) + '_' + str(j) + '_' + str(k) # Note: i,j,k corresponds to y,x,z and not x,y,z
						print(file_id)
						I_BF_left = cv2.imread(RAWFolder + tname[tnameInd] + '/' + CSIDfull[fInd] + '/0/' + file_id + '_' + 'BF_LED_matrix_left_half.bmp')
						I_BF_right = cv2.imread(RAWFolder + tname[tnameInd] + '/' + CSIDfull[fInd] + '/0/' + file_id + '_' + 'BF_LED_matrix_right_half.bmp')

						# make images float and 0-1
						I_BF_left = I_BF_left.astype('float') / 255
						I_BF_right = I_BF_right.astype('float') / 255
						# generate dpc
						I_DPC = generate_dpc(I_BF_left, I_BF_right)
						#I_DPC_u16 = (I_DPC*65535).astype(np.uint16)
						filename = UID[UIDInd] + '_' + tname[tnameInd] + '_' + CSID[fInd] + '_' + file_id + '_DPC'
						cv2.imwrite(DPCFolder + tname[tnameInd] + '/' + CSIDfull[fInd] + '/' + filename + '.' + image_format, I_DPC * 255)
						#cv2.imwrite(pathParent + '/ImagesVanc/' + UID[UIDInd] + '/DPC/' + tname[tnameInd] + '/' + CSIDfull[fInd] + '/' + filename + '.' + image_format, I_DPC * 255,[cv2.IMWRITE_PNG_COMPRESSION, 0])
						#cv2.imwrite(pathParent + '/ImagesVanc/' + UID[UIDInd] + '/DPC/' + tname[tnameInd] + '/' + CSIDfull[fInd] + '/' + filename + '.' + image_format, I_DPC_u16)
						#print(I_DPC_u16.dtype)

						# Calculate and print elapsed time per iteration
						del_t = time.time() - t2
						str_time2 = "{:.3f}".format(del_t)
				strET2 = 'Elapsed Time =' + str_time2 + ' seconds'
				print(strET2)
				strLog = strLog + strET2 + '\n'

			else:
				# Run using i,j,k values from coordinates.csv file if acquisition did not complete
				strN = 'Note: Acquisition not completed for this case, and total number of images (' + str(len(iList)) + ') is less than target value (' + str(Nx*Ny) + ')'
				print(strN)
				strLog = strLog + strN + '\n'
				for ii in range(len(iList)):

					if (iList[ii] % 2 == 0): # Account for mismatch of j between coordinates.csv and saved images
						file_id = str(iList[ii]) + '_' + str(jList[ii]) + '_' + str(kList[ii])	# For even j, no mismatch
					else:
						file_id = str(iList[ii]) + '_' + str(Nx - jList[ii] - 1) + '_' + str(kList[ii]) # For odd j, images save in descending order of j, while .csv files it saves in ascending order
					print(file_id)
					I_BF_left = cv2.imread(RAWFolder + tname[tnameInd] + '/' + CSIDfull[fInd] + '/0/' + file_id + '_' + 'BF_LED_matrix_left_half.bmp')
					I_BF_right = cv2.imread(RAWFolder + tname[tnameInd] + '/' + CSIDfull[fInd] + '/0/' + file_id + '_' + 'BF_LED_matrix_right_half.bmp')

					# make images float and 0-1
					I_BF_left = I_BF_left.astype('float') / 255
					I_BF_right = I_BF_right.astype('float') / 255
					# generate dpc
					I_DPC = generate_dpc(I_BF_left, I_BF_right)
					#I_DPC_u16 = (I_DPC * 65535).astype(np.uint16)

					filename = UID[UIDInd] + '_' + tname[tnameInd] + '_' + CSID[fInd] + '_' + file_id + '_DPC'
					cv2.imwrite(DPCFolder + tname[tnameInd] + '/' + CSIDfull[fInd] + '/' + filename + '.' + image_format, I_DPC * 255)
					#cv2.imwrite(pathParent + '/ImagesVanc/' + UID[UIDInd] + '/DPC/' + tname[tnameInd] + '/' + CSIDfull[fInd] + '/' + filename + '.' + image_format, I_DPC * 255, [cv2.IMWRITE_PNG_COMPRESSION, 0])
					#cv2.imwrite(pathParent + '/ImagesVanc/' + UID[UIDInd] + '/DPC/' + tname[tnameInd] + '/' + CSIDfull[fInd] + '/' + filename + '.' + image_format, I_DPC_u16) # 65535 for 16 bit, 255 for 8 bit
					#print(I_DPC_u16.dtype)

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
	# logFileName = ' '.join([str(item) for item in UID])
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
