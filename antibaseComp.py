from pyteomics import mzxml
import math
import numpy as np
import sys
import csv
from datatoJson import *
import glob
import os
import webbrowser
import time

def importAdductMasses(antiBase_file):
	""" This function is meant to import all of the data from the antiBase file with
	all of the m/z values with all of their adducts """
	with open(antiBase_file) as file:
		reader = csv.reader(file)
		adduct_titles = next(reader)
		adduct_titles = adduct_titles[6:]
		antiBase_dict = {}
		
		for row in reader:
			add = True
			temp_array = []
			for x in range(6, len(row)):
				temp_array.append(row[x])
			
			# Some of the entries in antiBase have no chemical compound name or mass (look at row 17445 for example)
			# This is just meant to throw away those values
			for y in temp_array:
				if not y:
					add = False  
			if add:
				antiBase_dict[row[0]] = [float(z) for z in temp_array]
	return antiBase_dict, adduct_titles

def preprocess_sample(sample, antiBase_dict, adduct_titles, cutoff_percent, scan_level):
	""" Reads mzXML files in and compares their m/z values to those from antiBase file """
	print("Reading " + sample)
	scans = []
	r = mzxml.read(sample)
	while True:
		try:
			scans.append(r.next())
		except:
			break

	print(str(len(scans)) + " scans found in " + sample)
	base_peaks = {}
	all_peaks = []
	filtered_spectra = {}

	# if you choose to compare at the level of MSI, this script will only take peaks which have  
	# intensities that are a significant portion of the TIC (significance is decided by variable cutoff_percent)
	if scan_level == 1: 
		
		for scan in scans:
			if scan['msLevel'] == 1:
				num= int(scan['num'])
				RT= scan['retentionTime']
				intensity_array_MSI = scan['intensity array']
				mzs_MSI= scan['m/z array']
				#print(mzs_MSI)
				total_TIC = scan['totIonCurrent']
				
				filtered_spectra_properties_list = []
				for i in range(0, len(mzs_MSI)):
					if intensity_array_MSI[i]/total_TIC >= cutoff_percent:
						filtered_spectra_properties_list.append(float(mzs_MSI[i]))
						filtered_spectra_properties_list.append(RT)
						filtered_spectra[num] = filtered_spectra_properties_list
	# if you compare at the level of MS-II, this only takes peaks which generate MS2's 
	if scan_level == 2: 
		for scan in scans:
			if scan['msLevel'] == 2:
				num= int(scan['num'])
				RT= scan['retentionTime']
				parent_mz = scan['precursorMz'][0]['precursorMz']
				filtered_spectra_properties_list = []
				filtered_spectra_properties_list.append(parent_mz)
				filtered_spectra_properties_list.append(RT)
				filtered_spectra[num] = filtered_spectra_properties_list
	scans = None
	
	print("Scans have been scanned")
	return antiBase_dict, filtered_spectra

def antiBase_search(antiBase_dict, filtered_spectra, from_csv = True ): 
	print("now going through antiBase")
	output_mapping_dict = {}
	# compares spectra in mzXML to antibase m/z's. Spectra are considered similar if the error is less than 0.1 ppm
	for antiBase_key in antiBase_dict.keys():
		antiBase_adduct_map = []
		scan_name_map = []
		scan_RT_map = []
		scan_mz_map = []
		scan_ppm_map = []
		counter = 0
		for spectra_key in filtered_spectra.keys():
			if from_csv:
				spec_num = filtered_spectra[spectra_key]
				spec_RT = None
			else:
				spec_num = filtered_spectra[spectra_key][0]
				spec_RT = filtered_spectra[spectra_key][1]
			for j in range(0, len(antiBase_dict[antiBase_key])):
		
				base_num = antiBase_dict[antiBase_key][j]
				diff = abs(base_num - spec_num)/spec_num

				if diff < 0.0000001:
					counter += 1
					antiBase_adduct_map.append(adduct_titles[j])
					scan_name_map.append(spectra_key)
					scan_RT_map.append(spec_RT)
					scan_mz_map.append(spec_num)
					scan_ppm_map.append(diff * 10**8)
		#  only allows matches that have a hydrogen adduct or sodium adduct to continue
		if "Ion_Mass_M+H" in antiBase_adduct_map or "Ion_Mass_M+Na" in antiBase_adduct_map:
			output_mapping_dict[antiBase_key] = [antiBase_adduct_map,scan_name_map, scan_RT_map,scan_mz_map,scan_ppm_map]

	return output_mapping_dict		

def filter_sample(output_mapping_dict):
	""" This function filters out similar hits that have the same RT value, as it is probable that things 
	with the same RT value are the same thing """
	print("Now filtering")
	filtered_output_mapping_dict = {}
	for key in output_mapping_dict.keys():
		uniqueRT = 0  
		RT_list = output_mapping_dict[key][2]
		scan_mz_map = output_mapping_dict[key][3]
		pop_list = []
		for i in range(0, len(RT_list)):
			if abs(RT_list[i] - uniqueRT)/RT_list[i] <= 0.1:
				pop_list.append(i)
			else:
				uniqueRT = RT_list[i]

		filtered_antiBase_adduct_map = [output_mapping_dict[key][0][a] for a in range(0,len(scan_mz_map)) if a not in pop_list]
		filtered_scan_name_map = [output_mapping_dict[key][1][b] for b in range(0,len(scan_mz_map)) if not b in pop_list]
		filtered_scan_RT_map = [output_mapping_dict[key][2][c] for c in range(0,len(scan_mz_map)) if c not in pop_list]
		filtered_scan_mz_map = [output_mapping_dict[key][3][d] for d in range(0,len(scan_mz_map)) if d not in pop_list]
		filtered_scan_ppm_map = [output_mapping_dict[key][4][e] for e in range(0,len(scan_mz_map)) if e not in pop_list]
		filtered_output_mapping_dict[key] = [filtered_antiBase_adduct_map, filtered_scan_name_map, filtered_scan_RT_map, filtered_scan_mz_map, filtered_scan_ppm_map]

	filtered_output_mapping_dict= {keys:values for keys,values in filtered_output_mapping_dict.items() if filtered_output_mapping_dict[keys][0] and len(filtered_output_mapping_dict[keys][0])>=1}
	output_mapping_dict = None
	return filtered_output_mapping_dict

start_time = time.time()
antiBase_dict,adduct_titles = importAdductMasses("new_antiBase_file.csv")
for fname in os.listdir('.'):
	if fname.endswith('.mzXML'):
		antiBase_dict, filtered_spectra = preprocess_sample(fname, antiBase_dict,adduct_titles, 0.1,2)
		k = antiBase_search(antiBase_dict,filtered_spectra)
		k = filter_sample(k)
		makeJson(k)
	elif fname.endswith('.csv') and fname.startswith("avg"):
		filtered_spectra = {}

		with open(fname) as mz_file: 
			reader = csv.reader(mz_file)
			next(reader)
			for row in reader:
				#print(row)
				#input("...")
				filtered_spectra[row[0]] = float(row[1])
		k = antiBase_search(antiBase_dict,filtered_spectra)

		makeJson(k)
print("--- %s seconds ---" % (time.time() - start_time))
os.system("python3 -m http.server")


