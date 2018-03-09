from pyteomics import mzxml
import math
import numpy as np
import sys
import csv
from datatoJson import *

def importAdductMasses(antiBase_file):
	with open(antiBase_file) as file:
		reader = csv.reader(file)
		adduct_titles = next(reader)
		adduct_titles = adduct_titles[6:]
		print(adduct_titles)
		antiBase_dict = {}
		
		for row in reader:
			add = True
			temp_array = []
			for x in range(6, len(row)):
				temp_array.append(row[x])
			for y in temp_array:
				if not y:
					add = False  
			if add:
				antiBase_dict[row[0]] = [float(z) for z in temp_array]
	return antiBase_dict, adduct_titles

def preprocess_sample(sample, antiBase_dict, adduct_titles, cutoff_percent, scan_level):
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
	output_mapping_dict = {}
	print("Scans have been scanned")
	print("now going through antiBase")

	for antiBase_key in antiBase_dict.keys():
		antiBase_adduct_map = []
		scan_name_map = []
		scan_RT_map = []
		scan_mz_map = []
		counter = 0
		for spectra_key in filtered_spectra.keys():

			for j in range(0, len(antiBase_dict[antiBase_key])):
				base_num = antiBase_dict[antiBase_key][j]
				spec_num = filtered_spectra[spectra_key][0]
				spec_RT = filtered_spectra[spectra_key][1]

				diff = abs(base_num - spec_num)/spec_num
				if diff < 0.0000001:
					counter += 1
					antiBase_adduct_map.append(adduct_titles[j])
					scan_name_map.append(spectra_key)
					scan_RT_map.append(spec_RT)
					scan_mz_map.append(spec_num)
					#print(adduct_titles[j])
					#input("...")
					#print("Ion_Mass_M+Na" in antiBase_adduct_map)
					#input("...")
		if "Ion_Mass_M+H" in antiBase_adduct_map or "Ion_Mass_M+Na" in antiBase_adduct_map:
			output_mapping_dict[antiBase_key] = [antiBase_adduct_map,scan_name_map, scan_RT_map,scan_mz_map]

	return output_mapping_dict		

def filter_sample(output_mapping_dict):
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
				#pop_list.append(i)
				uniqueRT = RT_list[i]

		filtered_antiBase_adduct_map = [output_mapping_dict[key][0][a] for a in range(0,len(scan_mz_map)) if a not in pop_list]
		filtered_scan_name_map = [output_mapping_dict[key][1][b] for b in range(0,len(scan_mz_map)) if not b in pop_list]
		filtered_scan_RT_map = [output_mapping_dict[key][2][c] for c in range(0,len(scan_mz_map)) if c not in pop_list]
		filtered_scan_mz_map = [output_mapping_dict[key][3][d] for d in range(0,len(scan_mz_map)) if d not in pop_list]
		filtered_output_mapping_dict[key] = [filtered_antiBase_adduct_map, filtered_scan_name_map, filtered_scan_RT_map, filtered_scan_mz_map]
		#print(filtered_output_mapping_dict[key])

	filtered_output_mapping_dict= {keys:values for keys,values in filtered_output_mapping_dict.items() if filtered_output_mapping_dict[keys][0] and len(filtered_output_mapping_dict[keys][0])>=1}
	output_mapping_dict = None
	return filtered_output_mapping_dict
import time
start_time = time.time()
antiBase_dict,adduct_titles = importAdductMasses("new_antiBase_file.csv")
k = preprocess_sample("010818_BH14_pos.mzXML", antiBase_dict,adduct_titles, 0.1,2)
k = filter_sample(k)
makeJson(k)
print("--- %s seconds ---" % (time.time() - start_time))
