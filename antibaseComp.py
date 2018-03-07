from pyteomics import mzxml
import math
import numpy as np
import sys
import csv

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

def preprocess_sample(sample, antiBase_dict, adduct_titles, cutoff_percent):
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
	for scan in scans:
		if scan['msLevel'] == 1:
			num= int(scan['num'])
			RT= scan['retentionTime']#float(scan['retentionTime'].replace("P","").replace("T","").replace("S",""))
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
	scans = None
	output_mapping_dict = {}
	print("Scans have been scanned")
	print("now going through antiBase")

	for antiBase_key in antiBase_dict.keys():
		antiBase_adduct_map = []
		scan_name_map = []
		scan_RT_map = [] #initialize array with impossible value
		counter = 0
		for spectra_key in filtered_spectra.keys():

			for j in range(0, len(antiBase_dict[antiBase_key])):
				base_num = antiBase_dict[antiBase_key][j]
				spec_num = filtered_spectra[spectra_key][0]
				spec_RT = filtered_spectra[spectra_key][1]
				#print(base_num, spec_num)
				#print(type(base_num), type(spec_num))
				diff = abs(base_num - spec_num)/spec_num
				if diff < 0.0001:
					# if not scan_RT_map:
						#print("new hit")
					counter += 1
					antiBase_adduct_map.append(adduct_titles[j])
					scan_name_map.append(spectra_key)
					scan_RT_map.append(spec_RT)
					# else:
					# 	#for ret in scan_RT_map:
					# 	#	if abs(ret - spec_RT)/spec_RT <= 0.001:
					# 			#print("new hit")
					# 			counter += 1
					# 			antiBase_adduct_map.append(adduct_titles[j])
					# 			scan_name_map.append(spectra_key)
					# 			scan_RT_map.append(spec_RT)
					#input("....")
		output_mapping_dict[antiBase_key] = [antiBase_adduct_map,scan_name_map, spec_RT]
	for key in output_mapping_dict.keys():
		uniqueRT = 0  
		RT_list = output_mapping_dict[key][2]
		for i in range(0, len(RT_list)):
			if 
				
			#if counter >=2:
			#	break

	return output_mapping_dict
import time
start_time = time.time()
antiBase_dict,adduct_titles = importAdductMasses("new_antiBase_file.csv")
k= preprocess_sample("010818_BH14_pos.mzXML", antiBase_dict,adduct_titles, 0.1)
print("--- %s seconds ---" % (time.time() - start_time))
