from pyteomics import mzxml
import math
import numpy as np
import sys
import csv

def importAdductMasses(antiBase_file):
	with open(antiBase_file) as file:
		reader = csv.reader(file)
		adduct_titles = next(reader)
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
			

				#temp_array = [float(row[x]) if row[x] != None for x in range(6, len(row)) ]
			#print(temp_array)
			#input("...")
			#for j in temp_array:
			#	print(float(j.replace(" ", "")))
			#print(temp_array)
			if add:
				antiBase_dict[row[0]] = [float(z) for z in temp_array]
				#input("...")
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
			intensity_array_MSI= scan['intensity array']
			mzs_MSI= scan['m/z array']
			total_TIC = scan['totIonCurrent']
			
			filtered_spectra_properties_list = []
			for i in range(0, len(intensity_array_MSI)):
				if intensity_array_MSI[i]/total_TIC >= cutoff_percent:
					filtered_spectra_properties_list.append(float(mzs_MSI[i]))
					filtered_spectra_properties_list.append(intensity_array_MSI[i])
					filtered_spectra[num] = filtered_spectra_properties_list
	scans = None
	output_mapping_dict = {}
	print("Scans have been scanned")
	print("now going through antiBase")

	for spectra_key in filtered_spectra.keys():
		for antiBase_key in antiBase_dict.keys():
			#print(filtered_spectra[spectra_key])
			#print(antiBase_dict[antiBase_key])
			#input("...")
			counter = 0
			antiBase_map = []
			for j in range(0, len(antiBase_dict[antiBase_key])):
				base_num = antiBase_dict[antiBase_key][j]
				spec_num = filtered_spectra[spectra_key][0]
				#print(base_num, spec_num)
				#print(type(base_num), type(spec_num))
				diff = abs(base_num - spec_num)/spec_num
				if diff < 0.00001:
					#print(diff)
					#print("Hit") 
					counter += 1
					antiBase_map.append(adduct_titles[j])
					#input("....")
				if counter >=2:
					print("Multiple adducts found")
					output_mapping_dict[spectra_key] = antiBase_map
					antiBase_map = None
					break
			#if counter >=2:
			#	break

	#return output_mapping_dict
import time
start_time = time.time()
antiBase_dict,adduct_titles = importAdductMasses("new_antiBase_file.csv")
preprocess_sample("010818_BH14_pos.mzXML", antiBase_dict,adduct_titles, 0.1)
print("--- %s seconds ---" % (time.time() - start_time))
