from pyteomics import mzxml
import math
import numpy as np
import sys
import csv

def importAdductMasses(antiBase_file):
	with open(antiBase_file) as file:
		reader = csv.reader(file)
		adduct_titles = reader.next()
		antiBase_dict = {}
		for row in reader:
			temp_array = [row[x] for x in range(6, len(row))]
			antiBase_dict[row[0]] = temp_array
	return antiBase_dict, adduct_title

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
					filtered_spectra_properties_list.append(mzs_MSI[i])
					filtered_spectra_properties_list.append(intensity_array_MSI[i])
					filtered_spectra[num] = filtered_spectra_properties_list
	scans = None
	output_mapping_dict = {}
	for spectra_key in filtered_spectra.keys():
		for antiBase_key in antiBase_dict.keys():
			counter = 0
			antiBase_map = []
			for j in range(0, len(antiBase_dict[antiBase_key])):
				if math.abs(antiBase_dict[antiBase_key][j] - filtered_spectra[spectra_key][0])/max(antiBase_dict[antiBase_key][j], filtered_spectra[spectra_key][0])< 0.0001: 
					counter += 1
					antiBase_map.append(adduct_titles[j])
			if counter >=2:
				output_mapping_dict[spectra_key] = antiBase_map
				break

	return output_mapping_dict
antiBase_dict = importAdductMasses("new_antiBase_file.csv")
