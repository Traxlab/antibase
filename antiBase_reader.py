import csv 
import sys
import pandas as pd

def processAntiBase(antiBase_file, new_antiBase_file):
	""" This function will take in the Antibase file and output a new file in which the value separated by commas 
		will be separated by columns instead. These comma separated values include molecular weight, 
		name of compound, and something else"""
	row_info = []
	with open(antiBase_file) as file: 
		reader = csv.reader(file) 
		next(reader)
		for row in reader: 
			row_info.append(row[0])

	csv_vals = []
	for comma_values in row_info:
		temp_array = []
		temp_array = comma_values.split(",")
		csv_vals.append(temp_array)

	with open(new_antiBase_file, "w") as new_file:
		writer = csv.writer(new_file)
		writer.writerow(["Name", "RT", "Molecular_Weight"])
		writer.writerows(csv_vals)
	f1 = pd.read_csv(antiBase_file)
	f1 = f1.drop(["name"],axis =1)
	f2 = pd.read_csv(new_antiBase_file)
	f2 = pd.concat([f2,f1], axis =1)
	f2.to_csv(new_antiBase_file, sep = ",")


def adductCalculator(new_antiBase_file, adduct_mass,adduct_divide_factor, adduct_Name):
	### This will add adducts to the molecular weight of all compounds in antiBase"""
	df = pd.read_csv(new_antiBase_file)


	adduct_series = (df.Mass * adduct_divide_factor)+ adduct_mass
	adduct_df = adduct_series.to_frame()
	#print(adduct_df)
	#adduct_df.rename(index = str, { "Molecular_Weight" : "Molecular_Weight"+adduct_Name})
	adduct_df.columns = ["Ion_Mass_" +adduct_Name]
	#print(adduct_df)

	df = pd.concat([df,adduct_df],axis = 1)
	cols=pd.Series(df.columns)
	#for dup in df.columns.get_duplicates(): cols[df.columns.get_loc(dup)]=[dup+'_'+adduct_Name if d_idx!=0 else dup for d_idx in range(df.columns.get_loc(dup).sum())]


	df.columns=cols
	df.to_csv(new_antiBase_file, sep = ",", index = False)



#processAntiBase("antibase_no_nanicuse_JPG_agilent.csv","new_antiBase_file.csv")
adduct_masses = [1.007276,8.33459,15.76619,22.989218,1.007276,9.52055,11.998247,19.985217,21.52055,22.989218,42.033823,62.547097,1.007276,18.033823,22.989218,33.033489,38.963158,42.033823,44.97116,61.06534,64.015765,76.91904,79.02122,83.06037,84.05511,1.007276,18.033823,22.989218,38.963158,42.033823,64.015765]
adduct_divide_factor = [0.33,0.33,0.33,0.33,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2]
adduct_names = ["M+3H","M+2H+Na","M+H+2Na","M+3Na","M+2H","M+H+NH4","M+H+Na","M+H+K","M+ACN+2H","M+2Na","M+2ACN+2H","M+3ACN+2H","M+H","M+NH4","M+Na","M+CH3OH+H","M+K","M+ACN+H","M+2Na-H","M+IsoProp+H","M+ACN+Na","M+2K-H","M+DMSO+H","M+2ACN+H","M+IsoProp+Na+H","2M+H","2M+NH4","2M+Na","2M+K","2M+ACN+H","2M+ACN+Na"]
for i in range(0,len(adduct_masses)):
	adductCalculator("new_antiBase_file.csv", adduct_masses[i], adduct_divide_factor[i], adduct_names[i])
