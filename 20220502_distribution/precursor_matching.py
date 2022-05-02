# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 10:57:47 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

spectra_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\SG-CHH-Fraction-UnoptimizedMethod_NewLCgradient_1_ms2_output_list.txt" #path to directory containing spectra
fragment_list_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\by_ions.csv" #path to directory with list of target fragment ions
target_list_import = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\precursor_list.csv"#path to directory with inclusion lists
final_dir =r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\output_directory" #path to directory for final, processed data

sample_name = 'SG-CHH-Fraction-UnoptimizedMethod_NewLCgradient_1'

error_precursor = 10 #+/- ppm, for precursor
error_fragment = 0.02 #+/- Da, for fragment ion, charge state 1

precursor_charges = [7]
fragment_charges = [1,2,3,4,5]

h_mass = 1.00784

print('loading files', datetime.now())

spectra_read= pd.read_csv(spectra_import, sep=" ",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','precursor_charge','empty'])

spectra_value = pd.DataFrame()
spectra_value['Fragment m/z'] = spectra_read['m/z']
spectra_value['resolution'] = spectra_read['resolution']
spectra_value['charge'] = spectra_read['charge']
spectra_value['intensity'] = spectra_read['intensity']
spectra_value['Precursor'] = spectra_read['MS2']
spectra_value['Scan #'] = spectra_read['scan_number']
spectra_value['Precursor_Charge'] = spectra_read['precursor_charge']

#spectra_value = spectra_value(spectra_value[spectra_value['Scan #']>40000].index) #remove first 40000 scans
spectra_value = spectra_value.drop(spectra_value[spectra_value['Precursor_Charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['intensity']<100].index) #remove intensity values less than 100

target_list = pd.read_csv(target_list_import)
MS1_sequence_storage = []
MS1_file_type_storage = []
MS1_scan_storage = []
MS1_theo_precursor_store = []
MS1_actual_precursor_store = []
MS1_precursor_theo_charge = []
MS1_precursor_actual_charge = []
MS2_fragment_charge = []
MS2_fragment_mz = []
MS1_intensity = []
MS1_resolution = []
error_storage = []

for b in precursor_charges:
        precursor_target = target_list[str(b)].values.tolist()
        spectra_value_filtered = spectra_value[spectra_value['Precursor_Charge'] == b]
        for c in precursor_target:
            spectra_value_filtered['ppm error'] =  ((abs(spectra_value_filtered['Precursor'] - c))/c) * 1E6
            spectra_value_filtered_err = spectra_value_filtered[spectra_value_filtered['ppm error'] <= error_precursor]
            spectra_value_filtered_err['Theoretical precursor'] = c
            spectra_value_filtered_err['Theoretical precursor charge'] = b
            if len(spectra_value_filtered_err) > 0:
                scan = spectra_value_filtered_err['Scan #'].values.tolist()
                theo_precursor = spectra_value_filtered_err['Theoretical precursor'].values.tolist()
                actual_precursor = spectra_value_filtered_err['Precursor'].values.tolist()
                precursor_error = spectra_value_filtered_err['ppm error'].values.tolist()
                precursor_theretical_charge = spectra_value_filtered_err['Theoretical precursor charge'].values.tolist()
                actual_precursor_charge = spectra_value_filtered_err['Precursor_Charge'].values.tolist()
                fragment_charge = spectra_value_filtered_err['charge'].values.tolist()
                fragment_mz = spectra_value_filtered_err['Fragment m/z'].values.tolist()
                intensity = spectra_value_filtered_err['intensity'].values.tolist()
                resolution = spectra_value_filtered_err['resolution'].values.tolist()
                
                for e in scan:
                    MS1_scan_storage.append(e)
                for f in theo_precursor:
                    MS1_theo_precursor_store.append(f)
                for g in actual_precursor:
                    MS1_actual_precursor_store.append(g)
                for h in precursor_error:
                    error_storage.append(h)
                for i in precursor_theretical_charge:
                    MS1_precursor_theo_charge.append(i)
                for j in actual_precursor_charge:
                    MS1_precursor_actual_charge.append(j)
                for k in fragment_charge:
                    MS2_fragment_charge.append(k)
                for l in fragment_mz:
                    MS2_fragment_mz.append(l)
                for m in intensity:
                    MS1_intensity.append(m)
                for n in resolution:
                    MS1_resolution.append(n)

precursor_matches = pd.DataFrame()

precursor_matches['Scan #'] = MS1_scan_storage
precursor_matches['Theoretical precursor'] = MS1_theo_precursor_store
precursor_matches['Actual precursor'] = MS1_actual_precursor_store
precursor_matches['ppm error'] = error_storage
precursor_matches['Precursor theoretical charge'] = MS1_precursor_theo_charge
precursor_matches['Precursor actual charge'] = MS1_precursor_actual_charge
precursor_matches['Fragment m/z'] = MS2_fragment_mz
precursor_matches['Fragment ion charge'] = MS2_fragment_charge
precursor_matches['Fragment intensity'] = MS1_intensity
precursor_matches['Fragment resolution'] = MS1_resolution

file_name = sample_name + 'precursor_matches.txt'
file_out_path = final_dir + '\\' + file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        precursor_matches.to_csv(filec,index=False)