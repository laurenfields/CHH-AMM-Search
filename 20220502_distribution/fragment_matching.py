# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 13:34:48 2022

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
precursor_matches = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\output_directory\SG-CHH-Fraction-UnoptimizedMethod_NewLCgradient_1precursor_matches.txt"#path to directory with inclusion lists
final_dir =r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\output_directory" #path to directory for final, processed data

sample_name = 'SG-CHH-Fraction-UnoptimizedMethod_NewLCgradient_1'

error_precursor = 10 #+/- ppm, for precursor
error_fragment = 0.02 #+/- Da, for fragment ion, charge state 1

precursor_charges = [7]
fragment_charges = [1,2,3,4,5]

h_mass = 1.00784

precursor_matches = pd.read_csv(precursor_matches, sep=",",skiprows=[0], names= ["Scan #", "Theoretical precursor", "Actual precursor", "ppm error","Precursor theoretical charge",
                                                                                 'Precursor actual charge','Fragment m/z','Fragment ion charge','Fragment intensity','Fragment resolution'])

fragment_database = pd.read_csv(fragment_list_import)

spectra_read= pd.read_csv(spectra_import, sep=" ",skiprows=[0], names= ["m/z", "resolution", "charge", "intensity","MS2",'scan_number','precursor_charge','empty'])
spectra_value = pd.DataFrame()
spectra_value['Fragment m/z'] = spectra_read['m/z']
spectra_value['resolution'] = spectra_read['resolution']
spectra_value['charge'] = spectra_read['charge']
spectra_value['intensity'] = spectra_read['intensity']
spectra_value['Precursor'] = spectra_read['MS2']
spectra_value['Scan #'] = spectra_read['scan_number']
spectra_value['Precursor_Charge'] = spectra_read['precursor_charge']
#spectra_value = spectra_value.drop(spectra_value[spectra_value['Scan #']<40000].index) #remove first 40000 scans
spectra_value = spectra_value.drop(spectra_value[spectra_value['Precursor_Charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['charge']==0].index) #remove charges equal to 0
spectra_value = spectra_value.drop(spectra_value[spectra_value['intensity']<100].index) #remove intensity values less than 100

frag_possible_seqeuence_store = []
frag_scan_store = []
frag_theo_precursor_store = []
frag_act_precursor_store = []
frag_precursor_err_store = []
frag_precursor_act_charge_store = []
frag_precursor_theo_charge_store = []
frag_mz_store = []
frag_ion_charge_store = []
frag_intensity_store = []
frag_resoltuion_store = []
frag_act_M_store = []
frag_theo_mz_store = []
frag_theo_M_store = []
Da_err_store = []

for a in fragment_charges: 
    theo_fragment = fragment_database[str(a)].values.tolist() #theoretical fragments for each charge
    seq_exp = precursor_matches[precursor_matches['Fragment ion charge'] == a] #actual fragments for each charge
    for u in theo_fragment:
        theo_frag_M = (u * a) - (h_mass * a)
        seq_exp['fragment M'] = (seq_exp['Fragment m/z'] * a) - (h_mass * a)
        seq_exp['theoretical fragment'] = u
        seq_exp['theoretical fragment M'] = theo_frag_M
        seq_exp['Da error'] = abs(seq_exp['fragment M'] - theo_frag_M)
        seq_exp_filter = seq_exp.sort_values(by='Da error')
        seq_exp_filter = seq_exp_filter[seq_exp_filter['Da error'] <= error_fragment]           
        if len(seq_exp_filter) > 0:
                        scan = seq_exp_filter['Scan #'].values.tolist()
                        theo_precursor = seq_exp_filter['Theoretical precursor'].values.tolist()
                        actual_precursor = seq_exp_filter['Actual precursor'].values.tolist()
                        precursor_error = seq_exp_filter['ppm error'].values.tolist()
                        precursor_theretical_charge = seq_exp_filter['Precursor theoretical charge'].values.tolist()
                        actual_precursor_charge = seq_exp_filter['Precursor actual charge'].values.tolist()
                        fragment_charge = seq_exp_filter['Fragment ion charge'].values.tolist()
                        fragment_mz = seq_exp_filter['Fragment m/z'].values.tolist()
                        intensity = seq_exp_filter['Fragment intensity'].values.tolist()
                        resolution = seq_exp_filter['Fragment resolution'].values.tolist()
                        frag_act_M = seq_exp_filter['fragment M'].values.tolist()
                        theo_frag = seq_exp_filter['theoretical fragment'].values.tolist()
                        theo_frag_M = seq_exp_filter['theoretical fragment M'].values.tolist()
                        Da_err = seq_exp_filter['Da error'].values.tolist()

                        for w in scan:
                            frag_scan_store.append(w)
                        for x in theo_precursor:
                            frag_theo_precursor_store.append(x)
                        for y in actual_precursor:
                            frag_act_precursor_store.append(y)
                        for z in precursor_error:
                            frag_precursor_err_store.append(z)
                        for aa in precursor_theretical_charge:
                            frag_precursor_theo_charge_store.append(aa)
                        for ab in actual_precursor_charge:
                            frag_precursor_act_charge_store.append(ab)
                        for ac in fragment_charge:
                            frag_ion_charge_store.append(ac)
                        for ad in fragment_mz:
                            frag_mz_store.append(ad)
                        for ae in intensity:
                            frag_intensity_store.append(ae)
                        for af in resolution:
                            frag_resoltuion_store.append(af)
                        for ag in frag_act_M:
                            frag_act_M_store.append(ag)
                        for ah in theo_frag:
                            frag_theo_mz_store.append(ah)
                        for ai in theo_frag_M:
                            frag_theo_M_store.append(ai)
                        for aj in Da_err:
                            Da_err_store.append(aj)
            
fragment_matches = pd.DataFrame()
fragment_matches['Scan #']= frag_scan_store
fragment_matches['Theoretical Precursor']= frag_theo_precursor_store
fragment_matches['Actual Precursor']= frag_act_precursor_store
fragment_matches['Precursor error (ppm)']= frag_precursor_err_store
fragment_matches['Precursor Actual Charge']= frag_precursor_act_charge_store
fragment_matches['Precursor Theoretical Charge']= frag_precursor_theo_charge_store
fragment_matches['Actual Fragment m/z']= frag_mz_store
fragment_matches['Actual Fragment Charge']= frag_ion_charge_store
fragment_matches['Fragment Intensity']= frag_intensity_store
fragment_matches['MS2 Resolution']= frag_resoltuion_store
fragment_matches['Fragment Actual M']= frag_act_M_store
fragment_matches['Fragment Theoretical m/z']= frag_theo_mz_store
fragment_matches['Fragment Theoretical M']= frag_theo_M_store
fragment_matches['Fragment error (Da)']= Da_err_store   

file_name = sample_name + 'fragment_matches.csv'
file_out_path = final_dir + '\\' + file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        fragment_matches.to_csv(filec,index=False)