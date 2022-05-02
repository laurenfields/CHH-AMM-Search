# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 11:56:58 2022

@author: lawashburn
"""

import os
import csv
import pandas as pd
import numpy as np
from datetime import datetime
now = datetime.now()

fragment_matches = pd.read_csv(r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\output_directory\SG-CHH-Fraction-UnoptimizedMethod_NewLCgradient_1fragment_matches.csv")
by_ion_list = pd.read_csv(r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\by_ions.csv")
final_dir = r"C:\Users\lawashburn\Documents\Nhu_Prescursor_Matching\20220417_oldprecursor\output_directory"

sample_name = 'SG-CHH-Fraction-UnoptimizedMethod_NewLCgradient_1'
fragment_charges = [1,2,3,4,5]

scan_no = []
t_p = []
a_p = []
p_err = []
p_a_z = []
p_t_z = []
a_f_mz = []
a_f_z = []
f_intensity = []
resolution_archive = []
f_a_M = []
f_t_mz = []
f_t_M = []
f_err = []
species_archive = []
f_t_mz2 = []
f_t_z = []

for a in fragment_charges:
    fragments_filtered = pd.DataFrame()
    fragments_filtered['species'] = by_ion_list['peptide']
    fragments_filtered[str(a)] = by_ion_list[str(a)]
    
    fragments_merged = fragment_matches.merge(fragments_filtered,left_on='Fragment Theoretical m/z',right_on=[str(a)])
    fragments_merged['Theoretical ion charge'] = a
    scan = fragments_merged['Scan #'].values.tolist()
    t_prec = fragments_merged['Theoretical Precursor'].values.tolist()
    a_prec = fragments_merged['Actual Precursor'].values.tolist()
    prec_err = fragments_merged['Precursor error (ppm)'].values.tolist()
    a_prec_z = fragments_merged['Precursor Actual Charge'].values.tolist()
    t_prec_z = fragments_merged['Precursor Theoretical Charge'].values.tolist()
    a_frag_mz = fragments_merged['Actual Fragment m/z'].values.tolist()
    a_frag_z = fragments_merged['Actual Fragment Charge'].values.tolist()
    intensity = fragments_merged['Fragment Intensity'].values.tolist()
    resolution = fragments_merged['MS2 Resolution'].values.tolist()
    a_frag_M = fragments_merged['Fragment Actual M'].values.tolist()
    t_frag_mz = fragments_merged['Fragment Theoretical m/z'].values.tolist()
    t_frag_M = fragments_merged['Fragment Theoretical M'].values.tolist()
    frag_err = fragments_merged['Fragment error (Da)'].values.tolist()
    species = fragments_merged['species'].values.tolist()
    t_frag_mz_2 = fragments_merged[str(a)].values.tolist()
    t_z_z = fragments_merged['Theoretical ion charge'].values.tolist()
    
    for b in scan:
        scan_no.append(b)
    
    for c in t_prec:
        t_p.append(c)
    
    for d in a_prec:
        a_p.append(d)
        
    for e in prec_err:
        p_err.append(e)
    
    for f in a_prec_z:
        p_a_z.append(f)
    
    for g in t_prec_z:
        p_t_z.append(g)
    
    for h in a_frag_mz:
        a_f_mz.append(h)
    
    for i in a_frag_z:
        a_f_z.append(i)
    
    for j in intensity:
        f_intensity.append(j)
    
    for k in resolution:
        resolution_archive.append(k)
    
    for l in a_frag_M:
        f_a_M.append(l)
    
    for m in t_frag_mz:
        f_t_mz.append(m)
    
    for n in t_frag_M:
        f_t_M.append(n)
    
    for o in frag_err:
        f_err.append(o)
    
    for p in species:
        species_archive.append(p)
    
    for q in t_frag_mz_2:
        f_t_mz2.append(q)
    
    for t in t_z_z:
        f_t_z.append(t)

match_table = pd.DataFrame()
match_table['Species'] = species_archive
match_table['Scan #'] = scan_no
match_table['theoretical precursor'] = t_p
match_table['actual precursor'] = a_p
match_table['theoretical precursor charge'] = p_t_z
match_table['actual precursor charge'] = p_a_z
match_table['precursor error (ppm)'] = p_err
match_table['theoretical fragment m/z'] = f_t_mz
match_table['actual fragment m/z'] = a_f_mz
match_table['theoretical fragment charge'] = f_t_z
match_table['actual fragment charge'] = a_f_z
match_table['theoretical fragment M'] = f_t_M
match_table['actual fragment M'] = f_a_M
match_table['fragment error (Da)'] = f_err
match_table['fragment intensity'] = f_intensity
match_table['MS2 resolution'] = resolution_archive

match_table = match_table[match_table['theoretical fragment charge'] == match_table['actual fragment charge']]  

species_reference = match_table['Species'].values.tolist()

species_filtered = []
for z in species_reference:
    if z not in species_filtered:
        species_filtered.append(z)
       
for y in species_filtered:
    scan_summary_no = []
    scan_summary_instances = []
    
    species_matches = match_table[match_table['Species'] == y]  
    species_matches_no_dups = species_matches.drop_duplicates(subset = ['Scan #','theoretical precursor','actual precursor','theoretical precursor charge','actual precursor charge',
                                                                'precursor error (ppm)','theoretical fragment m/z','actual fragment m/z','theoretical fragment charge',
                                                                'actual fragment charge','theoretical fragment M','actual fragment M','fragment error (Da)',
                                                                'fragment intensity','MS2 resolution'])
    if len(species_matches_no_dups)>0:
        file_name = sample_name + '_' + y + '_fragment_matches.csv'
        file_out_path = final_dir + '\\' + file_name
        with open(file_out_path,'w',newline='') as filec:
                writerc = csv.writer(filec)
                species_matches_no_dups.to_csv(filec,index=False)
        
        scans_present = species_matches_no_dups['Scan #'].values.tolist()
        
        scans_filtered = []
        for u in scans_present:
            if u not in scans_filtered:
                scans_filtered.append(u)
        
        for v in scans_filtered:
            scan_count = scans_present.count(v)
            scan_summary_no.append(v)
            scan_summary_instances.append(scan_count)
    
    scan_report = pd.DataFrame()
    scan_report['Scan #'] = scan_summary_no
    scan_report['# instances'] = scan_summary_instances
    
    file_name = sample_name + '_' + y + '_scan_instances.csv'
    file_out_path = final_dir + '\\' + file_name
    with open(file_out_path,'w',newline='') as filec:
            writerc = csv.writer(filec)
            scan_report.to_csv(filec,index=False)           