# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 11:54:58 2022

@author: lawashburn
"""

import pandas as pd
import csv
import os

zero_dir = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\06assigned_precursor_matches\0assign\working_directory"
normal_df = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\06assigned_precursor_matches\non0\Untarget_Brain1_20ppmfragment_matches.csv"
output_dir = r"C:\Users\lawashburn\Documents\HyPep1.0\HyPep_Simple_ASMS_Results\20220429\06assigned_precursor_matches\combined"
sample_name = 'Untarget_Brain1'

def get_file_names_with_strings(str_list):
    full_list = os.listdir(zero_dir)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]
    return final_list

file_query = '_zero_reassign_fragment_matches.csv'
fragment_list = (get_file_names_with_strings([file_query]))

all_fragments = pd.DataFrame()

for a in fragment_list:
    a_path = zero_dir + '\\' + a
    a_df = pd.read_csv(a_path)
    
    all_fragments = all_fragments.append(a_df)

normal_df = pd.read_csv(normal_df)
all_fragments = all_fragments.append(normal_df)

print(all_fragments)

file_name = sample_name + '_all_fragments.csv'
file_out_path = output_dir + '\\' + file_name
with open(file_out_path,'w',newline='') as filec:
        writerc = csv.writer(filec)
        all_fragments.to_csv(filec,index=False)