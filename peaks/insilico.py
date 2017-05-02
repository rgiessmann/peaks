# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 17:57:07 2016

@author: rgiessmann
"""


# generate a large data set with random errors

import numpy as np
import peaks.footprint
import pandas
import csv
footprint = peaks.footprint.Footprinter()


## PARAMETERS FOR VALIDATION ##
NP_RANDOM_SEED            = 123
THRESHOLD_KD_NO_FOOTPRINT = 1000 #uM
NOISE_PEAK_ABSOLUTE_SD    = 10   #A.U.
NOISE_PEAK_RELATIVE_SD    = 1e-3 #relative A.U.
COEFFICIENT_SD_ACCEPTED   = 2.7
##                           ##


## INTERNAL CONSTANTS ##
NAME_sample_id = "Sample File Name"
NAME_sample_file = "GetFromFile"
##                    ##

def __init__():
    np.random.seed(NP_RANDOM_SEED)
    

def generate_dataset(kdmatrix_filename):
    kdmatrix = read_kdmatrix(kdmatrix_filename)
    kdvalues = get_kd_values(kdmatrix)
    ref_trace = generate_averaged_negative_control(kdvalues)
    create_traces_and_peaks_with_noise(sample_files, ref, kdvalues)
    

def read_kdmatrix(filename_kdmatrix):
    matrix = pandas.read_csv(filename_kdmatrix)
    return matrix
    
def get_kd_values(matrix):
    kd_series = matrix["KD mean"]
    return kd_series
    
def generate_averaged_negative_control(kd_series):
    how_many_peaks = len(kd_series)
    
    ref = footprint.Trace(file_name = "ref", dye_color = "B", Ltot_conc = 0, Rtot_conc = rtot, peaks=[])

    for i in range(how_many_peaks):
        peak_height = np.random.random_integers(50,1000)
        peak_bp = i
        ref.peaks.append(footprint.Peak(peak_bp, peak_height))
    
    return ref
    
def generate_sample_traces(filename_inputtraces, ref):
    trace_list = []
    
    sample_files = read_sample_traces(filename_inputtraces)
    
    for entry in sample_files:
        entry[NAME_sample_id]
        
        
    
    return trace_list
    
def read_sample_traces(filename_inputtraces):
    
    with open(filename_inputtraces) as g:
        csv_reader = csv.DictReader(g)
        sample_files = []
        for row in csv_reader:
            sample_files.append(row)
    return sample_files



def write_out_sample_traces(filename_inputtraces, trace_list):
    sample_files = read_sample_traces(filename_inputtraces)
    
    all_files = [entry[NAME_sample_file] for entry in sample_files]
    
    for f in all_files:
        !rm {f}
        !touch {f}
        
        
    
def apply_noise_peakheight(height):
    absolute_noise = np.random.normal(loc=0, scale=NOISE_PEAK_ABSOLUTE_SD)
    relative_noise = height * np.random.normal(loc=0, scale=NOISE_PEAK_RELATIVE_SD)
    
    height = height + absolute_noise + relative_noise
    
    return height
    
def give_footprinted_peak_height(ref_height, lconc, kd_value):
    correct_height = ref_height * (1 - lconc/(lconc+kd_value))
    return correct_height
        
def create_traces_and_peaks_with_noise(sample_files, ref, kd_series): 
    
    trace_list = []
    
    for data in sample_files:

        filename = data["Sample File Name"]
        dye = data["Dye"]
        ltot = data["Ltot"]
        rtot = data["Rtot"]
        
        trace_list.append(footprint.Trace(file_name = filename, dye_color = dye, Ltot_conc = ltot, Rtot_conc = rtot, peaks=[]))
    
        for ref_peak, kd_value in zip(ref.peaks, kd_series):
            if kd_value < THRESHOLD_KD_NO_FOOTPRINT:
                peak_height = give_footprinted_peak_height(ref_peak.peak_height, ltot, kd_value)
            else:
                peak_height = ref_peak.peak_height
 
            peak_height = apply_noise_peakheight(peak_height)
            peak_bp     = ref_peak.size_bp

            trace_list[-1].peaks.append(footprint.Peak(peak_bp, peak_height))
        
    return trace_list
    
def analyze_success(kd_matrix_target, kd_matrix_found):
    success_list = []
    difference_list = []
    sd_list = []
    
    for expected, found in zip(kd_matrix_target, kd_matrix_found):
        diff_mean = abs(expected - found["KD mean"])
        if diff_mean < found["KD SD"] * COEFFICIENT_SD_ACCEPTED:
            success_list.append(True)
        else:
            success_list.append(False)
        difference_list.append(diff_mean)
        sd_list.append(found["KD SD"])
        
    big_dict = {
        "success" : success_list,
        "diff"    : difference_list,
        "sd"      : sd_list
    }
    
    return big_dict
    
        