# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 17:57:07 2016

@author: rgiessmann
"""


# generate a large data set with random errors

import numpy as np
import pandas
import csv

import footprint
footprinter = footprint.Footprinter()


## PARAMETERS FOR VALIDATION ##
NP_RANDOM_SEED            = 123
THRESHOLD_KD_NO_FOOTPRINT = 1000 #uM
NOISE_PEAK_ABSOLUTE_SD    = 50   #A.U.
NOISE_PEAK_RELATIVE_SD    = 0    #relative A.U.
COEFFICIENT_SD_ACCEPTED   = 2.7
##                           ##


## INTERNAL CONSTANTS ##
NAME_sample_id = "Sample File Name"
NAME_sample_file = "GetFromFile"
##                    ##

def __init__():
    np.random.seed(int(NP_RANDOM_SEED))
    print("Checking seeding with {}".format(int(NP_RANDOM_SEED)))
    print("numpy.random.normal gives:")
    print(np.random.normal())
    print("---")

def get_data(kdmatrix_filename, inputtraces_filename):
    kdmatrix = read_kdmatrix(kdmatrix_filename)
    kdvalues = get_kd_values(kdmatrix)
    ref_trace = generate_averaged_negative_control(kdvalues)
    sample_files = read_sample_traces(inputtraces_filename)
    trace_list = create_traces_and_peaks_with_noise(sample_files, ref_trace, kdvalues)
    return trace_list
    
    
def generate_ref_trace(kdmatrix_filename):
    kdmatrix = read_kdmatrix(kdmatrix_filename)
    kdvalues = get_kd_values(kdmatrix)
    ref_trace = generate_averaged_negative_control(kdvalues)
    return ref_trace

def read_kdmatrix(filename_kdmatrix):
    matrix = pandas.read_csv(filename_kdmatrix)
    return matrix
    
def get_kd_values(matrix):
    kd_series = matrix["KD mean"]
    return kd_series
    
def generate_averaged_negative_control(kd_series):
    how_many_peaks = len(kd_series)
    
    ref = footprint.Trace(file_name = "ref", dye_color = "B", Ltot_conc = 0, Rtot_conc = 0, peaks=[])

    for i in range(how_many_peaks):
        peak_height = np.random.random_integers(50,1000)
        peak_bp = i
        ref.peaks.append(footprint.Peak(peak_bp, peak_height))
    return ref
    
def generate_sample_traces(filename_inputtraces, ref):
    ## NONFUNCTIONAL
    
    trace_list = []
    
    sample_files = read_sample_traces(filename_inputtraces)
    
    for entry in sample_files:
        entry[NAME_sample_id]
        
        
    
    return trace_list
    
def read_sample_traces(filename_inputtraces):
    df = pandas.DataFrame.from_csv(filename_inputtraces)
    sample_files = []
    for index, row in df.iterrows():
        sample_files.append(row.to_dict())
    return sample_files



def write_out_sample_traces(filename_inputtraces, trace_list):
    ## NONFUNCTIONAL
    
    sample_files = read_sample_traces(filename_inputtraces)
    
    all_files = [entry[NAME_sample_file] for entry in sample_files]
    
    for f in all_files:
        pass
        
        
    
def apply_noise_peakheight(height):
    if NOISE_PEAK_ABSOLUTE_SD > 0:
        absolute_noise = np.random.normal(loc=0, scale=NOISE_PEAK_ABSOLUTE_SD)
    else:
        absolute_noise = 0

    if NOISE_PEAK_RELATIVE_SD > 0:
        relative_noise = height * np.random.normal(loc=0, scale=NOISE_PEAK_RELATIVE_SD)
    else:
        relative_noise = 0

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
    
def analyze_success(kd_values_target, kd_matrix_found):
    expect_list = []
    found_list =  []
    success_list = []
    difference_list = []
    sd_list = []
    
    for expected, found in zip(kd_values_target, kd_matrix_found):
        
        if expected > THRESHOLD_KD_NO_FOOTPRINT:
            expect_list.append(False) ## expect: negative
        else:
            expect_list.append(True) ## expect: positive

        if found[1] > THRESHOLD_KD_NO_FOOTPRINT:
            found_list.append(False) ## found: definitely negative
        else:
            found_list.append(True) ## found: positive
 
        diff_mean = abs(expected - found[1])
            
        if diff_mean < found[2] * COEFFICIENT_SD_ACCEPTED:
            success_list.append(True) ## criterion true: value within uncertainty
        else:
            success_list.append(False) ## value without uncertainy
        difference_list.append(diff_mean)
        sd_list.append(found[2])
        
    big_dict = {
        "expected_kd" : list(kd_values_target),
        "found_kd" : list([found[1] for found in kd_matrix_found]),
        "sd_found"      : sd_list,

        "expected_footprint" : expect_list,
        "found_footprint" : found_list,
        "value_within_uncertainty" : success_list,
        "difference_of_means"    : difference_list,

    }
    
    df = pandas.DataFrame.from_dict(big_dict)
    return df
       
def return_binary_classification_success(kd_values_target, p_matrix_found, alpha=0.05):
    expect_list = []
    found_list =  []
    success_list = []
    difference_list = []
    sd_list = []
    
    for expected, (i, found) in zip(kd_values_target, p_matrix_found.iterrows() ):
        
        if expected > THRESHOLD_KD_NO_FOOTPRINT:
            expect_list.append(False) ## condition negative ( = no footprint)
        else:
            expect_list.append(True)  ## condition positive ( = sensitive)

        if found["pvalue"] > alpha:
            found_list.append(False) ## prediction negative ( = no footprint)
        else:
            found_list.append(True)  ## prediction positive ( = sensitive)
 
        
    big_dict = {
        "condition" : expect_list,
        "prediction" : found_list,
    }
    df = pandas.DataFrame.from_dict(big_dict)
    return df

def analyze_binary_classification_success(df):
    confusion_matrix = {}

    sel_condition_positive  = ( df["condition"]  == True  )
    sel_condition_negative  = ( df["condition"]  == False )
    sel_prediction_positive = ( df["prediction"] == True  )
    sel_prediction_negative = ( df["prediction"] == False )

    true_positive = df[sel_condition_positive][sel_prediction_positive]
    true_negative = df[sel_condition_negative][sel_prediction_negative]
    false_positive = df[sel_condition_negative][sel_prediction_positive]
    false_negative = df[sel_condition_positive][sel_prediction_negative]


    total_population = float(len(df))

    P = condition_positive = float(len(df[sel_condition_positive]))
    N = condition_negative = float(len(df[sel_condition_negative]))
    prediction_positive = float(len(df[sel_prediction_positive]))
    prediction_negative = float(len(df[sel_prediction_negative]))



    ## inner fields
    TP = float(len(true_positive))
    FP = float(len(false_positive))
    TN = float(len(true_negative))
    FN = float(len(false_negative))

    ## cross products below
    try:
        PPV = TP / prediction_positive
    except ZeroDivisionError:
        PPV = float("nan")

    try:
        FDR = FP / prediction_positive
    except ZeroDivisionError:
        FDR = float("nan")

    try:
        FOR = FN / prediction_negative
    except ZeroDivisionError:
        FOR = float("nan")

    try:
        NPV = TN / prediction_negative
    except ZeroDivisionError:
        NPV = float("nan")


    ## cross products right
    try:
        TPR = TP / condition_positive
    except ZeroDivisionError:
        TPR = float("nan")

    try:
        FPR = FP / condition_negative
    except ZeroDivisionError:
        FPR = float("nan")

    try:
        FNR = FN / condition_positive
    except ZeroDivisionError:
        FNR = float("nan")

    try:
        TNR = TN / condition_negative
    except ZeroDivisionError:
        TNR = float("nan")


    ## below cross products right 
    try:
        LR_plus = TPR / FPR
    except ZeroDivisionError:
        LR_plus = float("nan")

    try:
        LR_minus = FNR / TNR
    except ZeroDivisionError:
        LR_minus = float("nan")
    
    try:
        DOR = LR_plus / LR_minus
    except ZeroDivisionError:
        DOR = float("nan")




    ## edge bottom left
    try:
        ACC = accuracy = (TP + TN ) / total_population
    except ZeroDivisionError:
        ACC = float("nan")


    ## edge top right
    try:
        prevalence = condition_positive / total_population
    except ZeroDivisionError:
        prevalence = float("nan")

    ## statistics
    try:
        F1 = 2 * PPV * TPR / (PPV + TPR)
    except ZeroDivisionError:
        F1 = float("nan")

    try:
        MCC = (TP*TN - FP*FN) / np.sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) )
    except ZeroDivisionError:
        MCC = float("nan")

    try:
        BM = TPR + TNR - 1
    except ZeroDivisionError:
        BM = float("nan")

    try:
        MK = PPV + NPV - 1
    except ZeroDivisionError:
        MK = float("nan")


    result_dict = {}
    result_dict.update( { "TP": TP } )
    result_dict.update( { "FP": FP } )
    result_dict.update( { "TN": TN } )
    result_dict.update( { "FN": FN } )
    result_dict.update( { "MCC": MCC } )


    return result_dict
        
