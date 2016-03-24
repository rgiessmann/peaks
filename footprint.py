#!/usr/bin/env python

import logging
import sys
import getopt

def main(argv=""):
    config_file = None
    try:
        opts, remaining_args = getopt.getopt(argv,"hc:i:",["help","config-file=","input-file="])
    except getopt.GetoptError:
       print('You provided unusual arguments. Call me with -h to learn more.')
       sys.exit(2)
    for opt, arg in opts:
       if opt in ('-h', '--help'):
           ## TODO
           print("...help")
       if opt in ('-c','--config-file'):
           config_file = arg
       if opt in ('-i','-input-file'):
           input_file = arg
    input_file += remaining_args

    ## TODO: check input_files; split   

    return


class Trace:
    def __init__(self, file_name, dye_color, Ltot_concentration):
        self.file_name = file_name
        self.dye_color = dye_color
        self.Ltot_concentration = Ltot_concentration
        return
    def __repr__(self):
        return repr({"file_name" : self.file_name, "dye_color" : self.dye_color, "Ltot_concentration" : self.Ltot_concentration})

class Peak:
    def __init__(self, peak_height):
        self.peak_height = peak_height
        return
    def PeakHeight(self):
        return self.peak_height
    def __repr__(self):
        return
    

def get_data(parameters):
    ## WARNING : this is a non-functional skeleton function

    #read csv (?) file
    #if color blue and if peak size between ... and ..., copy size, name and area to peak_list

    #trace_list = ... #contains factor for each trace, footprint peaks, fractional occupancies, ligand and receptor concentrations, kd values
    #peak_list = ... #contains Peaks and areas, new calculated areas

    ## TODO: read in the data from config and input files
    ## trace_list = ... #contains factor for each trace, footprint peaks, fractional occupancies, ligand and receptor concentrations, kd values
    ## peak_list = ... #contains Peaks and areas, new calculated areas

    ## 1. create minimal data objects with classes
    trace_list = [
    Trace(file_name = "01-18-16-11-27 AM.fsa", dye_color = "B", Ltot_conc = 0),
    Trace(file_name = "01-18-16-35-11 AM.fsa", dye_color = "B", Ltot_conc = 5),
    ]
    
    peak_list = [
    Peak(),
    Peak(),    
    ]

    return trace_list, peak_list

def cluster_peaks(peak_list):
    for peak in peak_list:
        peak.cluster = 1
    return 

def calculate_deviance_for_all_peaks(from_bp, to_bp, trace, ref):
    deviance_for_all_peaks = 0        
    for peak_cluster in set([peak.cluster for peak in peak_list]):
        trace_peak = [peak for peak in trace if peak.cluster == peak_cluster]
        ref_peak = [peak for peak in ref if peak.cluster == peak_cluster]
        if trace_peak and ref_peak:
            deviance_for_all_peaks += ref_peak.PeakHeight() - trace_peak.PeakHeight()
    #deviance_for_all_peaks = 1
    return deviance_for_all_peaks

def determine_factor_numerically(ref, trace):
    # Robert's approach
    #define reference
    #loop1 begin
     #for each trace
     #loop2 begin
       #use factors 0 to 3.5 in 0.01 steps , calculate new peak areas from peak areas
       #use calculate_deviance_for_all_peaks with trace and ref
       #list deviance_new (including factor)
       #compare deviance_new with deviance_old, if better deviance_new -> deviance_old, else delete deviance_new
     #end loop2
     #put factor from deviance_old to trace_list (for right trace)
    #end loop1

    optimal_factor = 1
    return optimal_factor

def determine_factor_single_peak():
    #define reference
        #loop1 begin
     #for each trace
     #loop2 begin
       #calculate factor for each Peak to bring to size of reference peak
       #calculate new areas for every peak with factor 
       #calculate deviance between each peak and reference
       #list deviance_new (including factor)
       #compare deviance_new with deviance_old, if better deviance_new -> deviance_old, else delete deviance_new
     #end loop2
     #put factor from deviance_old to trace_list (for right trace)
    #end loop1
    print("")
    return


def correct_peaks_with_factor(trace, factor):
    #read and append peak_list
    #multiply area of each peak for all traces with right factor from trace_list 
    #add new area to peak_list
    print("")
    return

def which_peaks_differ(threshold=0.10):
    #begin loop
      #for each peak 
      #compare peak areas of traces with different concentration
      #if difference >0.1 add peak to trace list
    #end loop
    print("")
    return peak_list

def calculate_fractional_occupancies(peak_list):
    print("")
    return fractional_occupancies

def calculate_free_ligand_concentration():
    print("")
    return fractional_occupancies


def fit_data_determine_kd():
    print("")
    return kd_values

def plot_data():
    print("")
    return plots

