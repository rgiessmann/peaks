#!/usr/bin/env python
#asdf
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
...

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
    def __init__(self):
        return
    def __repr__(self):
        return
    

def get_data(parameters):
    ## WARNING : this is a non-functional skeleton function

    ## TODO: read in the data from config and input files

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

def calculate_deviance_for_all_peaks(from_bp, to_bp, trace, ref):
    deviance_for_all_peaks = 1
return deviance_for_all_peaks

def determine_factor_numerically(ref, trace):
    # Robert's approach
    optimal_factor = 1
return optimal_factor

def determine_factor_single_peak():
    optimal_factor=1
return optimal_factor


def correct_peaks_with_factor(trace, factor):
    ...
return

def which_peaks_differ(threshold=0.10):
    ...
return peak_list

def calculate_fractional_occupancies(peak_list):
    ...
return fractional_occupancies

def calculate_free_ligand_concentration():
    ...
return fractional_occupancies


def fit_data_determine_kd():
    ...
return kd_values

def plot_data():
    ...
return plots

