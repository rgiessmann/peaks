#!/usr/bin/env python
#asdf
import logging                                
import sys                                
import getopt                                

def main(argv=""):   
  try:
    opts, remaining_args = getopt.getopt(argv,"dhia:gf:o:k",["help","iterate","iterate-from=","iterate-to=","iter-step=","align-to=","down-pose","both-poses","generate-all-dummies","generate-dummies-from=","generate-dummies-to=","finegrain-step=","keep-dummies","dna-leading-strand=","dna-lagging-strand=","dna-leading-start=","dna-lagging-start=","bp-counting-parallel","minor-groove-offset=","output-prefix="])
  except getopt.GetoptError:
    print 'You provided unusual arguments. Call me with -h to learn more.'
    sys.exit(2)
  for opt, arg in opts:
    if opt in ('-h', '--help'):
...
return


## data structure

trace_list = ...
peak_list = ...

def get_data(parameter):
    
return

def calculate_deviance_for_all_peaks(from, to, trace, ref):
    ...
return deviance_for_all_peaks

def determine_factor_numerically(ref, trace, ):
    # Robert's approach
return optimal_factor

def determine_factor_single_peak():
    
return


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

