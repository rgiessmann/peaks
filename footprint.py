#!/usr/bin/env python

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

trace_list = ... #contains factor for each trace, footprint peaks, fractional occupancies, ligand and receptor concentrations, kd values
peak_list = ... #contains Peaks and areas, new calculated areas

def get_data(parameter):
    #read csv (?) file
    #if color blue and if peak size between ... and ..., copy size, name and area to peak_list
return

def calculate_deviance_for_all_peaks(from, to, trace, ref):
    ...
return deviance_for_all_peaks

def determine_factor_numerically(ref, trace, ):
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
return


def correct_peaks_with_factor(trace, factor):
    #read and append peak_list
    #multiply area of each peak for all traces with right factor from trace_list 
    #add new area to peak_list
return

def which_peaks_differ(threshold=0.10):
    #begin loop
      #for each peak 
      #compare peak areas of traces with different concentration
      #if difference >0.1 add peak to trace list
    #end loop
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

