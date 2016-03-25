#!/usr/bin/env python

import logging
import sys
import getopt

   

def main(argv=""):

    ## standard settings
    config_file = None
    input_file = None
    peak_list = None
    trace_list = None 

        
    ## argument / start option recognition
    try:
        opts, remaining_args = getopt.getopt(argv,"hc:i:",["help","config-file=","input-file="])
        ## DEBUG: 
        # print(opts,remaining_args)
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
    if remaining_args != []:
        print("You provided too many options. Call me with -h to learn more.")

    ## TODO: check input_files; split   

    ## WARNING: this is step-wise implementing and testing the whole script
    trace_list, peak_list = get_data(None)
    cluster_peaks(trace_list)
    print(calculate_deviance_for_all_peaks(trace_list[0], trace_list[1]))

    ## DEBUG
    #print([peak.cluster for peak in peak_list])

    return


class Trace:
    def __init__(self, file_name, dye_color, Ltot_conc, peaks=[]):
        self.file_name = file_name
        self.dye_color = dye_color
        self.Ltot_conc = Ltot_conc
        self.peaks = peaks
        return
    #def __repr__(self):
    #    return repr({"file_name" : self.file_name, "dye_color" : self.dye_color, "Ltot_conc" : self.Ltot_conc, "peaks" : self.peaks})

class Peak:
    def __init__(self, peak_height):
        self.peak_height = peak_height
        return
    #def __repr__(self):
    #    return repr({"peak_height" : self.peak_height})
    

def get_data(parameters=None):
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
    Trace(file_name = "01-18-16-11-27 AM.fsa", dye_color = "B", Ltot_conc = 0, peaks=[
    Peak(10),
    Peak(20),
    Peak(30)
    ]),
    Trace(file_name = "01-18-16-35-11 AM.fsa", dye_color = "B", Ltot_conc = 5, peaks=[
    Peak(11),
    Peak(21),
    ]),
    ]
    
    ## DEBUG    
    #print(trace_list)
    
    return trace_list, None

def cluster_peaks(trace_list):
    for trace in trace_list:
        i = 1
        for peak in trace.peaks:
            peak.cluster = i
            i += 1
    return 

def calculate_deviance_for_all_peaks(trace, ref, from_bp=20, to_bp=130):
    deviance_for_all_peaks = 0        
    for peak_cluster in set([peak.cluster for peak in ref.peaks]):
        ## DEBUG
        #print(peak_cluster)

        trace_peak = [peak for peak in trace.peaks if peak.cluster == peak_cluster]
        ref_peak = [peak for peak in ref.peaks if peak.cluster == peak_cluster]
        if len(trace_peak)==1 and len(ref_peak)==1:
            deviance_for_all_peaks += ref_peak[0].peak_height - trace_peak[0].peak_height
        else:
            ## WARNING
            print("WARNING!")

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
    #read and append trace list
    #begin loop
      #for each footprinted peak
      #divide peak area with area of biggest of the 0M peaks (at same bp)
      #add result to trace_list (fR)
    #end loop
    print("")
    return fractional_occupancies

def calculate_free_ligand_concentration():
    #read and append trace list
    #begin loop
      #for each L(total)
      #calculate: L(free)= L(total)-fR(1)*R(total)-fR(2)*R(total)-...-fR(n)*R(total) bzw. L(free)= L(total)-R(total)*(fR(1)+fR(2)+...+fR(n))
      #add result to trace list (L(free))
    #end loop
    print("")
    return fractional_occupancies


def fit_data_determine_kd():
    #read an append trace list
    #begin loop
      #for each footprinting site
      #fit fR(n)=L(free)/(Kd(n)+L(free))
      #add result to trace list
    #end loop
    print("")
    return kd_values

def plot_data():
    #plot fR vs L(free)
    #plot fit
    #show Kd and Basepair number
    print("")
    return plots


if __name__ == "__main__":
  main(sys.argv[1:])
