#!/usr/bin/env python

import logging
import sys
import getopt
import numpy
import scipy

   

def main(argv=""):

    ## standard settings
    config_file = None
    input_file = None
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
    trace_list = get_data(None)
    trace,ref=trace_list[0], trace_list[1]
    cluster_peaks(trace_list)
    print("uncorrected RMSD of peaks: "+str(calculate_deviance_for_all_peaks(ref,trace)))
    #for i,j in give_all_clustered_peaks(trace_list[0], trace_list[1]): print(i.peak_height,j.peak_height)
    print("optimal factor : "+str(determine_factor_numerically(trace_list[0], trace_list[1])))
    add_fractional_occupancies(ref,trace)
    print("Lfree : "+str(calculate_free_ligand_concentration(ref,trace)))
    ## DEBUG
    #print([peak.cluster for peak in peak_list])

    print("done.")

    return


class Trace:
    def __init__(self, file_name, dye_color, Ltot_conc, Rtot_conc, peaks=[]):
        self.file_name = file_name
        self.dye_color = dye_color
        self.Ltot_conc = Ltot_conc
        self.Rtot_conc = Rtot_conc
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
    Trace(file_name = "01-18-16-11-27 AM.fsa", dye_color = "B", Ltot_conc = 5, Rtot_conc = 0.1, peaks=[
    Peak(10),
    Peak(40),
    Peak(60)
    ]),
    Trace(file_name = "01-18-16-35-11 AM.fsa", dye_color = "B", Ltot_conc = 0, Rtot_conc = 0.1, peaks=[
    Peak(20),
    Peak(40),
    Peak(60)    
    ]),
    ]
    
    ## DEBUG    
    #print(trace_list)
    
    return trace_list

def cluster_peaks(trace_list):
    ## TODO: make this skeleton function come alive

    for trace in trace_list:
        i = 1
        for peak in trace.peaks:
            peak.cluster = i
            i += 1
    return 

def calculate_deviance_for_all_peaks(ref, trace, from_bp=20, to_bp=130):
    '''calculates the RMSD for peaks that were identified as clustered in trace _ref_, compared to _trace_, in the range (from_bp, to_bp)'''

    deviance_for_all_peaks = 0        
    n=0
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
        #print(ref_peak,trace_peak)        
        
        ## SIMPLE VERSION        
        #deviance_for_all_peaks += (ref_peak.peak_height - trace_peak.peak_height)**2
        
        ## FANCY VERSION
        if ref_peak.peak_height > trace_peak.peak_height:
            deviance_for_all_peaks += (ref_peak.peak_height - trace_peak.peak_height)**2
        else:
            ## deviance_for_all_peaks += 0
            ## equals 
            pass

        n+=1

    rmsd = numpy.sqrt(deviance_for_all_peaks/n)

    return rmsd
    
def give_all_clustered_peaks(ref,trace):
    for peak_cluster in set([peak.cluster for peak in ref.peaks]):
        trace_peak = [peak for peak in trace.peaks if peak.cluster == peak_cluster]
        ref_peak = [peak for peak in ref.peaks if peak.cluster == peak_cluster]        
        if len(trace_peak)==1 and len(ref_peak)==1:
            yield (ref_peak[0],trace_peak[0])
                        

def determine_factor_numerically(ref, trace):
    ## TODO: implement real optimizer via scipy.optimize()    
    
    optimal_factor = 1
    rmsd_old = calculate_deviance_for_all_peaks(ref,trace)

    ## store original data
    for peak in trace.peaks:
        peak.peak_height_original = peak.peak_height


    # Robert's approach
    #define reference
    ## -> already defined by function call

    #loop1 begin
     #for each trace
    ## -> actually, we will loop only through this given set of two traces, the other looping is done in the parent function

     #loop2 begin
       #use factors 0 to 3.5 in 0.01 steps , calculate new peak areas from peak areas
    for factor in numpy.arange(0,3.5,0.01):

        correct_peaks_with_factor(trace,factor)
            
        #use calculate_deviance_for_all_peaks with trace and ref
        rmsd_new = calculate_deviance_for_all_peaks(ref,trace)
                
        ## restore all peak heights to original height
        for peak in trace.peaks:
            peak.peak_height = peak.peak_height_original 
        
        ## DEBUG        
        #list deviance_new (including factor)
        #print(str(factor)+" : "+str(rmsd_new))

        #compare deviance_new with deviance_old, if better deviance_new -> deviance_old, else delete deviance_new

        if rmsd_new < rmsd_old:
            rmsd_old = rmsd_new
            optimal_factor = factor
            
            ##DEBUG
            #print("new best")

        #end loop2
        #put factor from deviance_old to trace_list (for right trace)
        ## -> this factor is simply returned, let parent function decide what to do with it!        
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
    print("")
    return


def correct_peaks_with_factor(trace, factor):
    #read and append peak_list
    #multiply area of each peak for all traces with right factor from trace_list 
    ## -> handled with in parent function
    #add new area to peak_list
    
    for peak in trace.peaks:
        peak.peak_height = peak.peak_height * factor
    return trace

def which_peaks_differ(threshold=0.10):
    #begin loop
      #for each peak 
      #compare peak areas of traces with different concentration
      #if difference >0.1 add peak to trace list
    #end loop
    print("")
    return peak_list

def add_fractional_occupancies(ref,trace):
    #read and append trace list
    #begin loop
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
        #for each footprinted peak
        trace_peak.fractional_occupancy = trace_peak.peak_height / ref_peak.peak_height        
        #divide peak area with area of biggest of the 0M peaks (at same bp)
        #add result to trace_list (fR)
    
    #end loop
    print("")
    return 

def calculate_free_ligand_concentration(ref,trace):
    #read and append trace list
    #begin loop
      #for each L(total)
      #calculate: L(free)= L(total)-fR(1)*R(total)-fR(2)*R(total)-...-fR(n)*R(total) bzw. L(free)= L(total)-R(total)*(fR(1)+fR(2)+...+fR(n))
    sum_fractional_occupancies = 0
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
        #for each footprinted peak
        if trace_peak.fractional_occupancy < 1:
            sum_fractional_occupancies += trace_peak.fractional_occupancy
        
    ## DEBUG
    print(sum_fractional_occupancies)        
        
    Lfree_conc = trace.Ltot_conc - trace.Rtot_conc*sum_fractional_occupancies
    #add result to trace list (L(free))
    #end loop
    print("")
    return Lfree_conc


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
