#!/usr/bin/env python

import logging
import sys
import getopt
import numpy
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
   

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

 ## -> find test cases in test/ directory!
    
    print("done.")

    return


class Trace:
    def __init__(self, file_name, dye_color, Ltot_conc, Rtot_conc, peaks=[]):
        self.file_name = file_name
        self.dye_color = dye_color
        self.Ltot_conc = float(Ltot_conc)
        self.Rtot_conc = float(Rtot_conc)
        self.peaks = peaks
        return
    def __repr__(self):
        return repr(vars(self))
    #    return repr({"file_name" : self.file_name, "dye_color" : self.dye_color, "Ltot_conc" : self.Ltot_conc, "peaks" : self.peaks})

class Peak:
    def __init__(self, size_bp, peak_height):
        self.size_bp = float(size_bp)
        self.peak_height = float(peak_height)
        return
    def __repr__(self):
        return repr(vars(self))
    #    return repr({"peak_height" : self.peak_height})

class Index:
    pass    

def list_traces(read_filelist="../HexA.csv"):
    import csv
    if type(read_filelist) is not list:
            read_filelist = [read_filelist]
    storage_traces = []
    for file in read_filelist:
        csv_reader = csv.reader(open(file))
        header = csv_reader.__next__()
        index = Index()
        index.file_name = header.index('Sample File Name')
        index.sample_name = header.index('Sample Name')
        
        for row in csv_reader:
            if row[index.file_name] not in [t[0] for t in storage_traces]:
                # Nope.
                storage_traces.append([row[index.file_name],row[index.sample_name]])

    print("Writing read traces to output_traces.csv...")
    w=csv.writer(open("output_traces.csv","w"))
    w.writerow(["Sample ID","Sample File Name", "Dye", "Ltot", "Rtot"])
    for row in storage_traces:
        w.writerow([row[1],row[0],"?","?","?"])
    return storage_traces

def get_data(read_filelist="../HexA.csv"):
    ## WARNING : this is a non-functional skeleton function
    ## TODO: read in the data from config and input files


    #read csv (?) file
    #if color blue and if peak size between ... and ..., copy size, name and area to peak_list

    #trace_list = ... #contains factor for each trace, footprint peaks, fractional occupancies, ligand and receptor concentrations, kd values
    #peak_list = ... #contains Peaks and areas, new calculated areas


    ## 1. create minimal data objects with classes
    if read_filelist == None:
        trace_list = [
        Trace(file_name = "01-18-16-11-27 AM.fsa", dye_color = "B", Ltot_conc = 5, Rtot_conc = 0.1, peaks=[
        Peak(1.07,10),
        Peak(1.98,39),
        Peak(3.01,61)
        ]),
        Trace(file_name = "01-18-16-35-11 AM.fsa", dye_color = "B", Ltot_conc = 0, Rtot_conc = 0.1, peaks=[
        Peak(1,20),
        Peak(2,40),
        Peak(3,60)    
        ]),
        ]
    else:
        import csv
        trace_list = []
        ## TODO: what entries to accept?
        config_file = "../input_traces.csv"
        with open(config_file) as g:
            csv_reader = csv.DictReader(g)
            sample_files = []
            for row in csv_reader:
                sample_files.append(row)
                  
        sample_file_names = []
        for i in sample_files:
            sample_file_names.append(i["Sample File Name"]) 
        
        if type(read_filelist) is not list:
            read_filelist = [read_filelist]

        storage_traces = []
        for file in read_filelist:
            with open(file) as f:            
                csv_reader = csv.reader(f)
                header = csv_reader.__next__()
                index = Index()
                index.peak_height = header.index("Height")
                index.size_bp = header.index("Size")
                index.file_name = header.index('Sample File Name')
                index.sample_name = header.index('Sample Name')
                index.dye = len(header)
                index.sample_peak = len(header)+1
                
                for row in csv_reader:
                    #rules for entry acceptance?
                    
                    if True: #split_combined_fields == 
                        row.extend(row[header.index('Dye/Sample Peak')].split(","))
                    if row[index.file_name] in sample_file_names and "B" in row[index.dye]:
                        # trace already in trace_list?                                        
                        if row[index.file_name] not in storage_traces:
                            # Nope.
                            ## load additional data
                            data = [sample for sample in sample_files if sample["Sample File Name"] == row[index.file_name]][0]
                            trace_list.append(
                            Trace(file_name = data["Sample File Name"], dye_color = data["Dye"], Ltot_conc = data["Ltot"], Rtot_conc = data["Rtot"], peaks=[])                    
                            )
                            storage_traces.append(row[index.file_name])
                        
                        if row[index.size_bp] == "" or row[index.peak_height] == "":
                            continue
                        
                        #to which trace?
                        t = [trace for trace in trace_list if trace.file_name == row[index.file_name]][0]
                        t.peaks.append(Peak(row[index.size_bp], row[index.peak_height]))
                        del(t)            

        for foo in sample_file_names:
            if foo not in storage_traces:
                print("Couldn't find trace "+str(foo)+" in the given files...")
            
                
                
    
    ## DEBUG    
    #print(trace_list)
    
    return trace_list

def cluster_peaks(ref,trace_list,accepted_offset=0.25):
    ## TODO: make this skeleton function come alive

    #for trace in trace_list:
    #    i = 1
    #    for peak in trace.peaks:
    #        peak.cluster = i
    #        i += 1

    i = 0            
    for ref_peak in ref.peaks:
        i += 1
        ref_peak.cluster = i
        for trace in trace_list:
            for peak in trace.peaks:
                if ref_peak.size_bp - accepted_offset < peak.size_bp < ref_peak.size_bp + accepted_offset:
                    peak.cluster = i
    for trace in trace_list:
        for peak in trace.peaks:
            if "cluster" not in vars(peak):
                peak.cluster = 0
    
    return 

def calculate_deviance_for_all_peaks(ref, trace, weight_twoside=1,weight_oneside=0,from_bp=20, to_bp=130,):
    '''calculates the RMSD for peaks that were identified as clustered in trace _ref_, compared to _trace_, in the range (from_bp, to_bp)'''

    deviance_for_all_peaks = 0
    deviance_for_every_peak = 0        
    n=0
    m=0
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
        trace_peak=trace_peak[0]
        
        if ref_peak.peak_height > trace_peak.peak_height:
            # deviation for potentially footprinted peaks
            deviance_for_all_peaks += (ref_peak.peak_height - trace_peak.peak_height)**2
            n+=1
        else:
            # deviation for potentially hypersensitive peaks
            deviance_for_every_peak += (ref_peak.peak_height - trace_peak.peak_height)**2
            m+=1

    #rmsd = numpy.sqrt(deviance_for_all_peaks/n)
    rmsd_all = numpy.sqrt((deviance_for_all_peaks + deviance_for_every_peak)/(m+n))

    #print(rmsd_all)    

    ## decide what to evaluate here -- maybe something weighted, maybe not?
    value = weight_twoside*rmsd_all + weight_oneside*rmsd
    
    return value
    
def give_all_clustered_peaks(ref,trace_list):

    if type(trace_list) is not list:
        trace_list = [trace_list]
    for peak_cluster in set([peak.cluster for peak in ref.peaks]):
        ref_peak = [peak for peak in ref.peaks if peak.cluster == peak_cluster]        
        trace_peaks=[]   
        for trace in trace_list:
            foo=[peak for peak in trace.peaks if peak.cluster == peak_cluster]
            if len(foo) == 1:
                trace_peaks.append(foo[0])
        if len(ref_peak)==1:
            if len(trace_peaks) != len(trace_list):
                print("WARNING: Unable to find clustered peaks for all given traces...")
            yield (ref_peak[0],trace_peaks)
                        
def determine_factor_numerically(ref, trace,weight_twoside=1,weight_oneside=0):
    ## TODO: implement real optimizer via scipy.optimize()    
    
    optimal_factor = 1
    rmsd_old = calculate_deviance_for_all_peaks(ref,trace,weight_twoside,weight_oneside)

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
        rmsd_new = calculate_deviance_for_all_peaks(ref,trace,weight_twoside,weight_oneside)
                
        ## restore all peak heights to original height
        for peak in trace.peaks:
            peak.peak_height = peak.peak_height_original 
        
        ## DEBUG        
        #list deviance_new (including factor)
        print(str(factor)+" : "+str(rmsd_new))

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

def determine_factor_single_peak(ref, trace):
   
    optimal_factor = 1
    rmsd_old = calculate_deviance_for_all_peaks(ref,trace)

    ## store original data
    for peak in trace.peaks:
        peak.peak_height_original = peak.peak_height
        
    #define reference
        #loop1 begin
     #for each trace
     #loop2 begin
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
        factor =  ref_peak.peak_height / trace_peak.peak_height 
        correct_peaks_with_factor(trace,factor)
            
        #use calculate_deviance_for_all_peaks with trace and ref
        rmsd_new = calculate_deviance_for_all_peaks(ref,trace)
     
        ## restore all peak heights to original height
        for peak in trace.peaks:
            peak.peak_height = peak.peak_height_original 
        
        ## DEBUG        
        #list deviance_new (including factor)
        print(str(factor)+" : "+str(rmsd_new))

        #compare deviance_new with deviance_old, if better deviance_new -> deviance_old, else delete deviance_new

        if rmsd_new < rmsd_old:
            rmsd_old = rmsd_new
            optimal_factor = factor
            
    print("")
    return optimal_factor


def correct_peaks_with_factor(trace, factor):
    #read and append peak_list
    #multiply area of each peak for all traces with right factor from trace_list 
    ## -> handled with in parent function
    #add new area to peak_list
    
    for peak in trace.peaks:
        peak.peak_height = peak.peak_height * factor
    return trace

def which_peaks_differ(ref, trace):
    #begin loop
      #for each peak 
      #compare peak areas of traces with different concentration
      #if difference >0.1 add peak to trace list
    #end loop
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
         if ref_peak.peak_height - trace_peak.peak_height > 0.1*trace_peak.peak_height and ref.Ltot_conc != trace.Ltot_conc:
           trace_peak.footprinted_peak = 1
         else:
           trace_peak.footprinted_peak = 0
        
    print("")
    return 

def add_fractional_occupancies(ref,trace):
    #read and append trace list
    #begin loop
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
        #for each footprinted peak
        trace_peak.fractional_occupancy = 1 - trace_peak.peak_height / ref_peak.peak_height        
        #divide peak area with area of biggest of the 0M peaks (at same bp)
        #add result to trace_list (fR)
    fractional_occupancy = trace_peak.fractional_occupancy
    #end loop
    print("")
    return fractional_occupancy

def calculate_free_ligand_concentration(ref,trace):
    #read and append trace list
    #begin loop
      #for each L(total)
      #calculate: L(free)= L(total)-fR(1)*R(total)-fR(2)*R(total)-...-fR(n)*R(total) bzw. L(free)= L(total)-R(total)*(fR(1)+fR(2)+...+fR(n))
    sum_fractional_occupancies = 0
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
        #for each footprinted peak
        if trace_peak.fractional_occupancy > 0:
            sum_fractional_occupancies += trace_peak.fractional_occupancy
        
    ## DEBUG
    #print(sum_fractional_occupancies)        
        
    Lfree_conc = trace.Ltot_conc - trace.Rtot_conc*sum_fractional_occupancies
    #add result to trace list (L(free))
    #end loop
    print("")
    return Lfree_conc

def fitFunc(Lfree_conc, Kd):
    #for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace):
        #Func = (Lfree_conc)/(Lfree_conc * Kd)
    #Kd = 
    return (Lfree_conc)/(Lfree_conc + Kd)

def fit_data_determine_kd(ref, trace_list, Lfree_conc, Kd, fractional_occupancy):
    #read an append trace list
    #begin loop
    #fitFunc = (Lfree_conc)/(Lfree_conc + Kd)
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace_list):
           if trace_peak.footprinted_peak == 1:
              Lfree = numpy.array([Lfree_conc])
              fR = numpy.array([fractional_occupancy]) 
              Kd_values, Covar = curve_fit(fitFunc, Lfree, fR, Kd)
      #for each footprinting site
      #fit fR(n)=L(free)/(Kd(n)+L(free))
      #add result to trace list
    #end loop
    print("")
    return Kd_values, Covar

def plot_data(ref, trace_list, Lfree_conc, fractional_occupancy, kd_values, Covar):
    for ref_peak,trace_peak in give_all_clustered_peaks(ref,trace_list):
           #if trace_peak.footprinted_peak == 1:
              plt.ylabel('Fractional Occupancy', fontsize = 16)
              plt.xlabel('Free Ligand Concentration', fontsize = 16)
              plt.title(trace_peak.size_bp)
              #plt.text(60, .025, kd_values)
              #datapoints and errorbars
              plt.errorbar(Lfree_conc, fractional_occupancy, fmt = 'ro', yerr = 0.2)
              #sigma = [Covar[0,0], \
              #Covar[1,1], \
              #Covar[2,2] \
              # ]
              Lfree = numpy.array([Lfree_conc])
              fR = numpy.array([fractional_occupancy])
              Kd = kd_values
              plt.plot(Lfree, fR,\
              Lfree, fitFunc(Lfree_conc, Kd))
              #t, fitFunc(t, Kd[0] + sigma[0], Kd[1] - sigma[1], Kd[2] + sigma[2]),\
              #t, fitFunc(t, Kd[0] - sigma[0], Kd[1] + sigma[1], Kd[2] - sigma[2]) )
    #plot fR vs L(free)
    #plot fit
    #show Kd and Basepair number
    plt.show()          
    print("")
    return 

#def write_data_to_csv ()

  #return


if __name__ == "__main__":
    main(sys.argv[1:])
