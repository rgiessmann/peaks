#!/usr/bin/env python

import logging
import sys
import getopt
import numpy
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import copy
import pandas
import scipy.optimize

## turn on logging
logging.basicConfig(level=logging.CRITICAL)
global log
log = logging   

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
        log.fatal('You provided unusual arguments. Call me with -h to learn more.')
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
        log.error("You provided too many options. Call me with -h to learn more.")

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

class Peak:
    def __init__(self, size_bp, peak_height):
        self.size_bp = float(size_bp)
        self.peak_height = float(peak_height)
        return
    def __repr__(self):
        return repr(vars(self))

class Index:
    pass    

def list_traces(read_filelist="../HexA.csv"):
    """
    Produce a list of all traces found in read_filelist (which can be a string\
    or a list of strings). Saves this list of traces to file "output_traces.csv".
    """
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

def get_data(read_filelist="../HexA.csv", config_file="../input_traces.csv"):
    """
    Reads data from read_filelist, and returns an object containing all read
    information.
    
    This object consists of the following structure:
    
    [Trace(..., peaks = [Peak(...), Peak(...)], Trace(...)]
    
    See "Trace" and "Peak" to learn more about these classes. 
    """    

    import csv
    trace_list = []
    
    ## TODO: what entries to accept?
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
            header = next(csv_reader)
            index = Index()
            index.peak_height = header.index("Height")
            index.size_bp = header.index("Size")
            index.file_name = header.index('Sample File Name')
            index.sample_name = header.index('Sample Name')
            print(header)
            
            for row in csv_reader:
                #rules for entry acceptance?

                if "Dye/Sample Peak" in header: #split_combined_fields == 
                    row.extend(row[header.index('Dye/Sample Peak')].split(","))
                    index.dye = len(header)
                    index.sample_peak = len(header)+1
                else:
                    index.dye = header.index("Dye")
                    index.sample_peak = header.index("Sample Peak")
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
    
    return trace_list

def generate_averaged_negative_control(trace_list,accepted_offset=0.5, factor_method="peak", *args, **kwargs):
    
    ## clean up trace_list, in case it was processed already.
    trace_list_ref = copy.copy(trace_list)
    for t in trace_list_ref:
        for p in t.peaks:
            if "cluster" in vars(p):
                del(p.cluster)
                pass
    
    conc_0_traces = []
    for t in trace_list_ref:
        if t.Ltot_conc == 0:
            conc_0_traces.append(t)
    print("Found {} conc_0 traces.".format(len(conc_0_traces)))        
            
    if len(conc_0_traces) > 1:
        
        i = 0       
                                
        for trace in conc_0_traces:
            #print("Screening trace {}...".format(trace.file_name))
            for ref_peak in trace.peaks:                
                #print("+ Screening peak {}:".format(ref_peak))
                if "cluster" not in vars(ref_peak):
                    i += 1
                    ref_peak.cluster = i
                    for t in conc_0_traces:
                        #print("++ Screening trace {}...".format(t.file_name))
                        for peak in t.peaks:
                            #print("+++ Screening peak {}:".format(peak))
                            if abs(peak.size_bp - ref_peak.size_bp) < accepted_offset:
                                #print("++++ Found a peak closeby!")
                                if "cluster" not in vars(peak):
                                    #print("+++++ Identified a peak to cluster with.")
                                    peak.cluster = i
                                else:
                                    del(ref_peak.cluster)
                                    if peak != ref_peak:
                                        print("+++++ A peak is lying within the accepted offset for TWO clusters, not just one, as expected! This should not happen!")
                                        ref_peak.cluster = i
                                        print("{} : {}".format(trace.file_name, ref_peak))
                                        print("{} : {}".format(t.file_name, peak))
                                    else:
                                        #print("+++++ I found myself again...")
                                        ref_peak.cluster = i
        
        num_total_clusters = i
        print("Found {} potential clusters in total.".format(num_total_clusters))
        
        num_total_traces = float(len(conc_0_traces))

        
        ## gradually fit optimal heights
        for x in range(5):
            print("--- ROUND {} ---".format(x))    
            #print(conc_0_traces)
            ref = Trace("averaged_negative_control", "B", 0, 0, [])

            ## 1. create dummy peaks
            for i in range(1,num_total_clusters+1):
                #print(" ")
                print("Evaluate cluster {:3}".format(i)) 
                ## fill peak storage
                trace_storage = []
                for trace in conc_0_traces:
                    for peak in trace.peaks:
                        if peak.cluster == i:
                            trace_storage.append(peak)
                print(trace_storage)

                averaged_peak_height_mean = numpy.mean([p.peak_height for p in trace_storage])
                averaged_peak_size_bp_mean = numpy.mean([p.size_bp for p in trace_storage])

                ref.peaks.append(Peak(averaged_peak_size_bp_mean,averaged_peak_height_mean))


            #print(len(ref.peaks))
            ## 2. fit optimal factors
            cluster_peaks(ref, conc_0_traces, accepted_offset=accepted_offset)

            #print(conc_0_traces)
            
            if factor_method == "num":
                optimal_factors = determine_factor_numerically(ref, conc_0_traces, *args, **kwargs)
            elif factor_method == "peak":
                optimal_factors = determine_factor_single_peak(ref, conc_0_traces, *args, **kwargs)

            ## 3. correct with factor
            for index, trace in enumerate(conc_0_traces):
                correct_peaks_with_factor(trace,optimal_factors[index])
                #print index
                #print trace
        
        ref = Trace("averaged_negative_control", "B", 0, 0, [])        
        num_conc_0_traces = float(len(conc_0_traces))
        for i in range(1,num_total_clusters+1):
            trace_storage = []
            for trace in conc_0_traces:
                for peak in trace.peaks:
                    if peak.cluster == i:
                        trace_storage.append(peak)
            averaged_peak_height_mean = numpy.mean([p.peak_height for p in trace_storage])
            averaged_peak_height_n = len(trace_storage)
            averaged_peak_height_sd = numpy.std([p.peak_height for p in trace_storage])
            averaged_peak_size_bp_mean = numpy.mean([p.size_bp for p in trace_storage])
            averaged_peak_size_bp_n = len(trace_storage)
            averaged_peak_size_bp_sd = numpy.std([p.size_bp for p in trace_storage])
            
            ref.peaks.append(Peak(averaged_peak_size_bp_mean,averaged_peak_height_mean))
            ref.peaks[i-1].averaged_peak_height_n = averaged_peak_height_n
            ref.peaks[i-1].averaged_peak_height_nm = averaged_peak_height_n / num_conc_0_traces
            ref.peaks[i-1].averaged_peak_height_sd = averaged_peak_height_sd
            ref.peaks[i-1].averaged_peak_size_bp_n = averaged_peak_size_bp_n
            ref.peaks[i-1].averaged_peak_size_bp_sd = averaged_peak_size_bp_sd
            ref.peaks[i-1].averaged_peak_size_bp_nm = averaged_peak_size_bp_n / num_conc_0_traces
                


        

                            
            
    elif len(conc_0_traces) == 1:
        print("You called generate_averaged_negative_control although only 1 trace with Ltot = 0 is available. Returning the clustered trace.")        
        ## create new object
        c = conc_0_traces[0]    
        ref = Trace(c.file_name, c.dye_color, c.Ltot_conc, c.Rtot_conc)
        for p in c.peaks:
            ref.peaks.append(Peak(p.size_bp, p.peak_height))

        i = 0                    
        for ref_peak in ref.peaks:
            i += 1
            ref_peak.cluster = i

    else:
        print("The provided trace_list contains no negative control traces with Ltot = 0.")

    ## TODO: is this necessary?
    ## clean-up        
    del(trace_list_ref)
    
    ref.peaks = sorted(ref.peaks, key=lambda k: k.size_bp)
                
    return ref


def plot_height_size_overview_averaged_negative_control(ref, *args, **kwargs):
    fig, ax = plt.subplots(1, figsize=(15,10))
    plt.title("Reference trace - Overview Plot")
    plt.ylabel('averaged peak height (A.U.)', fontsize = 16)
    plt.xlabel('averaged peak size (bp)', fontsize = 16)
    if "ylim" in kwargs:
        plt.ylim(kwargs["ylim"])
    if "xlim" in kwargs:
        plt.xlim(kwargs["xlim"])    
   
    ## plot points
    heights = []
    heights_err = []
    
    sizes = []
    sizes_err = []
        
    for ref_peak in ref.peaks:
        height = ref_peak.peak_height
        height_err = ref_peak.averaged_peak_height_sd
        
        size = ref_peak.size_bp
        size_err = ref_peak.averaged_peak_size_bp_sd
        
        heights.append(height)
        heights_err.append(height_err)
        
        sizes.append(size)
        sizes_err.append(size_err)
        
    ax.errorbar(x=sizes, y=heights, xerr=sizes_err, yerr=heights_err, elinewidth=1, linewidth=0, marker='.')

    plt.show()
    
    return


def plot_height_size_clustering_success_averaged_negative_control(ref, *args, **kwargs):

    fig, axarr = plt.subplots(2, figsize=(15,10), sharex=True)
    axarr[0].set_title("Reference trace - Clustering Success Plot")
    axarr[0].set_ylabel('peak height clustering success ratio (n/N)')
    axarr[1].set_ylabel('peak size clustering success ratio (n/N)')
    axarr[1].set_xlabel('averaged peak size (bp)')
    fig.subplots_adjust(hspace=0)
    if "ylim" in kwargs:
        plt.ylim(kwargs["ylim"])
    if "xlim" in kwargs:
        plt.xlim(kwargs["xlim"])  
    
    ## plot points
    sizes = []
    
    heights_repr = []
    sizes_repr = []
        
    for ref_peak in ref.peaks:
        size = ref_peak.size_bp
        
        sizes.append(size)
        
        height_repr = ref_peak.averaged_peak_height_nm
        size_repr = ref_peak.averaged_peak_size_bp_nm
        
        heights_repr.append(height_repr)
        sizes_repr.append(size_repr)

    from matplotlib import cm
    
    _colors = [(a,a,a) for a in sizes_repr]
    axarr[0].scatter(x=sizes, y=heights_repr, marker='o')
    axarr[1].scatter(x=sizes, y=sizes_repr, marker='o')
    
    plt.show()
    
    return


def cluster_peaks(ref,trace_list,accepted_offset=0.25):
    """
    Compiles peaks from trace_list to clusters of peaks which match to one 
    particular peak in ref. This function directly modifies the objects given
    as arguments!
    
    accepted_offset is giving the allowed difference in size_bp of peaks to be
    matched.
    
    Peaks that don't match any peak from ref are numbered cluster 0.
    
    """
    
    i = 0            

    for ref_peak in ref.peaks:
        i += 1
        ref_peak.cluster = i
        for trace in trace_list:
            for peak in trace.peaks:
                if abs(peak.size_bp - ref_peak.size_bp) < accepted_offset:
                    peak.cluster = i
                    
    for trace in trace_list:
        for peak in trace.peaks:
            if "cluster" not in vars(peak):
                peak.cluster = 0
                
    return 


def calculate_deviance_for_all_peaks(ref, trace_list, weight_smaller=1,weight_bigger=1, weight_by_inverse_height=False, from_bp=20, to_bp=170):
    '''
    Calculates the area between peaks that were identified as clustered in ref, 
    compared to trace. Allows for comparison of two traces only.
    
    
    Options weight_smaller, weight_bigger are used to consider only peaks from 
    trace which are smaller or bigger, respectively, than the reference peak.
    
    Option //BROKE relative_mode allows to calculate the deviation in relative terms to
    the reference peak.
    
    ## TODO: 
    from_bp .. to_bp 
    
    ## TODO: calculate deviance for trace_list --> saves computing time
    '''

    warning_not_all_peaks_match = False
    warning_trace_range = False

    deviance_for_smaller_peaks = 0
    deviance_for_bigger_peaks = 0        
    num_smaller=0
    num_bigger=0
    
    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace_list):
        ## WORKAROUND for single trace mode
        # if there are no peaks clustered to the ref_peak, they cannot be included --> continue with next pair
        if trace_peaks == []:
            if warning_not_all_peaks_match == False:
                print("INFO: Not all peaks match -- omitting deviation for non-comparable peaks.")
                warning_not_all_peaks_match = True
            continue
        
        if ref_peak.size_bp < from_bp or ref_peak.size_bp > to_bp:
            if warning_trace_range == False:
                print("INFO: Trace contains peaks which are not included in deviation calculation.")
                warning_trace_range = True

        ## allows to calculate deviance for one trace only
        trace_peak=trace_peaks[0]
        
        if from_bp < ref_peak.size_bp < to_bp:
            if weight_by_inverse_height == True: 
                ## this mode calculates deviation as percentage point, with different
                ## weights for each peak
                weight = 1/ref_peak.peak_height
            else:
                ## similarly weighted, the difference between the two trace functions
                ## ~ area / integral between traces is calculated
                weight = 1
            
            if trace_peak.peak_height <= ref_peak.peak_height:
                # deviation for smaller, i.e. potentially footprinted peaks
                deviance_for_smaller_peaks += abs((ref_peak.peak_height - trace_peak.peak_height)*weight)
                num_smaller +=1
            else:
                # deviation for bigger, i.e. potentially hypersensitive peaks
                deviance_for_bigger_peaks += abs((ref_peak.peak_height - trace_peak.peak_height)*weight)
                num_bigger +=1

    weighted_deviation = (weight_smaller*deviance_for_smaller_peaks + weight_bigger*deviance_for_bigger_peaks)

    ## possible further outputs:
    ## num_smaller
    ## num_bigger
    ## -> mean deviation per peak ~ normalizing

    return weighted_deviation
    
    
    
def give_all_clustered_peaks(ref,trace_list):
    """
    Generates lists of clustered peaks from traces in trace_list, clustered to 
    the peaks in ref.
    
    Yields (ref_peak, [trace_peak, trace_peak, ...]) elements.
    """

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
            #if len(trace_peaks) != len(trace_list):
                #print("WARNING: Unable to find clustered peaks for all given traces...")
            yield (ref_peak[0],trace_peaks)


def determine_factor_numerically(ref, trace_list, weight_smaller=1, weight_bigger=1, relative_mode=False, from_bp=20, to_bp=170):
    """
    Determines the optimal factor for trace, when compared to ref. This function
    minimizes the deviation as calculated by calculate_deviance_for_all_peaks()
    with the given options.
    
    Returns the optimal factor.
    ## TODO: implement real optimizer via scipy.optimize()    
    """
    
    optimal_factors = [1 for trace in trace_list]
    rmsd_old = calculate_deviance_for_all_peaks(ref, trace_list, weight_smaller, weight_bigger, relative_mode, from_bp, to_bp)
    
    print("starting: no calibration --> factors {!s} ; deviation {:8f}".format(optimal_factors, rmsd_old))
    
    for index, trace in enumerate(trace_list):
        
        ## store original peak heights
        for peak in trace.peaks:
            if "peak_height_original" not in vars(peak):
                peak.peak_height_original = peak.peak_height
            else:
                ## was already modified by factor
                peak.peak_height_original = peak.peak_height_original 
                pass
                
    
        def cost_function(x):
            _factor = x
            
            correct_peaks_with_factor(trace,_factor)
            rmsd_new = calculate_deviance_for_all_peaks(ref,trace,weight_smaller,weight_bigger,relative_mode, from_bp, to_bp)
            return rmsd_new
        
        x0 = [0.]
        min_result = scipy.optimize.minimize(cost_function, x0, method='nelder-mead') #, options={'disp': True})
        optimal_factor = min_result.x[0]
        
        print("found optimal: deviation {:8f} --> factor {:4.2f} ".format(cost_function(optimal_factor), optimal_factor))

        optimal_factors[index] = optimal_factor
                
    ##
    for index, trace in enumerate(trace_list):
        correct_peaks_with_factor(trace_list[index],optimal_factors[index])
        print("calculating: trace # {:2} --> factor {:4.2f}".format(index, optimal_factors[index]))

    ## use calculate_deviance_for_all_peaks with trace and ref
    rmsd_new = calculate_deviance_for_all_peaks(ref,trace_list, weight_smaller, weight_bigger, relative_mode, from_bp, to_bp)

    print("found optimal: deviation {:8f} --> factors {!s} ".format(rmsd_new, optimal_factors))
    
    
    ## restore all peak heights to original height
    for trace in trace_list:
        for peak in trace.peaks:
            peak.peak_height = peak.peak_height_original 

        
    return optimal_factors


def determine_factor_single_peak(ref, trace_list, weight_smaller=1, weight_bigger=1, relative_mode=False, from_bp=20, to_bp=170):
    """
    Determines the optimal factor for trace, when compared to ref. This function
    minimizes the deviation as calculated by calculate_deviance_for_all_peaks()
    with the given options. However, only factors corresponding to existent
    peak pairs are considered.
    
    Returns the optimal factor.

    """
   
    optimal_factors = [1 for trace in trace_list]
    rmsd_old = calculate_deviance_for_all_peaks(ref, trace_list, weight_smaller, weight_bigger, relative_mode, from_bp, to_bp)
    which_peak = None
    
    print("starting: no calibration --> factors {!s} ; deviation {:8f}".format(optimal_factors, rmsd_old))


    ## store original peak heights
    for trace in trace_list:
        for peak in trace.peaks:
            if "peak_height_original" not in vars(peak):
                peak.peak_height_original = peak.peak_height
            else:
                ## was already modified by factor
                peak.peak_height_original = peak.peak_height_original 
                pass
                
        
    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace_list):

        ## WORKAROUND for single trace mode
        if trace_peaks == []:
            continue
        
        ## if not all peaks are clustered on this ref_peak, they cannot be used --> continue with next cycle
        if len(trace_peaks) != len(trace_list):
            continue

            
        testing_factors = []
        
        ## works with one trace at a time
        for index, trace_peak in enumerate(trace_peaks):
            ## have to rely on trace_peaks coming in same order as trace_list ...
            factor =  ref_peak.peak_height / trace_peak.peak_height_original
            correct_peaks_with_factor(trace_list[index],factor)
            testing_factors.append(factor)
            ## DEBUG
            print("testing: trace # {:2} on cluster {:3} @ {:8.1f} bp --> factor {:4.2f}".format(index, trace_peak.cluster, trace_peak.size_bp, factor))

        ## use calculate_deviance_for_all_peaks with trace and ref
        rmsd_new = calculate_deviance_for_all_peaks(ref,trace_list, weight_smaller, weight_bigger, relative_mode, from_bp, to_bp)

        

        ## DEBUG
        print("was testing: all traces on cluster {:3} @ {:8.1f} bp --> deviation {:8f}".format(ref_peak.cluster, ref_peak.size_bp, rmsd_new))
                    
        ## compare deviance_new with deviance_old, if better deviance_new -> deviance_old, else delete deviance_new
        if rmsd_new <= rmsd_old:
            rmsd_old = rmsd_new
            optimal_factors = testing_factors
            which_peak = ref_peak
        
    if which_peak == None:
        print("no peak would be better than current factorization --> factors {!s} ; deviation {:8f}".format( optimal_factors, rmsd_old))
    else:
         ## DEBUG        
        print("found optimal: cluster {:3} @ {:8.1f} bp ; deviation {:8f} --> factors {!s} ".format(which_peak.cluster, which_peak.size_bp, rmsd_old, optimal_factors))
        
    ## restore all peak heights to original height
    for trace in trace_list:
        for peak in trace.peaks:
            peak.peak_height = peak.peak_height_original 

    return optimal_factors 


def correct_peaks_with_factor(trace, factor):
    """
    Corrects all peak_heights in trace with factor. The multiplication is applied
    directly to the object given as argument.
    """

    for peak in trace.peaks:
        peak.peak_height = peak.peak_height_original * factor

    return


def mark_footprinted_peaks(ref, trace_list, threshold=0.1, mark_all=False):
    """
    Marks potentially footprinted peaks, i.e. peaks to be evaluated, by setting
    the flag footprinted_peak, if the peak_heights of ref and trace differ by
    threshold*ref_peak.peak_height (standard value = 10%) and if the ligand 
    concentrations of the compared traces differ from each other.
    
    Sets the flag to True, if the conditions above are fulfilled, otherwise
    to False.
    
    Modifies the object trace directly.
    """

    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace_list):

    ##DEBUG
    #print(ref_peak,trace_peaks)
    
    ## WORKAROUND for single trace mode
    # if there are no peaks clustered to the ref_peak, they cannot be marked --> continue with next cycle
        if trace_peaks == []:
            continue
        
        for trace_peak in trace_peaks:
            if ref_peak.peak_height - trace_peak.peak_height > threshold*ref_peak.peak_height:
                trace = [t for t in trace_list if trace_peak in t.peaks]
                if len(trace) == 1:
                    trace = trace[0]
                else:
                    logging.critical("!")
                    raise Exception
                if ref.Ltot_conc != trace.Ltot_conc:
                    trace_peak.footprinted_peak = True
                else:
                    trace_peak.footprinted_peak = False                
            else:
                trace_peak.footprinted_peak = False        
                
            if mark_all == True:
                trace_peak.footprinted_peak = True
            
    return 

def add_fractional_occupancies(ref,trace_list):
    """
    TODO:    ...
    
    Modifies the object trace directly.
    """

    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace_list):
        for trace_peak in trace_peaks:
            trace_peak.fractional_occupancy = 1 - trace_peak.peak_height / ref_peak.peak_height        

    return 


def calculate_free_ligand_concentration(ref,trace,mode="Ltot"):
    """
    TODO:    ...
    
    Calculates the free ligand concentration.
    
    Returns Lfree_conc.
    """
    
    ## // BROKEN
    #for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace):
        ## WORKAROUND for single trace mode
        # if there are no peaks clustered to the ref_peak, they don't need to be considered --> continue with next cycle
        #if trace_peaks == []:
        #    #print("Skipping in calc_Lfree_conc")
        #    continue
        ## works for one trace only
        #trace_peak = trace_peaks[0]
        ## DEBUG
        #print(trace_peak)
        ## DEBUG
        #trace_peak.Lfree_conc = trace.Ltot_conc - trace.Rtot_conc*trace_peak.fractional_occupancy
    
    if mode=="Ltot":
    ## to make this reproducible with old evaluations:    
        Lfree_conc = trace.Ltot_conc 
    else:
        ## no other modes than "Ltot" available!
        raise
    
    return Lfree_conc


def fitFunc_fR(Lfree_conc, Kd):
    ## underlying formula:
    ## y = fR
    ## x = Lfree_conc
    ## to be fitted: Kd
    
    fR = (Lfree_conc)/(Lfree_conc + Kd)
    
    return fR


def fit_data_determine_kd(ref, trace_list):
    """
    Determines the K_D by fitting of observed fractional occupancies against
    calculated free ligand concentration.
    
    """
    
    ## data has to be in the following format:
    ## xdata = [Lfree_conc1, Lfree_conc2, ...]
    ## ydata = [fractional_occupancy1, fractional_occupancy2, ...]    
    
    ## output:
    ## [ K_D(cluster1), K_D(cluster2), ... ]

    KD_matrix = []

    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace_list):

        xdata = []
        ydata = []

        ## DEBUG
        #print(ref_peak,trace_peak)            
        

        for trace_peak in trace_peaks:



            if trace_peak.footprinted_peak == True:
                xdata.append(calculate_free_ligand_concentration(ref, [trace for trace in trace_list if trace_peak in trace.peaks][0]))
                ydata.append(trace_peak.fractional_occupancy)
                
            ## TODO: shall we catch out clusters with no footprinted peaks at all?
            ## ...

        ## DEBUG
        #print(xdata,ydata)

                
        ## fitting...
        popt, pcov= curve_fit(fitFunc_fR, xdata, ydata)

        # compute SE, i.e. standard deviation errors 
        perr = numpy.sqrt(numpy.diag(pcov))


        ## save results...
        ## TODO: optimize this form.
        #KD_matrix.append(["cluster "+str(ref_peak.cluster), popt[0], perr[0]])
        KD_matrix.append(["cluster "+str(ref_peak.cluster), popt[0], perr[0], len(ydata), float(len(ydata))/float(len(trace_list))])
        print(KD_matrix[-1])


    return KD_matrix


def generate_xdata_ydata(ref,trace_list,cluster):
    """
    
    """    
    
    xdata = []
    ydata = []    
    
    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace_list):

        ## assumes that there is only one ref_peak with correct cluster number
        if ref_peak.cluster == cluster:
            for trace_peak in trace_peaks:        
                ## DEBUG
                #print(ref_peak,trace_peak)            
                        
                if trace_peak.footprinted_peak == True:
                    xdata.append(calculate_free_ligand_concentration(ref, [trace for trace in trace_list if trace_peak in trace.peaks][0]))
                    ydata.append(trace_peak.fractional_occupancy)
                    
            break
    return xdata, ydata
            

def plot_data(ref, trace_list, cluster, kd_matrix):

    fig, ax = plt.subplots(1)
    
    plt.title("Cluster " + str(cluster))
    plt.ylabel('Fractional Occupancy', fontsize = 16)
    plt.xlabel('Free Ligand Concentration', fontsize = 16)
    plt.ylim([-0.1, 1.1])
    plt.xlim([-1, 10.1])

    ## plot points
    xdata, ydata = generate_xdata_ydata(ref,trace_list,cluster)                
    ax.plot(xdata,ydata, 'o')


    ## plot line
    kd = [entry[1] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)]

    import numpy

    x = numpy.linspace(0,15,1000) # linearly spaced numbers
    y = fitFunc_fR(x,kd)
    ax.plot(x,y)
    
    ## show plot
    plt
    
    bp  = float([ref_peak.size_bp for ref_peak in ref.peaks if ref_peak.cluster == cluster][0])
    n   = int([entry[3] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)][0])
    n_m = float([entry[4] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)][0])
    kd  = float([entry[1] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)][0])
    std = float([entry[2] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)][0])
    
    textstr = """\
    pos = {:4.1f} bp
    n   = {}
    n/m = {:3.2}
    $K_D$ = {:3.2f} $\pm$ {:3.2f} $\mu$M
    """.format(bp, n, n_m, kd, std)


    # place a text box in upper left in axes coords
    props = dict(facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, family="monospace",
        verticalalignment='top', bbox=props)
    
    plt.show()
    
    return


def save_kd(kd_matrix, filename="kd-matrix.csv"):
    """
    Outputs the KD matrix in a comma-separated table (csv file). The
filename is
    as given (no cross-check for pre-existing files of the same name!), and
the
    KD matrix which is given, is written to the file.

    The columns are: "Cluster #" (which cluster is evaluated), "KD mean"
(the fitted
    KD value), "KD SD" (the estimated standard deviation of the fitted KD
value),
    "n" (number of data points used), "n/m" (the ratio of the number of
data points
    used to the number of all traces).
    """
    import csv
    with open(filename,"w") as f:
        w=csv.writer(f)
        w.writerow(["Cluster #", "KD mean", "KD SD", "n", "m/n"])
        for row in kd_matrix:
            w.writerow(row)
    return


def plot_peakscan(ref, trace_list, *args, **kwargs):

    _columns=["trace", "cluster","peak_height"]

    df = pandas.DataFrame(data=None, columns=_columns)
    
    iterlist = copy.copy(trace_list)
    iterlist.append(ref)
    
    for trace in iterlist:

        ddict = {"trace" : [trace.file_name for peak in trace.peaks], "cluster" : [float(peak.cluster) for peak in trace.peaks], "peak_height" : [peak.peak_height for peak in trace.peaks] }
        row = pandas.DataFrame.from_dict(data=ddict)

        zerorow1 = row.copy()
        zerorow1["cluster"] = zerorow1["cluster"] - 0.5
        zerorow1["peak_height"] = 0.0

        zerorow2 = row.copy()
        zerorow2["cluster"] = zerorow2["cluster"] + 0.5
        zerorow2["peak_height"] = 0.0

        row = row.append(zerorow1, ignore_index=True)
        row = row.append(zerorow2, ignore_index=True)

        row = row.sort_values(by="cluster")

        df = df.append(row, ignore_index=True)

        #print row
        ## WORKAROUND
        #break

    #df.head(n=15).plot.line(x="cluster", )
    #df.plot.line()

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(13,15))

    for label, df2 in df.groupby('trace'):
        if label=="averaged_negative_control":
            df2.plot(x="cluster", y="peak_height", kind="line", ax=ax, label=label, marker="s", linewidth=2)
        else:
            #print label
            df2.plot(x="cluster", y="peak_height", kind="line", ax=ax, label=label)
        
        #print df2.head(n=15)
        #break


    #ax.set_xlim([0, 5])
    leg = plt.legend(loc="upper right", ncol=2)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)

    if "ylim" in kwargs:
        plt.ylim(kwargs["ylim"])
    if "xlim" in kwargs:
        plt.xlim(kwargs["xlim"])    
        
    plt.show()
    
    return



if __name__ == "__main__":
    main(sys.argv[1:])
