#!/usr/bin/env python

import logging
import sys
import getopt
import numpy
import scipy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


## turn on logging
logging.basicConfig(level=logging.INFO)
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
            header = csv_reader.__next__()
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
                if ref_peak.size_bp - accepted_offset < peak.size_bp < ref_peak.size_bp + accepted_offset:
                    peak.cluster = i
                    
    for trace in trace_list:
        for peak in trace.peaks:
            if "cluster" not in vars(peak):
                peak.cluster = 0
                
    return 


def calculate_deviance_for_all_peaks(ref, trace, weight_smaller=1,weight_bigger=1, weight_by_inverse_height=False, from_bp=20, to_bp=125):
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

    deviance_for_smaller_peaks = 0
    deviance_for_bigger_peaks = 0        
    num_smaller=0
    num_bigger=0
    
    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace):
        ## WORKAROUND for single trace mode
        # if there are no peaks clustered to the ref_peak, they cannot be included --> continue with next pair
        if trace_peaks == []:
            print("WARNING: Not all peaks match -- omitting deviation for non-comparable peaks.")
            continue
        
        if from_bp < ref_peak.size_bp < to_bp:
            print("INFO: Trace contains peaks which are not included in deviation calculation.")

        ## allows to calculate deviance for one trace only
        trace_peak=trace_peaks[0]
        

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
            if len(trace_peaks) != len(trace_list):
                print("WARNING: Unable to find clustered peaks for all given traces...")
            yield (ref_peak[0],trace_peaks)


def determine_factor_numerically(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False):
    """
    Determines the optimal factor for trace, when compared to ref. This function
    minimizes the deviation as calculated by calculate_deviance_for_all_peaks()
    with the given options.
    
    Returns the optimal factor.

    ## TODO: implement real optimizer via scipy.optimize()    
    """
    
    optimal_factor = 1
    rmsd_old = calculate_deviance_for_all_peaks(ref,trace,weight_smaller,weight_bigger, relative_mode)

    ## store original peak heights
    for peak in trace.peaks:
        peak.peak_height_original = peak.peak_height

    for factor in numpy.arange(0,3.5,0.01):

        correct_peaks_with_factor(trace,factor)
            
        rmsd_new = calculate_deviance_for_all_peaks(ref,trace,weight_smaller,weight_bigger,relative_mode)
                
        ## restore all peak heights to original height
        for peak in trace.peaks:
            peak.peak_height = peak.peak_height_original 
        
        ## DEBUG
        #print(str(factor)+" : "+str(rmsd_new))

        #compare deviance_new with deviance_old, if better deviance_new -> deviance_old, else delete deviance_new
        if rmsd_new < rmsd_old:
            rmsd_old = rmsd_new
            optimal_factor = factor
            
    ## in the end: delete stored original peak heights
    for peak in trace.peaks:
        del(peak.peak_height_original)

    return optimal_factor


def determine_factor_single_peak(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False):
    """
    Determines the optimal factor for trace, when compared to ref. This function
    minimizes the deviation as calculated by calculate_deviance_for_all_peaks()
    with the given options. However, only factors corresponding to existent
    peak pairs are considered.
    
    Returns the optimal factor.

    """
   
    optimal_factor = 1
    rmsd_old = calculate_deviance_for_all_peaks(ref, trace, weight_smaller, weight_bigger, relative_mode)

    ## store original peak heights
    for peak in trace.peaks:
        peak.peak_height_original = peak.peak_height
        
    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace):

        ## WORKAROUND for single trace mode
        # if there are no peaks clustered to the ref_peak, they cannot be used --> continue with next cycle
        if trace_peaks == []:
            continue

        ## works with one trace only 
        trace_peak = trace_peaks[0]
        
        factor =  ref_peak.peak_height / trace_peak.peak_height 
        correct_peaks_with_factor(trace,factor)
            
        #use calculate_deviance_for_all_peaks with trace and ref
        rmsd_new = calculate_deviance_for_all_peaks(ref,trace, weight_smaller, weight_bigger, relative_mode)
     
        ## restore all peak heights to original height
        for peak in trace.peaks:
            peak.peak_height = peak.peak_height_original 
        
        ## DEBUG        
        #print(str(ref_peak.cluster) + " " + str(ref_peak.size_bp) + " -- " + str(factor)+" : "+str(rmsd_new))
        
        #compare deviance_new with deviance_old, if better deviance_new -> deviance_old, else delete deviance_new
        if rmsd_new < rmsd_old:
            rmsd_old = rmsd_new
            optimal_factor = factor
            
    return optimal_factor


def correct_peaks_with_factor(trace, factor):
    """
    Corrects all peak_heights in trace with factor. The multiplication is applied
    directly to the object given as argument.
    """
        
    for peak in trace.peaks:
        peak.peak_height = peak.peak_height * factor

    return


def mark_footprinted_peaks(ref, trace, threshold=0.1):
    """
    Marks potentially footprinted peaks, i.e. peaks to be evaluated, by setting
    the flag footprinted_peak, if the peak_heights of ref and trace differ by
    threshold*ref_peak.peak_height (standard value = 10%) and if the ligand 
    concentrations of the compared traces differ from each other.
    
    Sets the flag to True, if the conditions above are fulfilled, otherwise
    to False.
    
    Modifies the object trace directly.
    
    ## TODO: implement trace_list -> saves computing time
    """

    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace):

        ##DEBUG
        #print(ref_peak,trace_peaks)
        
        ## WORKAROUND for single trace mode
        # if there are no peaks clustered to the ref_peak, they cannot be marked --> continue with next cycle
        if trace_peaks == []:
            continue
        
        ## works for one trace only
        trace_peak = trace_peaks[0]
        
        if ref_peak.peak_height - trace_peak.peak_height > threshold*ref_peak.peak_height and ref.Ltot_conc != trace.Ltot_conc:
            trace_peak.footprinted_peak = True
        else:
            trace_peak.footprinted_peak = False        
            
    return 

def add_fractional_occupancies(ref,trace):
    """
    TODO:    ...
    
    Modifies the object trace directly.
    """

    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace):
        for trace_peak in trace_peaks:
            trace_peak.fractional_occupancy = 1 - trace_peak.peak_height / ref_peak.peak_height        

    return


def calculate_free_ligand_concentration(ref,trace):
    """
    TODO:    ...
    
    Calculates the free ligand concentration when considering all footprinted 
    peaks to diminish the pool of free ligand.
    
    Returns Lfree_conc.
    """

    sum_fractional_occupancies = 0
    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace):

        ## WORKAROUND for single trace mode
        # if there are no peaks clustered to the ref_peak, they don't need to be considered --> continue with next cycle
        if trace_peaks == []:
            #print("Skipping in calc_Lfree_conc")
            continue
        
        ## works for one trace only
        trace_peak = trace_peaks[0]

        ## DEBUG
        #print(trace_peak)
        
        if trace_peak.footprinted_peak:
            sum_fractional_occupancies += trace_peak.fractional_occupancy
    
    ## DEBUG
    #print(sum_fractional_occupancies)        
        
    Lfree_conc = trace.Ltot_conc - trace.Rtot_conc*sum_fractional_occupancies
    
    ## to make this reproducible with old evaluations:    
    #Lfree_conc = trace.Ltot_conc 
    
    return Lfree_conc


def fitFunc_fR(Lfree_conc, Kd):
    ## underlying formula:
    ## y = fR
    ## x = Lfree_conc
    ## to be fitted: Kd
    
    ## fR = (Lfree_conc)/(Lfree_conc + Kd)
    return (Lfree_conc)/(Lfree_conc + Kd)


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

        ## TODO: can we truly assume this all the time?!
        xdata.append(0) #ref.Ltot_conc)
        ydata.append(0) #ref_peak.fractional_occupancy)
        ## alternative:
        # xdata.append(ref.Ltot_conc)
        # ydata.append(ref_peak.fractional_occupancy)
        ## -> doesn't work, because fractional_occupancies are not set for ref, so far...

        for trace_peak in trace_peaks:

            ## DEBUG
            #print(ref_peak,trace_peak)            
            
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
        KD_matrix.append(["cluster "+str(ref_peak.cluster), popt[0], perr[0]])

    return KD_matrix


def generate_xdata_ydata(ref,trace_list,cluster):
    """
    
    """    
    
    xdata = []
    ydata = []    
    
    for ref_peak,trace_peaks in give_all_clustered_peaks(ref,trace_list):

        ## TODO: can we truly assume this all the time?!
        xdata.append(0) #ref.Ltot_conc)
        ydata.append(0) #ref_peak.fractional_occupancy)
        ## alternative:
        # xdata.append(ref.Ltot_conc)
        # ydata.append(ref_peak.fractional_occupancy)
        ## -> doesn't work, because fractional_occupancies are not set for ref, so far...


        ## assumes that there is only one ref_peak with correct cluster number
        if ref_peak.cluster == cluster:
            
            for trace_peak in trace_peaks:        
                ## DEBUG
                #print(ref_peak,trace_peak)            
                        
                if trace_peak.footprinted_peak == True:
                    xdata.append(calculate_free_ligand_concentration(ref, [trace for trace in trace_list if trace_peak in trace.peaks][0]))
                    ydata.append(trace_peak.fractional_occupancy)

            break
            ## breaks the give_all_clustered_peaks loop because we already found the correct cluster
            
    return xdata, ydata
            

def plot_data(ref, trace_list, cluster):

    plt.title("Cluster " + str(cluster))
    plt.ylabel('Fractional Occupancy', fontsize = 16)
    plt.xlabel('Free Ligand Concentration', fontsize = 16)
    plt.ylim([-0.1, 1.1])
    plt.xlim([-1, 10.1])

    ## plot points
    xdata, ydata = generate_xdata_ydata(ref,trace_list,cluster)                
    plt.plot(xdata,ydata, 'o')


    ## plot line
    kd_matrix = fit_data_determine_kd(ref, trace_list)
    kd = [entry[1] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)]

    import numpy

    x = numpy.linspace(0,15,1000) # linearly spaced numbers
    y = fitFunc_fR(x,kd)
    plt.plot(x,y)
    
    ## show plot
    plt.show()
    
    return 


#def write_data_to_csv(KD_matrix):
  #return


if __name__ == "__main__":
    main(sys.argv[1:])