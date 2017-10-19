# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 17:57:07 2016

@author: rgiessmann
"""


# generate a large data set with random errors

import numpy as np
import footprint



#==============================================================================
# generating traces
#==============================================================================

## !! WARNING: systematic error!! Ltot of trace is set to Lfree ... critical for low concentrations!!
def generate_data_set(rtot = 0.1, peaks_bp = np.arange(1,3,1), footprinted_peaks_bp_and_kd = [[1,1]], lfree_concs = np.arange(1,15,5), sd_bp = 0.1, sd_height = 5):
    """
    ## !! WARNING: systematic error!! Ltot of trace is set to Lfree ... critical for low concentrations!!
    """
    
    trace_list = []
    
    ## for control trace:
    trace_list.append(footprint.Trace(file_name = "0", dye_color = "B", Ltot_conc = 0, Rtot_conc = rtot, peaks=[]))
    ref=trace_list[0]
    
    for p in peaks_bp:
        peak_height = np.random.random_integers(50,6000)
        peak_bp = p
        ref.peaks.append(footprint.Peak(peak_bp, peak_height))
    
    
    ## for footprinting traces
    for l in lfree_concs:
        trace_list.append(footprint.Trace(file_name = str(l), dye_color = "B", Ltot_conc = l, Rtot_conc = rtot, peaks=[]))
    
        for p in peaks_bp:
            if p in [i[0] for i in footprinted_peaks_bp_and_kd]:
                kd = [i[1] for i in footprinted_peaks_bp_and_kd if p == i[0]]
                ref_peak_height = [pe.peak_height for pe in ref.peaks if pe.size_bp == p][0]
                factor = 1 - l/(l+kd)            
            else:
                ref_peak_height = [pe.peak_height for pe in ref.peaks if pe.size_bp == p][0]
                factor = 1
                
            peak_height = ( ref_peak_height * factor ) + np.random.normal(0,sd_height)
            peak_bp = np.random.normal(p, sd_bp)
            trace_list[-1].peaks.append(footprint.Peak(peak_bp, peak_height))
    
    ## DEBUG        
    #print(trace_list)
    
    return trace_list
