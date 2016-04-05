# -*- coding: utf-8 -*-

## workaround for relative imports...
import sys
from pathlib import Path # if you haven't already done so
root = Path(__file__).resolve().parents[1].as_posix()
# For older Python:
#   from os.path import dirname, realpath
#   root = dirname(dirname(realpath(__file__)))
sys.path.append(root)
## end workaround


import footprint


#==============================================================================
# Import test set data 
#==============================================================================

import generate_large_data_set
import numpy as np


## set conditions for data set generation


rtot = 0.1

## easy case
#peaks_bp = np.arange(1,3,1)
#footprinted_peaks_bp_and_kd = [[1,1]]
#lfree_concs = np.arange(1,15,5)

sd_bp = 0.1
sd_height = 5

## advanced case
#==============================================================================
peaks_bp = np.arange(1,20,1)
footprinted_peaks_bp_and_kd = [[3,0.3],[4,0.2],[5,0.1],[6,0.5],[7,0.8],[8,0.6],[9,0.2],[10,0.2],[11,0.3],[12,0.4],[13,0.5],[14,0.6],[15,0.7],[16,0.6],[17,0.5],[18,0.4]]#number varied
#lfree_concs = np.arange(0.1,15,3) ## TODO: multiple times
lfree_concs = np.random.random_integers(1,15,10) ## gets 30 random numbers in [1,15] (varied)
#==============================================================================


trace_list = generate_large_data_set.generate_data_set(rtot, peaks_bp, footprinted_peaks_bp_and_kd , lfree_concs, sd_bp, sd_height)

ref = trace_list[0]
traces = trace_list[1:]
trace = trace_list[1]


print("Successfully set up test case data.")

#==============================================================================
# Check if 0M traces are recognized correctly.
#==============================================================================

conc_0_traces = []
for t in trace_list:
    if t.Ltot_conc == 0:
        conc_0_traces.append(t)

def test_recognition_of_0M_traces():
    """
    >>> len(conc_0_traces) == 1
    True

    >>> conc_0_traces[0] == trace_list[0]
    ... ## should result in referencing to the second trace.
    True
    """
    
#==============================================================================
# cluster_peaks()
#==============================================================================

footprint.cluster_peaks(ref,traces)


def test_cluster_peaks():
    """
    >>> for p in ref.peaks: print(p.cluster,p.size_bp,p.peak_height)
    1 1.0 20.0
    2 2.0 40.0
    3 3.0 60.0

    >>> for t in trace.peaks: print(t.cluster,t.size_bp,t.peak_height)
    1 1.07 10.0
    2 1.98 39.0
    3 3.01 61.0
    """
    print('cluster')
    for p in ref.peaks: print(p.cluster,p.size_bp,p.peak_height)
    for t in trace.peaks: print(t.cluster,t.size_bp,t.peak_height)


#==============================================================================
# give_all_clustered_peaks()
#==============================================================================
footprint.give_all_clustered_peaks(ref, trace_list)


def test_give_all_clustered_peaks():
    """
    >>> for i,j in footprint.give_all_clustered_peaks(ref,traces): 
    ...     print(i.cluster,i.size_bp,[foo.size_bp for foo in j])
    1 1.0 1.07
    2 2.0 1.98
    3 3.0 3.01
    
    >>> print(len(j))
    1
    """
    print('give clustered peaks')
    for i,j in footprint.give_all_clustered_peaks(ref,traces):
        print(i.cluster,i.size_bp,[foo.size_bp for foo in j])


#==============================================================================
# calculate_deviance_for_all_peaks()
#==============================================================================

footprint.calculate_deviance_for_all_peaks(ref, trace, weight_smaller=1,weight_bigger=0,relative_mode=False, from_bp=0, to_bp=130,)


def test_calculate_deviance_for_all_peaks():
    """
    >>> # doctest: +ELLIPSIS
    ... footprint.calculate_deviance_for_all_peaks(ref, trace) 
    7.106...
    
    """
    print(footprint.calculate_deviance_for_all_peaks(ref,trace))


#==============================================================================
# correct_peaks_with_factor(trace, factor)
#==============================================================================
#def test_correct_peaks_with_factor_1(trace, factor):
   # """
   # >>> original_trace = trace
   # >>> footprint.correct_peaks_with_factor(trace, 1)
   # >>> footprint.correct_peaks_with_factor(trace, 2)
   # >>> footprint.correct_peaks_with_factor(trace, 0.5)
   # >>> original_trace == trace
    #True
   # """
    #footprint.correct_peaks_with_factor(trace, 1)
    #footprint.correct_peaks_with_factor(trace, 2)
    #footprint.correct_peaks_with_factor(trace, 0.5)


#==============================================================================
# determine_factor_numerically():
#==============================================================================
footprint.determine_factor_numerically(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False)


def test_determine_factor_numerically(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False):
    """
    >>> footprint.determine_factor_numerically(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False)
    1.01
    """
    print(footprint.determine_factor_numerically(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False))


#==============================================================================
# determine_factor_single_peak()
#==============================================================================
def test_determine_factor_single_peak(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False):
    """
    >>> footprint.determine_factor_single_peak(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False)
    1.0256410256410255
    """
    print(footprint.determine_factor_single_peak(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False))


#==============================================================================
# correct_peaks_with_factor
#==============================================================================
factor = footprint.determine_factor_numerically(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=False)
footprint.correct_peaks_with_factor(trace, factor)

def test_correct_peaks_with_factor():
    print('Corrected:')

for p in ref.peaks: print(p.cluster,p.size_bp,p.peak_height)
for t in trace.peaks: print(t.cluster,t.size_bp,t.peak_height)
    


#==============================================================================
#  mark_footprinted_peaks()
#==============================================================================

for t in traces:
    print(t)
    footprint.mark_footprinted_peaks(ref, t)

def test_mark_footprinted_peaks(ref, trace):
    """
    >>> footprint.mark_footprinted_peaks(ref, traces, threshold=0.1)
    >>> for ref_peak,trace_peaks in footprint.give_all_clustered_peaks(ref, traces):
    ...    print([trace_peak.footprinted_peak for trace_peak in trace_peaks])
    True
    False
    False
    """


#==============================================================================
# add_fractional_occupancies()
#==============================================================================

footprint.add_fractional_occupancies(ref,traces)

#for t in traces:
    #print([peak.fractional_occupancy for peak in t.peaks])
## -> doesn't work if peaks are not clustered together! -> 'Peak' object has no attribute 'fractional_occupancy'



#==============================================================================
# calculate_free_ligand_concentration()
#==============================================================================

def test_calculate_free_ligand_concentration():
    """
    >>> print(footprint.calculate_free_ligand_concentration(ref,trace))
    4.95
    """
print("Lfree : "+str(footprint.calculate_free_ligand_concentration(ref,trace)))


#==============================================================================
# fit_data_determine_kd()
#==============================================================================
footprint.fit_data_determine_kd(ref, traces)

def test_fit_data_determine_kd():
    """
    >>> footprint.fit_data_determine_kd(ref, traces)
    [['cluster 1', 4.9500000000000002, 0.0], ['cluster 2', 1.0, inf], ['cluster 3', 1.0, inf]]
    """

print(footprint.fit_data_determine_kd(ref, traces))




footprint.plot_data(ref, traces, 4)
footprint.plot_data(ref, traces, 5)
footprint.plot_data(ref, traces, 6)
#footprint.plot_data(ref, traces, 7)


## TODO: rest of tests.

#
#    factor = determine_factor_numerically(ref, trace, 1, 0)
#    print("optimal factor : "+str(factor))
#    correct_peaks_with_factor(trace,factor)
#
#    print("Corrected peak pairs: ")
#    for i,j in give_all_clustered_peaks(ref,trace): print(i.cluster,i.size_bp,i.peak_height,j.peak_height)
#
#    add_fractional_occupancies(ref,trace)
#    print("Lfree : "+str(calculate_free_ligand_concentration(ref,trace)))
#    ## DEBUG
#    #print([peak.cluster for peak in peak_list])



print("Done.")

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
