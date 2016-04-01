# -*- coding: utf-8 -*-

## workaround for relative imports...
import sys
import imp
from pathlib import Path # if you haven't already done so
root = Path(__file__).resolve().parents[1].as_posix()
# For older Python:
#   from os.path import dirname, realpath
#   root = dirname(dirname(realpath(__file__)))
sys.path.append(root)
## end workaround


import footprint
from footprint import Trace, Peak


#==============================================================================
# Import test set data 
#==============================================================================

## This is the data for this simple test case:
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

ref = trace_list[1]
target = trace_list[0]


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

    >>> conc_0_traces[0] == trace_list[1]
    ... ## should result in referencing to the second trace.
    True
    """
    
#==============================================================================
# cluster_peaks()
#==============================================================================

footprint.cluster_peaks(ref,[target])


def test_cluster_peaks():
    """
    >>> for p in ref.peaks: print(p.cluster,p.size_bp,p.peak_height)
    1 1.0 20.0
    2 2.0 40.0
    3 3.0 60.0

    >>> for t in target.peaks: print(t.cluster,t.size_bp,t.peak_height)
    1 1.07 10.0
    2 1.98 39.0
    3 3.01 61.0
    """
    for p in ref.peaks: print(p.cluster,p.size_bp,p.peak_height)
    for t in target.peaks: print(t.cluster,t.size_bp,t.peak_height)


#==============================================================================
# give_all_clustered_peaks()
#==============================================================================


def test_give_all_clustered_peaks():
    """
    >>> for i,j in footprint.give_all_clustered_peaks(ref,[target]): 
    ...     print(i.cluster,i.size_bp,j[0].size_bp)
    1 1.0 1.07
    2 2.0 1.98
    3 3.0 3.01
    
    >>> print(len(j))
    1
    """
    for i,j in footprint.give_all_clustered_peaks(ref,[target]):
        print(i.cluster,i.size_bp,j[0].size_bp)


#==============================================================================
# calculate_deviance_for_all_peaks()
#==============================================================================

def test_calculate_deviance_for_all_peaks():
    """
    >>> # doctest: +ELLIPSIS
    ... footprint.calculate_deviance_for_all_peaks(ref,target) 
    5.830...
    
    """
    print(footprint.calculate_deviance_for_all_peaks(ref,target))


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
    doctest.testmod()