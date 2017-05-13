# coding: utf-8

# In[1]:
# -------------------------------- #
### set plotting options
get_ipython().magic(u'matplotlib')
# -------------------------------- #


## SELECT YOUR LOGGING LEVEL ##
_logging_level = "DEBUG"
#_logging_level = "INFO"
#_logging_level = "WARNING"
#_logging_level = "CRITICAL"
## ------------------------  ##
import logging
_level = eval("logging."+_logging_level)
logging.basicConfig(level=_level)




# In[2]:
# -------------------------------- #
## import code for data processing
import peaks.footprint
footprint = peaks.footprint.Footprinter()
# -------------------------------- #


# In[20]:
# -------------------------------- #
### Import data of 'Peak Scanner 2'-Output and input_trace-file
trace_list = footprint.get_data("input_traces.csv")
# -------------------------------- #


# In[21]:
# -------------------------------- #
### Set options
## accepted_offset: peaks are clustered, if their size (=bp) lies within this difference (in bp)
## factor_method: if "num" look for optimal factor (free-floating), if "peak" find a footprinting-insensitive peak
## weight_smaller/weight_bigger: importance of misfit peaks which are smaller / bigger than reference peaks
## relative_mode: if "True" determine the importance of misfit peaks due to their relative size compared to the reference peak, if "False" consider absolute difference (in AU)
## from_bp/ to_bp: peaks in this interval of sizes (=bp) shall be considered for analysis.


_accepted_offset = 0.3
_factor_method   = "num"
_weight_smaller  = .5
_weight_bigger   = .5
_relative_mode   = False
_from_bp         = float("-inf")
_to_bp           = float("+inf")
_normalize_to    = 1000
_how_many_peaks_necessary = -1

### Generate the reference trace
ref = footprint.generate_averaged_negative_control(trace_list,accepted_offset=_accepted_offset,factor_method=_factor_method, weight_smaller=_weight_smaller, weight_bigger=_weight_bigger, relative_mode=_relative_mode, from_bp=_from_bp, to_bp=_to_bp, normalize_to=_normalize_to, how_many_peaks_necessary=_how_many_peaks_necessary)
# -------------------------------- #


# In[22]:
# -------------------------------- #
footprint.cluster_peaks(ref,trace_list,accepted_offset=_accepted_offset)
# -------------------------------- #
