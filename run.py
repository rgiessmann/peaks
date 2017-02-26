
# coding: utf-8

# In[1]:

get_ipython().magic(u'matplotlib nbagg')


# In[2]:

import peaks.footprint as footprint


# In[20]:

# -------------------------------- #
### Import data of 'Peak Scanner 2'-Output and input_trace-file
trace_list = footprint.get_data("work/magda/MS031_HexA.csv", "work/magda/input_traces_MS031_HexA.csv")
# -------------------------------- #


# In[21]:

reload(footprint)

length_substrate = 328

_from_bp = 50
_to_bp =   200

# -------------------------------- #
### Generate the reference trace
ref = footprint.generate_averaged_negative_control(trace_list,accepted_offset=0.5,factor_method="peak", weight_smaller=1, weight_bigger=1, relative_mode=False, from_bp=_from_bp, to_bp=_to_bp)
## Vary: accepted_offset -> see clustering; weight_smaller/weight_bigger: weighing of deviance for smaller and bigger peaks;
##       from_bp/ to_bp: size_bp area that should be considered. Leave out full fragment peaks, small size_bp peaks and saturated peaks
##       num: True -> determine_factor_numerically, False -> determine_factor_single_peak


# In[22]:

conc_0_traces = []
for t in trace_list:
    if t.Ltot_conc == 0:
        conc_0_traces.append(t)

_0 = conc_0_traces[0]
_1 = conc_0_traces[1]
_2 = conc_0_traces[2]

_2.peaks
        
footprint.cluster_peaks(ref, conc_0_traces)            


# In[23]:


for ref_peak in ref.peaks:
        print("cluster: {:3} ref_peak_size: {:8.2f} bp;   ref_peak_height: {:8.1f} AU;    n/N: [size: {:3.2f}, height: {:3.2f}]".format(ref_peak.cluster, ref_peak.size_bp, ref_peak.peak_height, ref_peak.averaged_peak_size_bp_nm,ref_peak.averaged_peak_height_nm ))
# -------------------------------- #


# In[24]:

get_ipython().magic(u'matplotlib inline')


# In[25]:

reload(footprint)
footprint.plot_height_size_clustering_success_averaged_negative_control(ref)


# In[26]:

reload(footprint)
footprint.plot_height_size_overview_averaged_negative_control(ref, ylim=(0,600))


# In[27]:

reload(footprint)
footprint.plot_peakscan(ref,[ref], ylim=(0,600))


# In[ ]:




# In[28]:

# -------------------------------- #
### Cluster the Peaks of the same bp_size of all traces 
footprint.cluster_peaks(ref,trace_list,accepted_offset=0.35)
## Vary accepted_offset to allow bigger or smaller differences in size_bp within each cluster

for ref_peak in ref.peaks:
    print("cluster: {:3} ;   peak size: {:8.1f} bp".format(ref_peak.cluster , ref_peak.size_bp))
# -------------------------------- #


# In[29]:

reload(footprint)


# In[ ]:




# In[30]:

# -------------------------------- #
### Determine factor for height adjustment and then adjust traces with it to one height
#factor = footprint.determine_factor_numerically(ref, trace, weight_smaller=1, weight_bigger=1, relative_mode=True, from_bp=20, to_bp=103)
optimal_factors = footprint.determine_factor_single_peak(ref, trace_list, weight_smaller=.1, weight_bigger=1.5, relative_mode=True, from_bp=20, to_bp=103)

for index, null in enumerate(trace_list):
    footprint.correct_peaks_with_factor(trace_list[index],optimal_factors[index])
    print("file: {:30.30} --> factor: {:4.2f}".format(trace_list[index].file_name, optimal_factors[index]))
# -------------------------------- #


# In[15]:

# -------------------------------- #
### Calculate fractional occupancy for all footprinted peaks, then calculate the ligand concentration (modes Ltot/ Lfree)   
footprint.mark_footprinted_peaks(ref, trace_list, threshold=0.01, mark_all=True) 
footprint.add_fractional_occupancies(ref,trace_list)

for trace in trace_list:
    Lfree = footprint.calculate_free_ligand_concentration(ref,trace,mode="Ltot")
    trace.Lfree = Lfree
    del(Lfree)   
    print("file: {:30.30} : Ltot: {:5.2f} , Lfree_conc: {:5.2f}".format(trace.file_name, trace.Ltot_conc, trace.Lfree))
# -------------------------------- #    


# In[16]:

# -------------------------------- #
### Calculate K_d by fitting  fractional_occupancy = L/(L + K_d) to data		
kd_matrix = footprint.fit_data_determine_kd(ref, trace_list)
# -------------------------------- #


# In[17]:

# -------------------------------- #
### Automatic plotting for every cluster
reload(footprint)
for cluster in [ref_peak.cluster for ref_peak in ref.peaks]:
    footprint.plot_data(ref, trace_list, cluster, kd_matrix)
# -------------------------------- #


# In[60]:

for trace in trace_list:
    print trace.Lfree
    print [peak for peak in trace.peaks if peak.cluster==29]


# In[ ]:




# In[ ]:

# -------------------------------- #
### Manual plotting for a specific cluster
#cluster = 1
#footprint.plot_data(ref, trace_list, cluster, kd_matrix)                                                                          
# -------------------------------- #


# In[28]:

# -------------------------------- #
# Save K_d-Matrix -> cluster, size_bp, K_d, Error & number of data points 
footprint.save_kd(kd_matrix, filename="kd-matrix_MS035_036_Ltot.csv")
# -------------------------------- #


# In[ ]:



