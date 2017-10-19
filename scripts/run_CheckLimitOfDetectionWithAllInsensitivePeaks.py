# coding: utf-8

import matplotlib 
matplotlib.use('Agg')

import json_tricks as pickler
import pandas
import numpy

class Container():
    def __init__(self):
        return

container = Container()

##
IN_FILENAME_PICKLED = "parameters_and_constants.json"
IN_FILENAME_INSILICO_EXPECTED = "kd-matrix-expected.csv"
IN_FILENAME_INSILICO_TEMPLATE = "input_traces.csv"

OUT_FILENAME_RESULTS_CSV = "summary_performance.csv"
##

with open(IN_FILENAME_PICKLED, "rb") as f:
    constants_and_values = pickler.load(f)
    ## WORKAROUND for OrderedDict import of json_tricks
    if pickler.__name__ == "json_tricks":
        constants_and_values = dict(constants_and_values)

import footprint
footprint = footprint.Footprinter()

_accepted_offset = constants_and_values["accepted_offset"]
_factor_method   = constants_and_values["factor_method"] 
_weight_smaller  = constants_and_values["weight_smaller"]
_weight_bigger   = constants_and_values["weight_bigger"]
_relative_mode   = constants_and_values["relative_mode"]
_from_bp         = constants_and_values["from_bp"]
_to_bp           = constants_and_values["to_bp"]


import insilico
insilico.NP_RANDOM_SEED          = constants_and_values["NP_RANDOM_SEED"]
insilico.NOISE_PEAK_ABSOLUTE_SD  = constants_and_values["NOISE_PEAK_ABSOLUTE_SD"]
insilico.NOISE_PEAK_RELATIVE_SD  = constants_and_values["NOISE_PEAK_RELATIVE_SD"]
insilico.PEAK_HEIGHT_FROM        = constants_and_values["PEAK_HEIGHT_FROM"]
insilico.PEAK_HEIGHT_TO          = constants_and_values["PEAK_HEIGHT_TO"]

insilico.__init__()



trace_list = insilico.get_data(IN_FILENAME_INSILICO_EXPECTED, IN_FILENAME_INSILICO_TEMPLATE)

ref = footprint.generate_averaged_negative_control(trace_list,accepted_offset=_accepted_offset,factor_method=_factor_method, weight_smaller=0.5, weight_bigger=0.5, relative_mode=_relative_mode, from_bp=_from_bp, to_bp=_to_bp)

footprint.cluster_peaks(ref,trace_list,accepted_offset=_accepted_offset)

#optimal_factors = footprint.determine_factor_numerically(ref, trace_list, weight_smaller=_weight_smaller, weight_bigger=_weight_bigger, relative_mode=_relative_mode, from_bp=_from_bp, to_bp=_to_bp)
#
#for index, null in enumerate(trace_list):
#    footprint.correct_peaks_with_factor(trace_list[index],optimal_factors[index])
#    print("file: {:30.30} --> factor: {:4.2f}".format(trace_list[index].file_name, optimal_factors[index]))
footprint.fit_all_traces_to_one(ref, trace_list, factor_method=_factor_method,  weight_smaller=_weight_smaller, weight_bigger=_weight_bigger, relative_mode=_relative_mode, from_bp=_from_bp, to_bp=_to_bp)


footprint.mark_footprinted_peaks(ref, trace_list, threshold=0.01, mark_all=True) 
footprint.add_fractional_occupancies(ref,trace_list)
for trace in trace_list:
    Lfree = footprint.calculate_free_ligand_concentration(ref,trace,mode="Ltot")
    trace.Lfree = Lfree
    del(Lfree)   
    print("file: {:30.30} : Ltot: {:5.2f} , Lfree_conc: {:5.2f}".format(trace.file_name, trace.Ltot_conc, trace.Lfree))

kd_matrix = footprint.fit_data_determine_kd(ref, trace_list)

footprint.save_kd(kd_matrix, filename="kd-matrix-found.csv")


###
statistics_dict = {}

df = pandas.read_csv("kd-matrix-found.csv")

mini = df["KD mean"].min()
statistics_dict.update( { "kd_min" : mini }  )


def listify_dict_items(foo):
    for key in foo.keys():
        foo[key] = [ foo[key] ]
    return

listify_dict_items(statistics_dict)

df_out = pandas.DataFrame.from_dict(statistics_dict)
df_out.to_csv("statistics.csv")
###

