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
insilico.NP_RANDOM_SEED = constants_and_values["NP_RANDOM_SEED"]
insilico.NOISE_PEAK_ABSOLUTE_SD    = constants_and_values["NOISE_PEAK_ABSOLUTE_SD"]
insilico.NOISE_PEAK_RELATIVE_SD   = constants_and_values["NOISE_PEAK_RELATIVE_SD"]

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

#for cluster in [ref_peak.cluster for ref_peak in ref.peaks]:
#    fig, ax = footprint.plot_data(ref, trace_list, cluster, kd_matrix)
#    filename = "figures/plot-cluster{}.png".format(cluster)
#    fig.savefig(filename)
#    pass

#footprint.save_kd(kd_matrix, filename=OUT_FILENAME_RESULTS_CSV)

kdmat_expected = insilico.read_kdmatrix(IN_FILENAME_INSILICO_EXPECTED)
kdvalues_expected = insilico.get_kd_values(kdmat_expected)
df =  insilico.analyze_success(kdvalues_expected, kd_matrix)
df.to_csv(OUT_FILENAME_RESULTS_CSV+".long.csv")


##


kd_df = footprint.kdmatrix_to_df(kd_matrix)
kdvalues_expected.name="kd_expected"

pmatrix = footprint.fit_data_peaks_sensitive(ref, trace_list)
df_original = None
for alpha in [0.05, 0.025, 0.01, 0.005, 0.001]:
    classification_df = insilico.return_binary_classification_success(kdvalues_expected, pmatrix, alpha)
    statistics_dict = insilico.analyze_binary_classification_success(classification_df)
    statistics_dict.update( { "alpha" : alpha } )


    df = pmatrix.join( kdvalues_expected).join(kd_df).join(classification_df)
    sel_condition_positive  = ( df["condition"]  == True  )
    sel_prediction_positive = ( df["prediction"] == True  )
    true_positives = df[sel_condition_positive][sel_prediction_positive]
    true_positives["abs_diff_of_means"] = abs(true_positives["kd_expected"]-true_positives["KD mean"])
    true_positives["value_within_uncertainty"] = true_positives["abs_diff_of_means"] < true_positives["KD SD"]*insilico.COEFFICIENT_SD_ACCEPTED
    TP_within_uncertainty = true_positives["value_within_uncertainty"].sum()
    TP_outside_uncertainty = len(true_positives["value_within_uncertainty"]) - true_positives["value_within_uncertainty"].sum()
    true_positives["uncertainty_CV"] = true_positives["abs_diff_of_means"] / true_positives["KD SD"]
    true_positives["uncertainty_CV"].describe()
    
    statistics_dict.update( { "TP_within_uncertainty" : TP_within_uncertainty } )
    statistics_dict.update( { "TP_outside_uncertainty" : TP_outside_uncertainty } )

    def listify_dict_items(foo):
        for key in foo.keys():
            foo[key] = [ foo[key] ]
        return

    listify_dict_items(statistics_dict)

    if df_original is None:
        df_original = pandas.DataFrame.from_dict(statistics_dict)
    else:
        df_o        = pandas.DataFrame.from_dict(statistics_dict)
        df_original = df_original.append(df_o)


## 

import scipy.optimize

def cost_function(x):
    alpha= x #x[0]
    classification_df = insilico.return_binary_classification_success(kdvalues_expected, pmatrix, alpha)
    statistics_dict = insilico.analyze_binary_classification_success(classification_df)
    MCC = statistics_dict["MCC"]
    cost = 1 - MCC
    if numpy.isnan(cost):
        cost = float("inf")
    ##DEBUG
    _str = "x: {}, cost: {}, MCC: {}".format(x,cost,MCC)
    print(_str)
    return cost

def cost_function_logspace(x):
    x = 10 ** x
    return cost_function(x)

_ranges = [ slice(-7, 0.1, 1e-1) ]
_result = scipy.optimize.brute( func=cost_function_logspace, ranges=_ranges )

alpha = 10 ** _result[0]
classification_df = insilico.return_binary_classification_success(kdvalues_expected, pmatrix, alpha)
statistics_dict = insilico.analyze_binary_classification_success(classification_df)
statistics_dict.update( { "alpha" : alpha } )

df = pmatrix.join( kdvalues_expected).join(kd_df).join(classification_df)
sel_condition_positive  = ( df["condition"]  == True  )
sel_prediction_positive = ( df["prediction"] == True  )
true_positives = df[sel_condition_positive][sel_prediction_positive]
true_positives["abs_diff_of_means"] = abs(true_positives["kd_expected"]-true_positives["KD mean"])
true_positives["value_within_uncertainty"] = true_positives["abs_diff_of_means"] < true_positives["KD SD"]*insilico.COEFFICIENT_SD_ACCEPTED
TP_within_uncertainty = true_positives["value_within_uncertainty"].sum()
TP_outside_uncertainty = len(true_positives["value_within_uncertainty"]) - true_positives["value_within_uncertainty"].sum()
true_positives["uncertainty_CV"] = true_positives["abs_diff_of_means"] / true_positives["KD SD"]
true_positives["uncertainty_CV"].describe()
    
statistics_dict.update( { "TP_within_uncertainty" : TP_within_uncertainty } )
statistics_dict.update( { "TP_outside_uncertainty" : TP_outside_uncertainty } )


listify_dict_items(statistics_dict)
df_optimized = pandas.DataFrame.from_dict(statistics_dict)




df_out = df_original.append(df_optimized)
df_out.to_csv("statistics.csv")


###

