#!/usr/bin/env python

import itertools
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
import statsmodels.formula.api


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

class Footprinter():
    global log

    def __init__(self):
        ## turn on logging
        logging.basicConfig(level=logging.DEBUG)
        self.log = logging.getLogger(__name__)
        return

    def main(self, argv=""):

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



    def list_traces(self, read_filelist="../HexA.csv"):
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
            header = next(csv_reader)
            index = Index()
            index.file_name = header.index('Sample File Name')
            index.sample_name = header.index('Sample Name')

            for row in csv_reader:
                if row[index.file_name] not in [t[0] for t in storage_traces]:
                    # Nope.
                    storage_traces.append([row[index.file_name],row[index.sample_name]])

        self.log.info("Writing read traces to output_traces.csv...")
        w=csv.writer(open("output_traces.csv","w"))
        w.writerow(["Sample ID","Sample File Name", "Dye", "Ltot", "Rtot"])
        for row in storage_traces:
            w.writerow([row[1],row[0],"?","?","?"])
        return storage_traces

    def get_data(self, config_file="input_traces.csv"):
        """Reads data from read_filelist, and returns an object containing all read
        information.

        This object consists of the following structure:

        [Trace(..., peaks = [Peak(...), Peak(...)], Trace(...)]
        See Trace and Peak to learn more about these classes. """

        NAME_sample_id = "Sample File Name"
        NAME_sample_file = "GetFromFile"

        import csv
        trace_list = []

        ## TODO: what entries to accept?
        with open(config_file) as g:
            csv_reader = csv.DictReader(g)
            sample_files = []
            for row in csv_reader:
                sample_files.append(row)

        storage_traces = []

        for i in sample_files:
            sample_name = i[NAME_sample_id]
            read_file = i[NAME_sample_file]

            self.log.info("Trying to add trace {} from file {}...".format(sample_name, read_file))

            with open(read_file) as f:
                csv_reader = csv.reader(f)
                header = next(csv_reader)
                index = Index()
                index.peak_height = header.index("Height")
                index.size_bp = header.index("Size")
                index.file_name = header.index('Sample File Name')

                for row in csv_reader:
                    if "Dye/Sample Peak" in header: #split_combined_fields ==
                        row.extend(row[header.index('Dye/Sample Peak')].split(","))
                        index.dye = len(header)
                        index.sample_peak = len(header)+1
                    else:
                        index.dye = header.index("Dye")
                        index.sample_peak = header.index("Sample Peak")

                    if row[index.file_name] == sample_name and "B" in row[index.dye]:
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

        for foo in [sample_file[NAME_sample_id] for sample_file in sample_files]:
            if foo not in storage_traces:
                self.log.critical("Couldn't find trace "+str(foo)+" in the given files...")

        return trace_list

    def generate_averaged_negative_control(self, trace_list, \
                                           accepted_offset=0.5, \
                                           factor_method="peak", \
                                           normalize_to=1000, \
                                           how_many_peaks_necessary=1, \
                                           *args, **kwargs):

        ## clean up trace_list, in case it was processed already.
        trace_list_ref = copy.deepcopy(trace_list)
        for t in trace_list_ref:
            for p in t.peaks:
                if "cluster" in vars(p):
                    del(p.cluster)
                    pass

        conc_0_traces = []
        for t in trace_list_ref:
            if t.Ltot_conc == 0:
                conc_0_traces.append(t)
        self.log.info("Found {} conc_0 traces.".format(len(conc_0_traces)))

        if len(conc_0_traces) > 1:
            i = 0

            for trace in conc_0_traces:
                self.log.debug("Screening trace {}...".format(trace.file_name))
                for ref_peak in trace.peaks:
                    self.log.debug("+ Screening peak {}:".format(ref_peak))
                    if "cluster" not in vars(ref_peak):
                        i += 1
                        ref_peak.cluster = i
                        for t in conc_0_traces:
                            self.log.debug("++ Screening trace {}...".format(t.file_name))
                            for peak in t.peaks:
                                self.log.debug("+++ Screening peak {}:".format(peak))
                                if abs(peak.size_bp - ref_peak.size_bp) < accepted_offset:
                                    self.log.debug("++++ Found a peak closeby!")
                                    if "cluster" not in vars(peak):
                                        self.log.debug("+++++ Identified a peak to cluster with.")
                                        peak.cluster = i
                                    else:
                                        del(ref_peak.cluster)
                                        if peak != ref_peak:
                                            self.log.critical("+++++ A peak is lying within the accepted offset for TWO clusters, not just one, as expected! This should not happen!")
                                            ref_peak.cluster = i
                                            self.log.critical("{} : {}".format(trace.file_name, ref_peak))
                                            self.log.critical("{} : {}".format(t.file_name, peak))
                                        else:
                                            ref_peak.cluster = i
            num_total_clusters = i
            self.log.info("Found {} potential clusters in total.".format(num_total_clusters))


            self.log.info("Trying to find clusters which shall be rejected...")

            for trace in conc_0_traces:
                index = conc_0_traces.index(trace)
                slice_conc_0_traces = conc_0_traces[:index] + conc_0_traces[index+1:]
                self.prune_tracepeaks_to_peaks_present_in_other_traces(trace, slice_conc_0_traces, how_many_peaks_necessary)

            storage_all_clusterids = []
            for trace in conc_0_traces:
                clusterids_in_this_trace = [peak.cluster for peak in trace.peaks if hasattr(peak, "cluster")]
                storage_all_clusterids.extend(clusterids_in_this_trace)
            storage_all_clusterids = set(storage_all_clusterids)

            # reset cluster index
            constant_shift = 100
            convert_lookup_list = list(enumerate(storage_all_clusterids,len(storage_all_clusterids)+constant_shift))

            for new_index, old_index in convert_lookup_list:
                for trace in conc_0_traces:
                    peaks_of_interest = [peak for peak in trace.peaks if getattr(peak, "cluster", None)==old_index]
                    for peak_of_interest in peaks_of_interest:
                        peak_of_interest.cluster = new_index

            for reset_index, (new_index, old_index) in enumerate(convert_lookup_list, 1):
                for trace in conc_0_traces:
                    peaks_of_interest = [peak for peak in trace.peaks if getattr(peak, "cluster", None)==new_index]
                    for peak_of_interest in peaks_of_interest:
                        peak_of_interest.cluster = reset_index


            num_total_clusters = len(storage_all_clusterids)
            self.log.info("... left with {} clusters.".format(num_total_clusters))

            num_total_traces = float(len(conc_0_traces))

            ## gradually fit optimal heights
            for x in range(5):
                self.log.debug("--- ROUND {} ---".format(x))
                ref = Trace("averaged_negative_control", "B", 0, 0, [])

                ## 1. create dummy peaks
                for i in range(1,num_total_clusters+1):
                    self.log.debug("Evaluate cluster {:3}".format(i))
                    ## fill peak storage
                    trace_storage = []
                    for trace in conc_0_traces:
                        for peak in trace.peaks:
                            if peak.cluster == i:
                                trace_storage.append(peak)

                    averaged_peak_height_mean = numpy.mean([p.peak_height for p in trace_storage])
                    averaged_peak_size_bp_mean = numpy.mean([p.size_bp for p in trace_storage])

                    ref.peaks.append(Peak(averaged_peak_size_bp_mean,averaged_peak_height_mean))


                ## 2. fit optimal factors
                self.cluster_peaks(ref, conc_0_traces, accepted_offset=accepted_offset)


                if factor_method == "num":
                    optimal_factors = self.determine_factor_numerically(ref, conc_0_traces, *args, **kwargs)
                elif factor_method == "peak":
                    optimal_factors = self.determine_factor_single_peak(ref, conc_0_traces, *args, **kwargs)

                ## 3. correct with factor
                for index, trace in enumerate(conc_0_traces):
                    self.correct_peaks_with_factor(trace,optimal_factors[index])

                ## normalize to largest peak

                #normalize_to  = max([peak.peak_height for peak in ref.peaks])
                #normalize_to /= 1000.0
                #for i in range(len(ref.peaks)):
                #    ref.peaks[i].peak_height              /= normalize_to


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
            self.log.info("You called generate_averaged_negative_control although only 1 trace with Ltot = 0 is available. Returning the clustered trace.")
            ## create new object
            c = conc_0_traces[0]
            ref = Trace(c.file_name, c.dye_color, c.Ltot_conc, c.Rtot_conc)
            for p in c.peaks:
                ref.peaks.append(Peak(p.size_bp, p.peak_height))
                ref.peaks[-1].averaged_peak_height_n = 1
                ref.peaks[-1].averaged_peak_height_nm = 1
                ref.peaks[-1].averaged_peak_height_sd = 0
                ref.peaks[-1].averaged_peak_size_bp_n = 1
                ref.peaks[-1].averaged_peak_size_bp_sd = 0
                ref.peaks[-1].averaged_peak_size_bp_nm = 1

            i = 0
            for ref_peak in ref.peaks:
                i += 1
                ref_peak.cluster = i

        else:
            self.log.critical("The provided trace_list contains no negative control traces with Ltot = 0.")

        ## normalize to largest peak
        self.normalize_trace_to_maxpeak(ref, normalize_to)

        ## TODO: is this necessary?
        ## clean-up
        del(trace_list_ref)

        ref.peaks = sorted(ref.peaks, key=lambda k: k.size_bp)

        return ref

    def collect_all_peaks(self, trace_list, \
                                       accepted_offset=0.5, \
                                       factor_method="peak", \
                                       normalize_to=1000, \
                                       how_many_peaks_necessary=1, \
                                       only_from_conc_0_traces = False,
                                       *args, **kwargs):

        ## clean up trace_list, in case it was processed already.
        trace_list_ref = copy.deepcopy(trace_list)
        for t in trace_list_ref:
            for p in t.peaks:
                if "cluster" in vars(p):
                    del(p.cluster)
                    pass

        conc_0_traces = []
        for t in trace_list_ref:
            if only_from_conc_0_traces:
                if t.Ltot_conc == 0:
                    conc_0_traces.append(t)
            else:
                conc_0_traces.append(t)
        self.log.info("Found {} conc_0 traces.".format(len(conc_0_traces)))

        if len(conc_0_traces) > 1:
            i = 0

            for trace in conc_0_traces:
                self.log.debug("Screening trace {}...".format(trace.file_name))
                for ref_peak in trace.peaks:
                    self.log.debug("+ Screening peak {}:".format(ref_peak))
                    if "cluster" not in vars(ref_peak):
                        i += 1
                        ref_peak.cluster = i
                        for t in conc_0_traces:
                            self.log.debug("++ Screening trace {}...".format(t.file_name))
                            for peak in t.peaks:
                                self.log.debug("+++ Screening peak {}:".format(peak))
                                if abs(peak.size_bp - ref_peak.size_bp) < accepted_offset:
                                    self.log.debug("++++ Found a peak closeby!")
                                    if "cluster" not in vars(peak):
                                        self.log.debug("+++++ Identified a peak to cluster with.")
                                        peak.cluster = i
                                    else:
                                        del(ref_peak.cluster)
                                        if peak != ref_peak:
                                            self.log.critical("+++++ A peak is lying within the accepted offset for TWO clusters, not just one, as expected! This should not happen!")
                                            ref_peak.cluster = i
                                            self.log.critical("{} : {}".format(trace.file_name, ref_peak))
                                            self.log.critical("{} : {}".format(t.file_name, peak))
                                        else:
                                            ref_peak.cluster = i
            num_total_clusters = i
            self.log.info("Found {} potential clusters in total.".format(num_total_clusters))


            self.log.info("Trying to find clusters which shall be rejected...")

            for trace in conc_0_traces:
                index = conc_0_traces.index(trace)
                slice_conc_0_traces = conc_0_traces[:index] + conc_0_traces[index+1:]
                self.prune_tracepeaks_to_peaks_present_in_other_traces(trace, slice_conc_0_traces, how_many_peaks_necessary)

            storage_all_clusterids = []
            for trace in conc_0_traces:
                clusterids_in_this_trace = [peak.cluster for peak in trace.peaks if hasattr(peak, "cluster")]
                storage_all_clusterids.extend(clusterids_in_this_trace)
            storage_all_clusterids = set(storage_all_clusterids)

            # reset cluster index
            constant_shift = 100
            convert_lookup_list = list(enumerate(storage_all_clusterids,len(storage_all_clusterids)+constant_shift))

            for new_index, old_index in convert_lookup_list:
                for trace in conc_0_traces:
                    peaks_of_interest = [peak for peak in trace.peaks if getattr(peak, "cluster", None)==old_index]
                    for peak_of_interest in peaks_of_interest:
                        peak_of_interest.cluster = new_index

            for reset_index, (new_index, old_index) in enumerate(convert_lookup_list, 1):
                for trace in conc_0_traces:
                    peaks_of_interest = [peak for peak in trace.peaks if getattr(peak, "cluster", None)==new_index]
                    for peak_of_interest in peaks_of_interest:
                        peak_of_interest.cluster = reset_index


            num_total_clusters = len(storage_all_clusterids)
            self.log.info("... left with {} clusters.".format(num_total_clusters))

            num_total_traces = float(len(conc_0_traces))

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
            self.log.info("You called generate_averaged_negative_control although only 1 trace with Ltot = 0 is available. Returning the clustered trace.")
            ## create new object
            c = conc_0_traces[0]
            ref = Trace(c.file_name, c.dye_color, c.Ltot_conc, c.Rtot_conc)
            for p in c.peaks:
                ref.peaks.append(Peak(p.size_bp, p.peak_height))
                ref.peaks[-1].averaged_peak_height_n = 1
                ref.peaks[-1].averaged_peak_height_nm = 1
                ref.peaks[-1].averaged_peak_height_sd = 0
                ref.peaks[-1].averaged_peak_size_bp_n = 1
                ref.peaks[-1].averaged_peak_size_bp_sd = 0
                ref.peaks[-1].averaged_peak_size_bp_nm = 1

            i = 0
            for ref_peak in ref.peaks:
                i += 1
                ref_peak.cluster = i

        else:
            self.log.critical("The provided trace_list contains no negative control traces with Ltot = 0.")

        ## TODO: is this necessary?
        ## clean-up
        del(trace_list_ref)

        ref.peaks = sorted(ref.peaks, key=lambda k: k.size_bp)

        return ref



    def normalize_trace_to_maxpeak(self, trace, normalize_to=1000):
        ## normalize to largest peak
        if normalize_to is None:
            pass
        else:
            normalize_to  = max([peak.peak_height for peak in trace.peaks]) / normalize_to
            for i in range(len(trace.peaks)):
                if hasattr(trace.peaks[i], "peak_height"):
                    trace.peaks[i].peak_height              /= normalize_to
                if hasattr(trace.peaks[i], "averaged_peak_height_sd"):
                    trace.peaks[i].averaged_peak_height_sd  /= normalize_to
        return

    def normalize_trace_to_sumpeak(self, trace, normalize_to=1000):
        ## normalize to largest peak
        if normalize_to is None:
            pass
        else:
            normalize_to  = sum([peak.peak_height for peak in trace.peaks]) / normalize_to
            for i in range(len(trace.peaks)):
                if hasattr(trace.peaks[i], "peak_height"):
                    trace.peaks[i].peak_height              /= normalize_to
                if hasattr(trace.peaks[i], "averaged_peak_height_sd"):
                    trace.peaks[i].averaged_peak_height_sd  /= normalize_to
        return

    def fit_all_traces_to_one(self, ref, trace_list, factor_method="num", *args, **kwargs):

        if factor_method   == "peak":
            optimal_factors = self.determine_factor_single_peak(ref, trace_list, *args, **kwargs)
        elif factor_method  == "num":
            optimal_factors = self.determine_factor_numerically(ref, trace_list, *args, **kwargs)

        self.log.debug(optimal_factors)
        self.correct_all_traces_with_factors(trace_list, optimal_factors)

        return


    def correct_all_traces_with_factors(self, trace_list, factor_list):
        # Correct traces with optimal factors
        optimal_factors = factor_list

        for index in range(len(trace_list)):
            self.correct_peaks_with_factor(trace_list[index],optimal_factors[index])
            self.log.debug("file: {:30.30} --> factor: {:4.2f}".format(trace_list[index].file_name, optimal_factors[index]))

        return


    def plot_height_size_overview_averaged_negative_control(self, ref, *args, **kwargs):
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

        return fig


    def plot_height_size_clustering_success_averaged_negative_control(self, ref, *args, **kwargs):

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

        return fig


    def cluster_peaks(self, ref, trace_list,accepted_offset=0.25):
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


    def calculate_deviance_for_all_peaks(self, ref, trace_list, weight_smaller=1,weight_bigger=1, relative_mode=False, from_bp=float("-inf"), to_bp=float("+inf")):
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

        weight_by_inverse_height = relative_mode

        warning_not_all_peaks_match = False
        warning_trace_range = False

        deviance_for_smaller_peaks = 0
        deviance_for_bigger_peaks = 0
        num_smaller=0
        num_bigger=0

        for ref_peak,trace_peaks in self.give_all_clustered_peaks(ref,trace_list):
            ## WORKAROUND for single trace mode
            # if there are no peaks clustered to the ref_peak, they cannot be included --> continue with next pair
            if trace_peaks == []:
                if warning_not_all_peaks_match == False:
                    self.log.info("Not all peaks match -- omitting deviation for non-comparable peaks.")
                    warning_not_all_peaks_match = True
                continue

            if ref_peak.size_bp < from_bp or ref_peak.size_bp > to_bp:
                if warning_trace_range == False:
                    self.log.info("Trace contains peaks which are not included in deviation calculation.")
                    warning_trace_range = True

            ## allows to calculate deviance for one trace only
            for trace_peak in trace_peaks:
                if from_bp < ref_peak.size_bp < to_bp:
                    if weight_by_inverse_height == True:
                        ## this mode calculates deviation as percentage point, with different
                        ## weights for each peak
                        weight = 1
                    elif weight_by_inverse_height == False:
                        ## similarly weighted, the difference between the two trace functions
                        ## ~ area / integral between traces is calculated
                        weight = 0
                    else:
                        weight = weight_by_inverse_height
                    counterweight = 1 - weight

                    deviation_abs = abs(ref_peak.peak_height - trace_peak.peak_height)
                    deviation_rel = deviation_abs / ref_peak.peak_height

                    deviation = weight*deviation_rel + counterweight*deviation_abs

                    if trace_peak.peak_height <= ref_peak.peak_height:
                        # deviation for smaller, i.e. potentially footprinted peaks
                        deviance_for_smaller_peaks += deviation
                        num_smaller +=1
                    else:
                        # deviation for bigger, i.e. potentially hypersensitive peaks
                        deviance_for_bigger_peaks += deviation
                        num_bigger +=1

        weighted_deviation = (weight_smaller*deviance_for_smaller_peaks + weight_bigger*deviance_for_bigger_peaks)

        ## possible further outputs:
        ## num_smaller
        ## num_bigger
        ## -> mean deviation per peak ~ normalizing

        return weighted_deviation



    def give_all_clustered_peaks(self, ref, trace_list, cluster=None):
        """
        Generates lists of clustered peaks from traces in trace_list, clustered to
        the peaks in ref.

        Yields (ref_peak, [trace_peak, trace_peak, ...]) elements.
        """

        if type(trace_list) is not list:
            trace_list = [trace_list]

        if cluster is None:
            clusters_to_search = set([peak.cluster for peak in ref.peaks])
        else:
            clusters_to_search = cluster

        for peak_cluster in clusters_to_search:
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


    def determine_factor_numerically(self, ref, trace_list, weight_smaller=1, weight_bigger=1, relative_mode=False, from_bp=float("-inf"), to_bp=float("+inf")):
        """
        Determines the optimal factor for trace, when compared to ref. This function
        minimizes the deviation as calculated by calculate_deviance_for_all_peaks()
        with the given options.

        Returns the optimal factor.
        ## TODO: implement real optimizer via scipy.optimize()
        """

        ## store original peak heights
        for trace in trace_list:
            for peak in trace.peaks:
                if "peak_height_original" not in vars(peak):
                    peak.peak_height_original = peak.peak_height
                else:
                    ## was already modified by factor
                    peak.peak_height_original = peak.peak_height_original
                peak.peak_height_tmp = peak.peak_height
                peak.peak_height = peak.peak_height_original

        optimal_factors = [1 for trace in trace_list]
        rmsd_old = self.calculate_deviance_for_all_peaks(ref, trace_list, weight_smaller, weight_bigger, relative_mode, from_bp, to_bp)

        self.log.info("starting: no calibration --> factors {!s} ; deviation {:8f}".format(optimal_factors, rmsd_old))

        for index, trace in enumerate(trace_list):

            def cost_function(x):
                _factor = x

                self.correct_peaks_with_factor(trace,_factor)
                rmsd_new = self.calculate_deviance_for_all_peaks(ref,trace,weight_smaller,weight_bigger,relative_mode, from_bp, to_bp)
                return rmsd_new

            x0 = [1.]
            min_result = scipy.optimize.minimize(cost_function, x0) #, method='nelder-mead', options={'disp': True})
            optimal_factor = min_result.x[0]

            self.log.info("found optimal: deviation {:8f} --> factor {:4.2f} ".format(cost_function(optimal_factor), optimal_factor))

            optimal_factors[index] = optimal_factor

        ##
        for index, trace in enumerate(trace_list):
            self.correct_peaks_with_factor(trace_list[index],optimal_factors[index])
            self.log.debug("calculating: trace # {:2} --> factor {:4.2f}".format(index, optimal_factors[index]))

        ## use calculate_deviance_for_all_peaks with trace and ref
        rmsd_new = self.calculate_deviance_for_all_peaks(ref,trace_list, weight_smaller, weight_bigger, relative_mode, from_bp, to_bp)

        self.log.info("found optimal: deviation {:8f} --> factors {!s} ".format(rmsd_new, optimal_factors))


        ## restore all peak heights to previous height
        for trace in trace_list:
            for peak in trace.peaks:
                peak.peak_height = peak.peak_height_tmp
                del(peak.peak_height_tmp)


        return optimal_factors

    def check_which_peaks_to_optimize(self, ref, trace_list, weight_smaller=1, weight_bigger=1, relative_mode=False, from_bp=float("-inf"), to_bp=float("+inf")):
        """

        """

        storage_dicts_to_panda = []

        for ref_peak,trace_peaks in self.give_all_clustered_peaks(ref,trace_list):
            if len(trace_peaks) != len(trace_list):
                self.log.debug("Peak at {:6.2f} bp is not found in all traces.".format(ref_peak.size_bp))
                self.log.debug("-- It is contained in {} traces, though. The missing ones are:".format(len(trace_peaks)))

                which_traces_are_contained = []

                for trace in trace_list:
                    for trace_peak in trace_peaks:
                        if trace_peak in trace.peaks:
                            which_traces_are_contained.append(trace)
                for trace in trace_list:
                    if not trace in which_traces_are_contained:
                        self.log.debug("-- {!s}".format(trace.file_name))
            else:
                which_traces_are_contained = trace_list[:]

            dict_to_panda = {}
            dict_to_panda.update({"cluster" : ref_peak.cluster})
            dict_to_panda.update({"size_bp" : ref_peak.size_bp})
            dict_to_panda.update({"number_of_traces_for_this_peak" : len(which_traces_are_contained)})
            dict_to_panda.update({"percentage_of_traces_for_this_peak" : float(len(which_traces_are_contained))/float(len(trace_list)) })


            for trace in trace_list:
                if trace in which_traces_are_contained:
                    dict_to_panda.update({trace.file_name : 1 })
                else:
                    dict_to_panda.update({trace.file_name : 0 })

            storage_dicts_to_panda.append(dict_to_panda)


        df = pandas.DataFrame.from_records(storage_dicts_to_panda)

        df_sum = df.sum()
        df_sum.name = "total"
        df = df.append(df_sum)

        df_relative = df_sum.copy()
        df_relative.name="total_percentage"
        df_relative /= float(len(ref.peaks))
        df = df.append(df_relative)

        dict_to_panda = {}
        for trace in trace_list:
            dict_to_panda.update({trace.file_name : trace.Ltot_conc })
        print(dict_to_panda)
        Ltots = pandas.DataFrame.from_records([dict_to_panda])
        Ltots.index = ["Ltot_conc"]
        df = df.append(Ltots)

        return df


    def check_which_trace_to_eliminate(self, ref, trace_list, weight_smaller=1, weight_bigger=1, relative_mode=False, from_bp=float("-inf"), to_bp=float("+inf")):
        which_traces_are_the_worst = numpy.asarray([0 for i in trace_list])

        for ref_peak, trace_peaks in self.give_all_clustered_peaks(ref,trace_list):

            which_traces_contained_this_peak = []
            for trace in trace_list:
                for trace_peak in trace_peaks:
                    if trace_peak in trace.peaks:
                        which_traces_contained_this_peak.append(trace)

            malus_array = numpy.asarray([0 for i in trace_list])
            for index, trace in enumerate(trace_list):
                if not trace in which_traces_contained_this_peak:
                    malus_array[index] = 1

            which_traces_are_the_worst += malus_array

        foo = zip(*sorted(list(enumerate(which_traces_are_the_worst)), key=lambda k: k[1], reverse=True))[0]
        sorted_tracelist_descending_order = [trace_list[index].file_name for index in foo]
        sorted_cumulated_malus_descending_order = [which_traces_are_the_worst[index] for index in foo]

        self.log.debug(sorted_tracelist_descending_order)
        self.log.debug(sorted_cumulated_malus_descending_order)

        return which_traces_are_the_worst


    def determine_factor_single_peak(self, ref, trace_list, weight_smaller=1, weight_bigger=1, relative_mode=False, from_bp=20, to_bp=170):
        """
        Determines the optimal factor for trace, when compared to ref. This function
        minimizes the deviation as calculated by calculate_deviance_for_all_peaks()
        with the given options. However, only factors corresponding to existent
        peak pairs are considered.

        Returns the optimal factor.

        """

        ## store original peak heights
        for trace in trace_list:
            for peak in trace.peaks:
                if "peak_height_original" not in vars(peak):
                    peak.peak_height_original = peak.peak_height
                else:
                    ## was already modified by factor
                    peak.peak_height_original = peak.peak_height_original
                peak.peak_height_tmp = peak.peak_height
                peak.peak_height = peak.peak_height_original

        optimal_factors = [1 for trace in trace_list]
        rmsd_old = self.calculate_deviance_for_all_peaks(ref, trace_list, weight_smaller, weight_bigger, relative_mode, from_bp, to_bp)
        which_peak = None

        self.log.info("starting: original peak heights = factors {!s} --> deviation {:8f}".format(optimal_factors, rmsd_old))

        for ref_peak,trace_peaks in self.give_all_clustered_peaks(ref,trace_list):

            ## WORKAROUND for single trace mode
            if trace_peaks == []:
                continue

            ## if not all peaks are clustered on this ref_peak, they cannot be used --> continue with next cycle
            if len(trace_peaks) != len(trace_list):
                self.log.debug("Skipping peak at {:6.2f} bp because it is not found in all traces.".format(ref_peak.size_bp))
                continue

            ## if ref_peak is inside range, follow the program = pass; else: omit this entry --> continue with next cycle.
            if from_bp < ref_peak.size_bp < to_bp:
                pass
            else:
                continue

            self.log.debug("Using peak at {:6.2f} bp ...".format(ref_peak.size_bp))


            testing_factors = []

            ## works with one trace at a time
            for index, trace_peak in enumerate(trace_peaks):
                ## have to rely on trace_peaks coming in same order as trace_list ...
                factor =  ref_peak.peak_height / trace_peak.peak_height_original
                self.correct_peaks_with_factor(trace_list[index],factor)
                testing_factors.append(factor)
                ## DEBUG
                self.log.debug("testing: trace # {:2} on cluster {:3} @ {:8.1f} bp --> factor {:4.2f}".format(index, trace_peak.cluster, trace_peak.size_bp, factor))

            ## use calculate_deviance_for_all_peaks with trace and ref
            rmsd_new = self.calculate_deviance_for_all_peaks(ref,trace_list, weight_smaller, weight_bigger, relative_mode, from_bp, to_bp)



            ## DEBUG
            self.log.debug("was testing: all traces on cluster {:3} @ {:8.1f} bp --> deviation {:8f}".format(ref_peak.cluster, ref_peak.size_bp, rmsd_new))

            ## compare deviance_new with deviance_old, if better deviance_new -> deviance_old, else delete deviance_new
            if rmsd_new <= rmsd_old:
                rmsd_old = rmsd_new
                optimal_factors = testing_factors
                which_peak = ref_peak

        if which_peak == None:
            self.log.info("no peak would be better than current factorization --> factors {!s} ; deviation {:8f}".format( optimal_factors, rmsd_old))
        else:
            self.log.info("found optimal: cluster {:3} @ {:8.1f} bp ; deviation {:8f} --> factors {!s} ".format(which_peak.cluster, which_peak.size_bp, rmsd_old, optimal_factors))

        ## restore all peak heights to original height
        for trace in trace_list:
            for peak in trace.peaks:
                peak.peak_height = peak.peak_height_tmp
                del(peak.peak_height_tmp)


        return optimal_factors


    def correct_peaks_with_factor(self, trace, factor):
        """
        Corrects all peak_heights in trace with factor. The multiplication is applied
        directly to the object given as argument.
        """

        for peak in trace.peaks:
            if not hasattr(peak, "peak_height_original"):
                peak.peak_height_original = peak.peak_height
            peak.peak_height = peak.peak_height_original * factor

        return


    def mark_footprinted_peaks(self, ref, trace_list, threshold=0.1, mark_all=False):
        """
        Marks potentially footprinted peaks, i.e. peaks to be evaluated, by setting
        the flag footprinted_peak, if the peak_heights of ref and trace differ by
        threshold*ref_peak.peak_height (standard value = 10%) and if the ligand
        concentrations of the compared traces differ from each other.

        Sets the flag to True, if the conditions above are fulfilled, otherwise
        to False.

        Modifies the object trace directly.
        """

        for ref_peak,trace_peaks in self.give_all_clustered_peaks(ref,trace_list):

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

    def add_fractional_occupancies(self, ref,trace_list):
        """
        TODO:    ...

        Modifies the object trace directly.
        """

        for ref_peak,trace_peaks in self.give_all_clustered_peaks(ref,trace_list):
            for trace_peak in trace_peaks:
                trace_peak.fractional_occupancy = 1 - trace_peak.peak_height / ref_peak.peak_height

        return


    def calculate_free_ligand_concentration(self, ref,trace,mode="Ltot"):
        """
        TODO:    ...

        Calculates the free ligand concentration.

        Returns Lfree_conc.
        """

        ## // BROKEN
        #for ref_peak,trace_peaks in self.give_all_clustered_peaks(ref,trace):
            ## WORKAROUND for singlecalculate_free_ligand_concentration trace mode
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


    def fitFunc_fR(self, Lfree_conc, Kd):
        ## underlying formula:
        ## y = fR
        ## x = Lfree_conc
        ## to be fitted: Kd

        fR = (Lfree_conc)/(Lfree_conc + Kd)

        return fR


    def fit_data_determine_kd(self, ref, trace_list):
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
        for cluster in [ref_peak.cluster for ref_peak in ref.peaks]:
            _result_dict = self.fit_data_for_one_cluster_kd(ref, trace_list, cluster)
            _ = self.append_to_kdmatrix(KD_matrix, _result_dict, cluster, trace_list)

        return KD_matrix

    def append_to_kdmatrix(self, KD_matrix, _result_dict, cluster, trace_list):
        ## TODO: optimize this form.
        _kd          = _result_dict["kd"]
        _err         = _result_dict["std"]
        _len_ydata   = _result_dict["len_ydata"]

        appendix = ["cluster "+str(cluster), _kd, _err, _len_ydata, float(_len_ydata)/float(len(trace_list))]
        self.log.debug(appendix)
        KD_matrix.append(appendix)

        return KD_matrix


    def fit_data_for_one_cluster_kd(self, ref, trace_list, cluster, excluded_points=[]):

        self.log.debug("Evaluating cluster {}...".format(cluster))
        ref_peak, trace_peaks = list(self.give_all_clustered_peaks(ref, trace_list, [cluster]))[0]

        if ref_peak==[] and trace_peaks==[]:
            self.log.critical("Got empty peak list for the upcoming cluster {}. Can't do that, sorry! :(".format(cluster))
        else:
            self.log.debug("Evaluating cluster {}...".format(ref_peak.cluster))
            self.log.debug("DEBUG DETAIL: {}".format([ref_peak,trace_peaks]))

        xdata, ydata = self.generate_xdata_ydata(ref,trace_list,cluster)
        ixy = list(enumerate(zip(xdata,ydata)))
        xy = [(x,y) for i,(x,y) in ixy if i not in excluded_points]
        xdata, ydata = zip(*list(xy))

        ## fitting...
        try:
            popt, pcov= curve_fit(self.fitFunc_fR, xdata, ydata)
            # compute SE, i.e. standard deviation errors
            perr = numpy.sqrt(numpy.diag(pcov))
            _kd = popt[0]
            _err = perr[0]
        except RuntimeError:
            _kd, _cov, _err = numpy.nan, numpy.nan, numpy.nan

        _return_dict = {}
        _return_dict["kd"]          = _kd
        _return_dict["std"]         = _err
        _return_dict["len_ydata"]   =  len(ydata)

        return _return_dict


    def fit_data_peaks_sensitive(self, ref, trace_list, clusters_to_check=None):
        """
        Determines the K_D by fitting of observed fractional occupancies against
        calculated free ligand concentration.

        """

        ## data has to be in the following format:
        ## xdata = [Lfree_conc1, Lfree_conc2, ...]
        ## ydata = [fractional_occupancy1, fractional_occupancy2, ...]

        ## output:
        ## [ K_D(cluster1), K_D(cluster2), ... ]

        sm = statsmodels.formula.api

        p_matrix = []

        for ref_peak, trace_peaks in self.give_all_clustered_peaks(ref, trace_list, clusters_to_check):

            if ref_peak==[] and trace_peaks==[]:
                self.log.critical("Got empty peak list for the upcoming cluster. Can't do that, sorry! :(")
                continue
            else:
                self.log.debug("Evaluating cluster {} for footprinting sensitivity...".format(ref_peak.cluster))
                #self.log.debug("DEBUG DETAIL: {}".format([ref_peak,trace_peaks]))

            xdata = []
            ydata = []

            for trace_peak in trace_peaks:
                if trace_peak.footprinted_peak == True:
                    xdata.append(self.calculate_free_ligand_concentration(ref, [trace for trace in trace_list if trace_peak in trace.peaks][0]))
                    ydata.append(trace_peak.fractional_occupancy)

                ## TODO: shall we catch out clusters with no footprinted peaks at all?
            df_dict = {
                       "x" : xdata,
                       "y" : ydata
                       }

            df = pandas.DataFrame.from_dict(df_dict)
            res = sm.ols(formula="y ~ x - 1", data=df).fit()

            appendix = []

            appendix.append("cluster {}".format(str(ref_peak.cluster)))
            appendix.append(res.pvalues["x"])
            appendix.append(res.params["x"])
            appendix.append(res.bse["x"])

            ## save results...
            ## TODO: optimize this form.

            self.log.debug(appendix)
            p_matrix.append(appendix)


        p_df = pandas.DataFrame.from_records(zip(*p_matrix), index=["cluster #", "pvalue", "coeff", "se"])
        return p_df.T


    def generate_xdata_ydata(self,ref,trace_list,cluster):

        xdata = []
        ydata = []
        i = 0

        for ref_peak,trace_peaks in self.give_all_clustered_peaks(ref,trace_list):

            ## assumes that there is only one ref_peak with correct cluster number
            if ref_peak.cluster == cluster:
                for trace_peak in trace_peaks:
                    if trace_peak.footprinted_peak == True:
                        xdata.append(self.calculate_free_ligand_concentration(ref, [trace for trace in trace_list if trace_peak in trace.peaks][0]))
                        ydata.append(trace_peak.fractional_occupancy)
                break
        return xdata, ydata


    def get_quality_characteristics(self, ref, trace_list, cluster, kd_matrix):
            qc = {}
            qc["bp"]                                      = float([ref_peak.size_bp for ref_peak in ref.peaks if ref_peak.cluster == cluster][0])
            qc["n"]                                       = int([entry[3] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)][0])
            qc["n/m"]                                     = float([entry[4] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)][0])
            qc["kd"]                                      = float([entry[1] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)][0])
            qc["std"]                                     = float([entry[2] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)][0])
            qc["std/kd"]                                  = qc["std"] / qc["kd"]


            qc["ref_height"]                              = float([ref_peak.peak_height for ref_peak in ref.peaks if ref_peak.cluster == cluster][0])
            qc["ref_height_sd"]                           = float([ref_peak.averaged_peak_height_sd for ref_peak in ref.peaks if ref_peak.cluster == cluster][0])
            qc["ref_peak_quality"]                        = float([ref_peak.averaged_peak_height_nm for ref_peak in ref.peaks if ref_peak.cluster == cluster][0])
            qc["peaks_of_this_cluster"]                   = list(self.give_all_clustered_peaks(ref, trace_list, [cluster]))[0][1]
            qc["original_peak_heights_of_this_cluster"]   = [peak.peak_height_original for peak in qc["peaks_of_this_cluster"]]
            qc["min_peak"]                                = min(qc["original_peak_heights_of_this_cluster"])
            qc["max_peak"]                                = max(qc["original_peak_heights_of_this_cluster"])
            return qc


    def save_quality_characteristics(self, ref, trace_list, kd_matrix, filename="quality-characteristics.csv"):
        import csv
        with open(filename,"w") as f:
            # _fieldnames = qc.keys()
            _fieldnames = ["cluster", "bp", "kd", "std", "min_peak", "max_peak", "n", "n/m", "std/kd"]

            w=csv.DictWriter(f, _fieldnames, extrasaction="ignore")
            w.writeheader()
            for cluster in [eval(i[0].replace("cluster ","")) for i in kd_matrix]:
                print(repr(cluster))
                _d = self.get_quality_characteristics(ref, trace_list, cluster, kd_matrix)
                _d["cluster"] = cluster
                w.writerow(_d)
        return

    def plot_data(self, ref, trace_list, cluster, kd_matrix, with_point_numbers=False, no_show=False, excluded_points=[]):

        if no_show:
            plt.ioff()

        fig, ax = plt.subplots(1)

        plt.title("Cluster " + str(cluster))
        plt.ylabel('Fractional Occupancy', fontsize = 16)
        plt.xlabel('Free Ligand Concentration', fontsize = 16)
        plt.ylim([-1.1, 1.1])
        plt.xlim([-1, 10.1])

        ## get data
        xdata, ydata = self.generate_xdata_ydata(ref,trace_list,cluster)
        ixy = list(enumerate(zip(xdata,ydata)))

        ## plot points
        for i,(x,y) in ixy:
            if i in excluded_points:
                ax.plot(x,y, marker='o', color="grey")
            else:
                ax.plot(x,y, marker='o', color="blue")

        ## plot point numbers, if whished for
        if with_point_numbers == True:
            for i,(x,y) in ixy:
                if i in excluded_points:
                    ax.text(x,y,i, bbox=dict(facecolor='grey', alpha=0.5))
                else:
                    ax.text(x,y,i, bbox=dict(facecolor='red', alpha=0.5))

        ## plot line
        kd = [entry[1] for entry in kd_matrix if entry[0] == "cluster "+str(cluster)]

        import numpy

        x = numpy.linspace(0,15,1000) # linearly spaced numbers
        y = self.fitFunc_fR(x,kd)
        ax.plot(x,y)

        ## show plot
        plt

        qc = self.get_quality_characteristics(ref, trace_list, cluster, kd_matrix)

        textstr = """\
        n   = {n}
        n/m = {n/m:3.2f}
        min_peakheight = {min_peak:5.1f} A.U.
        max_peakheight = {max_peak:5.1f} A.U.
        $K_D$ = {kd:3.2f} $\pm$ {std:3.2f} $\mu$M

        ref_pos = {bp:4.1f} bp
        ref_height = {ref_height:5.1f} $\pm$ {ref_height_sd:5.1f} A.U.
        ref_quality = {ref_peak_quality:3.2f}
        """
        textstr = textstr.format(**qc)
        textstr = "\n".join([foo.strip() for foo in textstr.splitlines()])


        self.log.debug(textstr)

        # place a text box in upper left in axes coords
        props = dict(facecolor='wheat', alpha=0.3)

        # place a text box in upper left in axes coords
        ax.text(0.95, 0.05, textstr, horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, fontsize=12, family="monospace",
            bbox=props, color=(0,0,0,0.5) )

        return fig, ax

    def kdmatrix_to_df(self,kd_matrix):
        df = pandas.DataFrame.from_records(kd_matrix, columns=["Cluster #", "KD mean", "KD SD", "n", "m/n"])
        return df



    def save_kd(self,kd_matrix, filename="kd-matrix.csv"):
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


    def plot_peakscan(self,ref, trace_list, *args, **kwargs):

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

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(13,15))

        for label, df2 in df.groupby('trace'):
            if label=="averaged_negative_control":
                df2.plot(x="cluster", y="peak_height", kind="line", ax=ax, label=label, marker="s", linewidth=2)
            else:
                df2.plot(x="cluster", y="peak_height", kind="line", ax=ax, label=label)



        leg = plt.legend(loc="upper right", ncol=2)
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)

        if "ylim" in kwargs:
            plt.ylim(kwargs["ylim"])
        if "xlim" in kwargs:
            plt.xlim(kwargs["xlim"])

        plt.show()

        return


    def prune_tracepeaks_to_peaks_present_in_other_traces(self, trace, trace_list, how_many=-1):
        cluster_list = list(self.give_all_clustered_peaks(trace, trace_list))

        if how_many < 0:
            how_many = len(trace_list)

        for trace_peak, other_traces in cluster_list:
            if len(other_traces) >= how_many:
                pass
            else:
                trace.peaks.remove(trace_peak)

        ## hunt for cluster==0 peaks
        def should_keep(trace_peak):
            if trace_peak.cluster != 0:
                return True
            else:
                return False
        for trace in trace_list:
            trace.peaks[:] = [trace_peak for trace_peak in trace.peaks if should_keep(trace_peak)]

        return



    def prune_tracepeaks_to_peaks_within_bp_limits(self, trace, from_bp=float("-inf"), to_bp=float("+inf") ):

        def should_keep(trace_peak):
            if from_bp <= trace_peak.size_bp <= to_bp:
                return True
            else:
                return False

        trace.peaks[:] = [trace_peak for trace_peak in trace.peaks if should_keep(trace_peak)]

        return

    def make_round_comparison(self, trace_list, refit=False, *args, **kwargs):
        l = len(trace_list)
        dev_result = numpy.ndarray((l,l))*numpy.nan
        for indexlist, trace_combination in zip( list(itertools.product(range(l), repeat=2)) , \
                                                 list(itertools.product(trace_list, repeat=2)) ):
            t,u = trace_combination
            t = copy.deepcopy(t)
            u = copy.deepcopy(u)
            self.prune_tracepeaks_to_peaks_present_in_other_traces(u, [t])
            if refit != False:
                self.log.debug("Refitting trace {} to {}".format(u.file_name, t.file_name))
                if refit == True or refit == "num":
                    optfactor = self.determine_factor_numerically(t, [u], *args, **kwargs)
                elif refit == "peak":
                    optfactor = self.determine_factor_single_peak(t, [u], *args, **kwargs)
                self.log.debug("Optimal factor: {}".format(optfactor[0]))
                self.correct_peaks_with_factor(u, optfactor[0])
            dev = self.calculate_deviance_for_all_peaks(t, u, *args, **kwargs)
            dev_result[indexlist] = dev
        return dev_result

    def tracelist_to_csv(self, trace_list, filename="dataframe.csv"):
        df_big = self.tracelist_to_dataframe(trace_list)
        df_big = self.dataframe_add_differences(df_big)
        df_big.to_csv(filename)
        return


    def tracelist_to_dataframe(self, trace_list):
        df_big=None
        for index, trace in enumerate(trace_list):
            bigdict={}
            clusterlist=[]
            heightlist=[]
            for peak in trace.peaks:
                if getattr(peak, "cluster", 0)!=0:
                    clusterlist.append(peak.cluster)
                    heightlist.append(peak.peak_height)
            bigdict.update({"cluster" : clusterlist, trace.file_name+" ({})".format(index) : heightlist})
            df_small = pandas.DataFrame.from_dict(bigdict)
            if not df_big is None:
                df_big = df_big.merge(df_small, on="cluster", how="outer")
            else:
                df_big=df_small
        return df_big

    def dataframe_add_differences(self, dataframe):
        ## calculate differences
        df_big = dataframe
        cols = list(dataframe.columns)
        cols = [col for col in cols if "(" in col and ")" in col]

        for indexlist, trace_combination in zip( list(itertools.combinations(range(len(cols)), r=2)) , \
                                                         list(itertools.combinations(cols, r=2)) ):
            tf,uf=trace_combination
            i1,i2=indexlist
            df_big["({}-{})".format(i1,i2)] = df_big[tf] - df_big[uf]
            df_big["({}-{})/{}".format(i1,i2,i1)] = ( df_big[tf] - df_big[uf] ) / df_big[tf]

        return df_big

    def unite_peaks_to_one(self, ref, trace_list, cluster_tuple, block_name=""):
        cluster_list = cluster_tuple
        if block_name=="":
            block_name = tuple(cluster_list)
        for trace in [ref] + trace_list:
            allhs = [(peak.peak_height, peak.size_bp) for peak in trace.peaks if getattr(peak, "cluster", None) in cluster_list]

            if len(allhs) > 0:
                pass
            else:
                ## = if no peaks are found at all in this trace, we can skip to the next trace and do not need to add a virtual peak
                continue

            allh, alls = zip(*allhs)
            h = sum(allh)
            s = numpy.mean(alls)
            p = Peak(s, h)

            p.cluster = block_name
            _1 = [ getattr(peak, "averaged_peak_size_bp_sd",None) for peak in trace.peaks if getattr(peak, "cluster", None) in cluster_list]
            _2 = [ getattr(peak,"averaged_peak_size_bp_nm",None) for peak in trace.peaks if getattr(peak, "cluster", None) in cluster_list]
            _3 = [ getattr(peak,"averaged_peak_height_n",None) for peak in trace.peaks if getattr(peak, "cluster", None) in cluster_list]
            _4 = [ getattr(peak,"averaged_peak_height_nm",None) for peak in trace.peaks if getattr(peak, "cluster", None) in cluster_list]
            _5 = [ getattr(peak,"averaged_peak_height_sd",None) for peak in trace.peaks if getattr(peak, "cluster", None) in cluster_list]
            _6 = [ getattr(peak,"averaged_peak_size_bp_n",None) for peak in trace.peaks if getattr(peak, "cluster", None) in cluster_list]
            _7 = [ getattr(peak,"peak_height_original",None) for peak in trace.peaks if getattr(peak, "cluster", None) in cluster_list]
            if None not in _1:
                    p.averaged_peak_size_bp_sd = sum(_1)
            if None not in _2:
                p.averaged_peak_size_bp_nm = numpy.mean(_2)
            if None not in _3:
                p.averaged_peak_height_n = sum(_3)
            if None not in _4:
                p.averaged_peak_height_nm = numpy.mean(_4)
            if None not in _5:
                p.averaged_peak_height_sd = sum(_5)
            if None not in _6:
                p.averaged_peak_size_bp_n = sum(_6)
            if None not in _7:
                p.peak_height_original = sum(_7)

            trace.peaks.append(p)
        return
