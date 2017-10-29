#!/usr/bin/env python
# coding: utf-8

import argparse
import io
import os
import pandas
import sys


def main(argv=None):
    if argv is None:
        argv = sys.argv

    ## prepare ArgumentParser
    parser = argparse.ArgumentParser(description="Reads a .csv file full of analyzed peaks and annotates it with raw peak heights from a collection of .csv files.", epilog="", add_help=True)

    parser.add_argument(dest="csv_files", metavar="<raw_peak_height_files>", nargs="+", help="the .csv file with raw peak heights and data point numbers.")
    parser.add_argument("-i", "--in", metavar="<analyzed_file>", dest="input", required=True, help="the .csv file containing the results from PeakScanner: analyzed peaks.")
    parser.add_argument("-o", "--out", metavar="<output_filename>", default="<analyzed_file>.annotated.csv", help="the filename of the annotated .csv file.")
    #    parser.add_argument("csv_files", metavar="<raw_peak_height_files>", nargs="*", help="additional .csv file with raw peak heights and data point numbers.")

    ## processes args from ArgumentParser
    args = parser.parse_args(argv[1:])

    input_file = args.input
    extracted_csv_files = args.csv_files
    output_file = args.out

    if output_file == parser.get_default("out"):
        output_file = input_file + ".annotated.csv"

    ## check validity of args
    assert(os.path.exists(input_file))
    #assert(not os.path.exists(output_file))

    ## data processing
    df = pandas.read_csv(input_file)
    df2 = pandas.DataFrame()

    columns_to_merge_on = ["Data Point", "Sample File Name"]

    print('''In total: {} entries to annotate for all sample files.'''.format(len(df)))

    for i in extracted_csv_files:
        assert(os.path.exists(i))
        sample_file_name = i[:-4]+".fsa"

        ## contains all Raw Peak Heights of one trace
        new_df = pandas.read_csv(i)
        new_df.columns=["Data Point", "Raw Height"]
        new_df["Sample File Name"] = sample_file_name

        _df = df.merge(new_df, on=columns_to_merge_on)

        ## quality checking
        how_many_should_i_annotate = sum(df["Sample File Name"] == sample_file_name)
        how_many_did_i_annotate = len(_df)
        print('''sample file {:^30}: number of entries to annotate: {:>4}; successfully annotated: {:>4}.'''.format(sample_file_name, how_many_should_i_annotate, how_many_did_i_annotate))

        df2 = df2.append(_df, ignore_index=True)


    intersection = df.join(df2, how="outer", rsuffix="_y")
    could_not_find_corresponding_raw_peak_height = intersection["Raw Height"].isnull()

    if could_not_find_corresponding_raw_peak_height.sum() > 0:
        print("Could not find corresponding raw peak heights for the following entries:")
        print(df[could_not_find_corresponding_raw_peak_height])
        print('''In total, I could not find corresponding raw peak height for {} entries.'''.format(sum(intersection["Raw Height"].isnull())))
    else:
        print("Found corresponding raw peak heights for all given entries.")


    ## write out file
    df2.to_csv(output_file)

    print("Have written annotated results to {}.".format(output_file))

if __name__ == "__main__":
    main()
