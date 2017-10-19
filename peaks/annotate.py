#!/usr/bin/env python
# coding: utf-8

import argparse
import io
import pandas
import os


def main():
    ## prepare ArgumentParser
    parser = argparse.ArgumentParser(description="Reads a .csv file full of analyzed peaks and annotates it with raw peak heights from a collection of .csv files.", epilog="")

    parser.add_argument("-i", "--in", metavar="<analyzed_file>", dest="input", required=True, help="the .csv file containing the results from PeakScanner: analyzed peaks.")
    parser.add_argument("csv_files", metavar="<raw_peak_height_file>", nargs="*", help="the .csv files with raw peak heights and data point numbers.")
    parser.add_argument("-o", "--out", metavar="<output_filename>", default="all_out.csv", help="the filename of the annotated .csv file.")

    ## processes args from ArgumentParser
    args = parser.parse_args()

    input_file = args.input
    extracted_csv_files = args.csv_files
    output_file = args.out

    #input_file = "../analyzed_BH429_all-RTG.csv"

    #extracted_csv_files = [
    #    "110-12-17-12-14 PM.csv",
    #    "210-12-17-1-01 PM.csv",
    #    "310-12-17-1-35 PM.csv",
    #    "410-12-17-2-10 PM.csv",
    #    "510-12-17-2-45 PM.csv",
    #    "610-12-17-3-19 PM.csv"
    #    ]

    df = pandas.read_csv(input_file)
    df2 = pandas.DataFrame()

    for i in extracted_csv_files:
        sample_file_name = i[:-4]+".fsa"

        new_df = pandas.read_csv(i)
        new_df.columns=["Data Point", "Raw Height"]
        new_df["Sample File Name"] = sample_file_name

        columns_to_merge_on = ["Data Point", "Sample File Name"]

        _df = df.merge(new_df, on=columns_to_merge_on)
        df2 = df2.append(_df, ignore_index=True)

    df2.to_csv(output_file)

if __name__ == "__main__":
    main()
