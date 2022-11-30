######################
#
# @author Yasmin Heimann
#
# @Description:
# The program parses an IEDB data file of MHC class I binding affinity to experimented
# epitopes, given as a user input
# The file is composed of 7 main parameters: ID, Reference to the article, Epitope with
# description of the type of cell it was experimented on, Antigen Processing,
# MHC allele(class I in this case), Assay description,
# and the Quantitative Measure which is the affinity of the MHC allele to the given epitope.
#
# The program extracts a filtered csv of duplication - meaning an assay
# of the same mhc - allele
#
######################

import sys
import pandas as pd
from numpy import nan
import DataHandler
from error_analyzer import analyze_error, analyze_method_error

###
# Script's indicators of the wanted state - sample, debug, filtering mood etc.
###
DEBUG = False
SAMPLE = False
# True to include rows with blank inequality column, or False to include only "="
WITH_NAN = True
# decide whether to extract the filtered data into a file
TO_FILE = False
# output error analysis and graphs
ANALYZE_ERROR = False
# filter type to filter the output data by
TYPES = ["FULL", "EC50", "YEAR", "HIGHEST", "MAJORITY_VOTE", "YEAR2000", "EC50_FULL", "YEAR1990", "YEAR1995"]
FILTER_TYPE = "YEAR"
# indicates if we want to extract the full data into the NN format or not
FULL_DATA = False
# Creates a data set using the dataHandler class
EXTRACT_LEARNING_DATA = True
# analyse the method parametes numerically
ANALYZE_METHOD_ERROR = False

COL_ORDER = ['MHC ligand ID', 'Reference ID', 'Epitope ID', 'Description', 'Allele Name', 'Units', 'Qualitative Measure', 'Measurement Inequality', 'Quantitative measurement', 'PubMed ID', 'Date', 'Assay Comments', 'Type', 'Authors', 'Journal', 'Title', 'Submission ID', 'Object Type', 'Starting Position', 'Ending Position', 'Non-peptidic epitope ChEBI ID', 'Antigen Name', 'Parent Protein', 'Parent Protein Accession', 'Organism Name', 'Parent Species', 'Parent Species Accession', 'Epitope Comments', 'Epitope Relationship', 'Object Type.1', 'Description.1', 'Starting Position.1', 'Ending Position.1', 'Non-peptidic object Accession', 'Synonyms', 'Antigen Name.1', 'Parent Protein.1', 'Organism Name.1', 'Parent Organism', 'Name', 'Host ID', 'Geolocation', 'MHC Types Present', 'Process Type', 'Disease State', 'Disease Stage', 'Processed Antigen Epitope Relation', 'Processed Antigen Object Type', 'Processed Antigen Description', 'Processed Antigen Starting Position', 'Processed Antigen Ending Position', 'Non-peptidic Processed Antigen ChEBI ID', 'Processed Antigen Source Molecule Name', 'Processed Antigen protein parent Name', 'Processed Antigen protein parent Accession', 'Processed Antigen Organism Name', 'Processed Antigen Organism Species', 'Processed Antigen Organism Species ID', 'In vitro administration type', 'Processed Antigen Epitope Relation.1', 'Processed Antigen Object Type.1', 'Processed Antigen Description.1', 'Processed Antigen Starting Position.1', 'Processed Antigen Ending Position.1', 'Non-peptidic Processed Antigen ChEBI ID.1', 'Processed Antigen Source Molecule Name.1', 'Protein Parent Name', 'Protein Parent Accession', 'Processed Antigen Organism Name.1', 'Immunogen Organism Species', 'Immunogen Organism Species ID', 'Processed Antigen Comments', 'Location of assay data in the manuscript', 'Method/Technique', 'Assay Group', 'Number of Subjects Tested', 'Number of Subjects Responded', 'Response Frequency', 'PDB ID', 'Cell Tissue Type', 'Cell Type', 'Cell Culture Conditions', 'Allele Evidence Code', 'MHC allele class']


##
# symbolic columns' names for easier use
##
EPITOPE_COL = 'Description'
MHC_COL = 'Allele Name'
AFFINITY_COL = 'Quantitative measurement'
MEASUREMENT_SIGN_COL = 'Measurement Inequality'
AUTHOR_COL = 'Authors'


def filter_data_with_nan(data, res):
    # for "=" and NaN
    data[MEASUREMENT_SIGN_COL].replace("", nan, inplace=True)
    # filter all the None values
    data.dropna(subset=[AFFINITY_COL], inplace=True)
    conditions = ["=", nan]
    comp_data = data.loc[~data[MEASUREMENT_SIGN_COL].isin(conditions)]
    data = data.loc[data[MEASUREMENT_SIGN_COL].isin(conditions)]
    # add to results how many assyas have the = and have a measurement
    res.append("number of assays with the \"=\" or NaN(\"\") inequality: " + str(len(data)))
    # 40488 for only "=" sign filtering
    return data, comp_data, res[0]


def filter_data_with_eq(data, res):
    comp_data = []
    # for only "=" :
    data = data[data[MEASUREMENT_SIGN_COL] == "="]
    # filter all the None values
    data = data.dropna(subset=[AFFINITY_COL])
    # add to results how many assyas have the = and have a measurement
    res.append("number of assays with the \"=\" inequality: " + str(len(data)))
    # 115603 for "=" and "" filtering
    return data, comp_data, res[0]


def filter_data(data):
    res = []
    # filter rows with no ("=" or ("=" and "")) sign for the loaded measurement
    if WITH_NAN:
        return filter_data_with_nan(data, res)
    else:
        return filter_data_with_eq(data, res)


def extract_data_to_file(dup_data):
    # extract the data to csv file with different columns order (MHC,epitope at the begining)
    dup_data = dup_data[COL_ORDER]
    dup_data = dup_data.sort_values(by=[EPITOPE_COL, MHC_COL])
    dup_data.to_csv('duplications_full.csv')


def process_data_set(g, comp_data, total_dups, dup_data):
    # analyze the error of the data
    if ANALYZE_ERROR:
        analyze_error(g)
    # create output data
    if EXTRACT_LEARNING_DATA:
        output = DataHandler.DataHandler(FILTER_TYPE, g, comp_data, total_dups)
        output.create_filtered_data()
    # extract the data frame into csv
    if TO_FILE:
        extract_data_to_file(dup_data)
    # analyze error causes by method type
    if ANALYZE_METHOD_ERROR:
        analyze_method_error(g)


def update_results(results, res1, len_dup_data, len_none_dup_data,total_dups):
    results.append(res1)
    results.append("number of rows with duplications: " + str(len_dup_data))  # 23871 // 72326
    results.append("number of duplications is: " + str(total_dups))  # 11115 // 32013
    print('none dups len: ', len_none_dup_data)
    print('dups len: ', len_dup_data)
    return results


def createDataSets(data, results):
    if FULL_DATA:
        output = DataHandler.DataHandler(FILTER_TYPE, None, data, 0)
        output.create_filtered_data()
    else:
        # filter data by inequality ("=" "NaN" "<" ">")
        data, comp_data, res1 = filter_data(data)

        # divide the data: get duplicated rows, grouped together
        dup_data = data[data.duplicated(subset=[EPITOPE_COL, MHC_COL], keep=False)]
        non_dups_data = data[~data.duplicated(subset=[EPITOPE_COL, MHC_COL], keep=False)]
        comp_data = pd.concat([comp_data, non_dups_data], ignore_index=True)
        groups = dup_data.groupby([EPITOPE_COL, MHC_COL])
        total_dups = groups.ngroups # find the number of duplications groups
        results = update_results(results, res1, len(dup_data), total_dups, len(non_dups_data))

        process_data_set(groups, comp_data, total_dups, dup_data)
    return results


def clean_invalid_rows(data, results):
    """
    Clean rows that are invalid for the learning process
    :param data: df
    :param results: list of final resluts as strings
    :return: the clean data and results
    """
    # filters rows with no affinity value of of improper length
    data.dropna(subset=[AFFINITY_COL], inplace=True)
    data = data.loc[(data[EPITOPE_COL].map(len) >= 8) & (data[EPITOPE_COL].map(len) <= 15)]
    results.append("number of assays with no affinity value: " + str(len(data)))  # 163,098
    return data, results


def extract_data(path1, path2, sample, results):
    """
    create a data frame out of one or two csv files
    :param path1: main csv
    :param path2: csv if exists, o.w. None
    :param sample: indicator says that path2=None
    :return: a dataframe of the full data
    """
    # read the data where the second row is the columns, as it is the more detailed titles
    if not sample:
        data1 = pd.read_csv(path1, header=1, low_memory=False)
        data2 = pd.read_csv(path2, header=1, low_memory=False)
        data = pd.concat([data1, data2])
    else:
        data = pd.read_csv(path1, header=1, low_memory=False)
    results.append("total number of assays: " + str(len(data)))
    return data, results


def results_to_file(results_file, results):
    """
    Extracts the results log into a file
    :param results_file: indicator wether to extract or not
    :param results: the log to save
    """
    if results_file is not None:
        for line in results:
            print(line)
            if not FULL_DATA:
                results_file.write(line)
                results_file.write("\n")
    else:  # sample
        for line in results:
            print(line)


def find_duplications(path1, path2, sample, results_file):
    """
    The function creates a data frame of a csv file from the IEDB database, and creates a file with
    duplicated assays (same MHC, same epitope)
    It creates data sets for further learning, and extracts the data set to an excel file
    @path1 all assays until 2008 (inc)
    @path2 all assays from 2009  (inc)
    @sample a boolean variable indicates this run is on a sample data for debugging (True), or False
            for a run on the real data
    """
    results = []  # a list of results as string lines
    data, results = extract_data(path1, path2, sample, results)
    data, results = clean_invalid_rows(data, results)

    # create sets of data with the filtered results
    results = createDataSets(data, results)

    # write results to file
    results_to_file(results_file, results)


if __name__ == '__main__':
    if SAMPLE:
        find_duplications("original_sample.csv", "", True, None)
    else:
        if WITH_NAN:
            output_name = "results_blanks_" + FILTER_TYPE + ".txt"
        else:
            output_name = "results_eq_" + FILTER_TYPE + ".txt"
        with open(output_name, "w") as file:
            find_duplications("mhc_ligand_table_export_from09.csv",
                              "mhc_ligand_table_export_to08.csv",
                              False, file)
