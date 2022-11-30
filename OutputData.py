######################
#
# @author Yasmin Heimann
#
# @Description:
# A class that handles the outputFile
######################

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from collections import Counter

##
# symbolic columns' names for easier use
##
EPITOPE_COL = 'Description'
MHC_COL = 'Allele Name'
OUTPUT_MHC_COL = "allele"
AFFINITY_COL = 'Quantitative measurement'
MES_INEQ_COL = 'Measurement Inequality'
MES_SOURCE_COL = 'Method/Technique'
CROSSV_CONST = 5
CROSSV_PREC = 0.05

OUTPUT_COLS = ["allele","peptide","measurement_value","measurement_inequality",
               "measurement_type","measurement_source","original_allele"]


class OutputData:

    def __init__(self):
        self.__final_data = None
        self.__firstReduced = True
        self.__alelle_specific_train_set = pd.DataFrame(columns=OUTPUT_COLS)  # df of the train sets of mhcs
        self.__alelle_specific_test_set = pd.DataFrame(columns=OUTPUT_COLS)  # df of the test
        self.__first = True
        self.__counter = 1

    def addData(self, group):
        """
        :param group: a groupby object/DF holding assays' information
        """
        print("adding", self.__counter)
        self.__counter += 1
        epitopes = group[EPITOPE_COL].tolist()
        mhcs = group[MHC_COL].tolist()
        vals = group[AFFINITY_COL].tolist()
        ineqs = group[MES_INEQ_COL].tolist()
        sources = group[MES_SOURCE_COL].tolist()
        mes_types = []
        for i, ineq in enumerate(ineqs):
            ineq = str(ineq)
            if ineq == 'nan':
                ineqs[i] = "="
                mes_types.append('quantitative')
            elif ineq == '=':
                mes_types.append('quantitative')
            else:
                mes_types.append('qualitative')
        mat = np.array([mhcs, epitopes, vals, ineqs, mes_types, sources,mhcs])
        if self.__first:
            self.__final_data = mat.transpose()
            self.__first = False
        else:
            self.__final_data = np.concatenate((self.__final_data, mat.transpose()))

    def addReducedData(self, data):
        """
        :param data: a np.matrix object holding the relevant info from the iedb assays
        """
        if self.__final_data is None:
            self.__final_data = data
        elif len(data) > 0:
            self.__final_data = np.concatenate((self.__final_data, data))

    def createOutput(self, name):
        self.__create_train_test_data()
        file_name = name + ".csv"
        # create test data
        self.__alelle_specific_test_set.to_csv("test_" + file_name, index=False)
        # create mhc test sets
        self.__alelle_specific_train_set.to_csv("train_" + file_name, index=False)

    def __create_allele_data(self, mhc, count):
        print("creating test and train for allele", mhc, "with", count)
        # create a filtered df
        if count == 1:  # only one peptide for this mhc, can be only learned about
            self.__alelle_specific_train_set.append(self.__final_data.loc[self.__final_data.allele == mhc],
                                                    ignore_index=True)
        else:
            prec = CROSSV_PREC
            if count < 5:  # take half of the data as test
                prec = 1 / count
            mhc_assays = self.__final_data.loc[self.__final_data.allele == mhc]
            train, test = train_test_split(mhc_assays, test_size=prec)
            self.__alelle_specific_train_set = self.__alelle_specific_train_set.append(train, ignore_index=True)
            self.__alelle_specific_test_set = self.__alelle_specific_test_set.append(test, ignore_index=True)

    def __create_train_test_data(self):
        # sample fifth of the data (20%)
        self.__final_data = pd.DataFrame(self.__final_data, columns=OUTPUT_COLS)
        all_mhcs = self.__final_data[OUTPUT_MHC_COL].tolist()
        distinct_mhcs_dict = Counter(all_mhcs)
        for mhc, count in distinct_mhcs_dict.items():
            self.__create_allele_data(mhc, count)

