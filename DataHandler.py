######################
#
# @author Yasmin Heimann
#
# @Description:
# A class that handles the outputFile
######################

import OutputData
import numpy as np
##
# symbolic columns' names for easier use
##
EPITOPE_COL = 'Description'
MHC_COL = 'Allele Name'
AFFINITY_COL = 'Quantitative measurement'
MES_INEQ_COL = 'Measurement Inequality'
MES_SOURCE_COL = 'Method/Technique'
YEAR_COL = 'Date'
EC50_COL = 'Assay Group'
ERROR = 0
GROUP = 1
PREC_TO_KEEP = 0.8
PREC_DROPPED = "20"
YEAR_THRESHOLD = 2000
MAJORITY_CONST = 500

OUTPUT_FILE = True


class DataHandler:
    file_name = "iedb_data_"

    def __init__(self, ft, dup_data, single_data, g_num):
        self.__filter_type = ft
        self.__dup_data = dup_data
        self.__groups_num = g_num
        self.__single_data = single_data  # a df
        self.__filtered_groups = None
        self.__outputFile = OutputData.OutputData()
        self.__count = 1
        self.__total_cond_drops = 0

    def create_filtered_data(self):
        if self.__filter_type in ["YEAR2000", "EC50_FULL", "YEAR1995", "YEAR1990"]:
            self.filter_by_full_data_condition()
            self.__outputFile.addData(self.__filtered_groups)
            self.__outputFile.createOutput(self.file_name + self.__filter_type)
            return
        # build the output data
        self.__outputFile.addData(self.__single_data)
        # choose the method to filter the groups by
        if self.__filter_type == "HIGHEST":
            self.__filter_highest()
            for g in self.__filtered_groups:
                self.__outputFile.addData(g)
        elif self.__filter_type == "YEAR" or self.__filter_type == "EC50":
            self.__filter_by_condition()
            self.__outputFile.addReducedData(self.__filtered_groups)
        elif self.__filter_type == "MAJORITY_VOTE":
            self.__filter_by_majority_vote()
            self.__outputFile.addReducedData(self.__filtered_groups)
        # if data is of type FULL, only first and last step will be processed
        # extract the output file
        if OUTPUT_FILE:
            self.__outputFile.createOutput(self.file_name + self.__filter_type)

    def __filter_highest(self):
        data_with_errors = []
        before_size_of_groups = 0
        self.__filter_type = self.__filter_type + PREC_DROPPED
        for name, g in self.__dup_data:
            before_size_of_groups += len(g)
            data_with_errors.append((self.__calculate_group_error(g), g))
        # sort by error, smallest to largest
        data_with_errors = sorted(data_with_errors, key=lambda x:x[ERROR])
        sorted_groups = [tup[GROUP] for tup in data_with_errors]
        max_index = int(PREC_TO_KEEP * self.__groups_num)
        self.__filtered_groups = sorted_groups[:max_index]
        after_size_of_groups = 0
        for g in self.__filtered_groups:
            after_size_of_groups += len(g)
        print("from label ", self.__filter_type, PREC_TO_KEEP," : ",
              before_size_of_groups - after_size_of_groups,
              " were dropped")

    def __calculate_group_error(self, group):
        """
        Calculates the maximal error between the group's assays' affinity value
        :param group: a groupby object holding the parameters of identical assays
        :return: the maximal error
        """
        affinity_vals = group[AFFINITY_COL].tolist()
        size = len(affinity_vals)
        errors = []
        # go through all the pairs in the group, and calculate the error
        for i in range(size):
            for j in range(i + 1, size):
                errors.append(abs(affinity_vals[i] - affinity_vals[j]))
        return max(errors)

    def __get_measurement_type(self, ineq):
        temp = str(ineq)
        if temp == 'nan' or temp == "=":
            return "=", 'quantitative'
        else:
            return ineq, 'qualitative'

    def __filter_group_by_year(self, g, epitopes, mhcs, vals, ineqs, sources):
        print("filter by year ", self.__count)
        self.__count += 1
        group_final = []
        years = g[YEAR_COL].tolist()
        # choose which element to add and which to drop
        for i, y in enumerate(years):
            if y >= YEAR_THRESHOLD:
                ineqs[i], ineq_type = self.__get_measurement_type(ineqs[i])
                group_final.append([mhcs[i], epitopes[i], vals[i], ineqs[i],
                                    ineq_type, sources[i], mhcs[i]])
            else:
                self.__total_cond_drops += 1
        return group_final

    def __filter_group_by_EC_val(self, g, epitopes, mhcs, vals, ineqs, sources):
        group_final = []
        concentraion_vals = g[EC50_COL].tolist()
        # choose which element to add and which to drop
        for i, conc in enumerate(concentraion_vals):
            if 'EC50' not in conc:
                # find quantitative or qualitative val
                ineqs[i], ineq_type = self.__get_measurement_type(ineqs[i])
                group_final.append([mhcs[i], epitopes[i], vals[i], ineqs[i],
                                    ineq_type, sources[i], mhcs[i]])
            else:
                self.__total_cond_drops += 1
        return group_final

    def __filter_by_condition(self):
        data_matrix = []
        for name, g in self.__dup_data:
            epitopes = g[EPITOPE_COL].tolist()
            mhcs = g[MHC_COL].tolist()
            vals = g[AFFINITY_COL].tolist()
            ineqs = g[MES_INEQ_COL].tolist()
            sources = g[MES_SOURCE_COL].tolist()
            group_final = []
            # choose which element to add and which to drop
            if self.__filter_type == 'EC50':
                group_final = self.__filter_group_by_EC_val(g, epitopes, mhcs, vals, ineqs, sources)
            elif self.__filter_type == 'YEAR':
                group_final = self.__filter_group_by_year(g, epitopes, mhcs, vals, ineqs, sources)
            if len(group_final) > 0:
                data_matrix.extend(group_final)
        print("from label ", self.__filter_type, " : ", self.__total_cond_drops,
              " were dropped")
        self.__filtered_groups = np.array(data_matrix, dtype=object)

    def filter_by_full_data_condition(self):
        size_before = len(self.__single_data)
        if self.__filter_type == "YEAR2000":
            self.__filtered_groups = self.__single_data[(self.__single_data[YEAR_COL] >= 2000)]
        if self.__filter_type == "YEAR1995":
            self.__filtered_groups = self.__single_data[(self.__single_data[YEAR_COL] >= 1995)]
        elif self.__filter_type == "YEAR1990":
            self.__filtered_groups = self.__single_data[(self.__single_data[YEAR_COL] >= 1990)]
        elif self.__filter_type == "EC50_FULL":
            self.__filtered_groups = self.__single_data[~self.__single_data[EC50_COL].str.contains('EC50')]
        print("from label ", self.__filter_type, " : ", size_before - len(self.__filtered_groups),
                                                                     " were dropped")

    def __choose_values_barrier(self, vals):
        vals = np.array(vals)
        bool_vals = vals <= MAJORITY_CONST
        # sums up all the True labels, and check if they are the majority (vals <= const)
        if (vals <= MAJORITY_CONST).sum() / len(vals) >= 0.5:
            return bool_vals, True
        else:
            return bool_vals, False

    def __create_majority_final_group(self, g, vals, bool_vals, chosen_bool_val):
        # extract the values
        epitopes = g[EPITOPE_COL].tolist()
        mhcs = g[MHC_COL].tolist()
        ineqs = g[MES_INEQ_COL].tolist()
        sources = g[MES_SOURCE_COL].tolist()

        group_final = []
        # choose which element to add and which to drop
        for i, bool_val in enumerate(bool_vals):
            if bool_val == chosen_bool_val:
                ineqs[i], ineq_type = self.__get_measurement_type(ineqs[i])
                group_final.append([mhcs[i], epitopes[i], vals[i], ineqs[i],
                                    ineq_type, sources[i], mhcs[i]])
        return group_final

    def __filter_by_majority_vote(self):
        data_matrix = []
        dropped = 0
        for name, g in self.__dup_data:
            if len(g) > 2:
                original_size = len(g)
                # filter by majority condition
                vals = g[AFFINITY_COL].tolist()
                bool_vals, chosen_bool_val = self.__choose_values_barrier(vals)
                group_final = self.__create_majority_final_group(g, vals, bool_vals, chosen_bool_val)
                dropped += (original_size - len(group_final))
                data_matrix.extend(group_final)
        self.__filtered_groups = np.array(data_matrix, dtype=object)
        print("from label ", self.__filter_type, " : ", dropped,
              " were dropped")
