import numpy as np
import pandas as pd
from sklearn import metrics
from scipy import stats
from collections import Counter
from matplotlib import pyplot as plt
import operator

PREDICT_COL = 'predicted labels'
TRUE_COL = 'true labels'
MINIMAL_VALUE = 1E-30
NORMALIZATION = np.log(50000)
MHC_COL = 'mhc allele'
MOST_COMMON_NUM = [5, 10, 15, 20, 25, 30, 35, 40]

ALLELE_INDEX = 0
TRAIN_INDEX = 1
TEST_INDEX = 2


class Error_Calc:
    def __init__(self, data_paths_list, labels_list):
        self.__data_list = [pd.read_csv(path) for path in data_paths_list]
        self.__labels = labels_list

    def addDataPath(self, path, label):
        self.__data_list.append(pd.read_csv(path))
        self.__labels.append(label)

    def __normalize_values(self, l1):
        return list(map(lambda val: (np.log(val + 1) / NORMALIZATION), l1))

    def __data_loss(self, data, type):
        predicts, true = self.__normalize_values(np.array(data[PREDICT_COL])), \
                         self.__normalize_values(np.array(data[TRUE_COL]))
        if type == 'Pearson':
            return stats.pearsonr(predicts, true)
        if type == "MSE":
            self.__greatest_n_loss(data, type)
            return metrics.mean_squared_error(true, predicts)

    def loss(self, type):
        if type == 'Pearson':
            for i, data in enumerate(self.__data_list):
                print("Data type: ", self.__labels[i], " Model Pearson correlation: ",
                      self.__data_loss(data, type))
        if type == "MSE":
            for i, data in enumerate(self.__data_list):
                print("Data type: ", self.__labels[i], " Model MSE is: ",
                      self.__data_loss(data, type))

    def multiple_greatest_n_loss(self, label):
        for i, data in enumerate(self.__data_list):
            print("Data type: ", self.__labels[i], " Model MSE is: ",
                  self.__greatest_n_loss(data, label))

    def __greatest_n_loss(self, data, label):
        for n in MOST_COMMON_NUM:
            # extract the n most common mhcs
            mhcs = data[MHC_COL].to_list()
            most_common = [mhc for mhc, count in Counter(mhcs).most_common(n)]
            filtered_data = data[data[MHC_COL].isin(most_common)]
            predicts, true = self.__normalize_values(np.array(filtered_data[PREDICT_COL])), \
                             self.__normalize_values(np.array(filtered_data[TRUE_COL]))
            loss = metrics.mean_squared_error(true, predicts)
            print("Data type:", label, "Model G", n, "-MSE is:", loss)


class AlleleSpecificError:
    def __init__(self, labels):
        self.__labels = labels
        self.__trains, self.__tests, self.__results = [], [], []
        for lbl in labels:
            self.__trains.append(pd.read_csv('C:\\Users\\Yasmin\\PycharmProjects\\MHCBinding\\'
                                             + 'train_results_' + lbl + '.csv'))
            self.__tests.append(pd.read_csv('C:\\Users\\Yasmin\\PycharmProjects\\MHCBinding\\'
                                             + 'test_results_' + lbl + '.csv'))

    def create_errors(self, to_csv):
        for i in range(len(self.__labels)):
            print("creating allele error for label", self.__labels[i])
            self.__create_lbl_error(self.__trains[i], self.__tests[i], self.__labels[i], to_csv)

    def __create_lbl_error(self, train, test, label, to_csv):
        distinct_mhcs = list(set(train[MHC_COL].tolist()))
        res = []
        hist_res = [[], [], []]
        for allele in distinct_mhcs:
            tr_size, train_error = self.__calc_allele_error(train, allele)
            ts_size, test_error = self.__calc_allele_error(test, allele)
            res.append([allele, tr_size, train_error, ts_size, test_error])
            hist_res[ALLELE_INDEX].append(allele)
            hist_res[TRAIN_INDEX].append(train_error)
            hist_res[TEST_INDEX].append(test_error)
        a, tr, ts = np.array(hist_res[ALLELE_INDEX]), np.array(hist_res[TRAIN_INDEX]), \
                    np.array(hist_res[TEST_INDEX])
        argsort_indices = ts.argsort()
        a, tr, ts = a[argsort_indices], tr[argsort_indices], ts[argsort_indices]
        self.__results.append([a, tr, ts])
        if to_csv:
            df = pd.DataFrame(res, columns=['allele', 'train size', 'train error',
                                                    'test size', 'test error'])
            df.to_csv(label + "_allele_specific_errors.csv")

    def __calc_allele_error(self, data, allele):
        allele_data = data[data[MHC_COL] == allele]
        predictions, trues = np.array(allele_data[PREDICT_COL]), np.array(allele_data[TRUE_COL])
        loss = metrics.mean_squared_error(self.__normalize_values(trues),
                                          self.__normalize_values(predictions))
        return len(allele_data), loss

    def __normalize_values(self, l1):
        return list(map(lambda val: (np.log(val + 1) / NORMALIZATION), l1))

    def create_alleles_hist(self, threshold, display):
        loops = int(1 / threshold)
        for i, res in enumerate(self.__results):
            div_size = int(threshold * len(res[ALLELE_INDEX]))
            if display == 1:
                self.__manage_reg(loops, div_size, i, res)
            else:
                self.__manage_multiple_display(loops, div_size, i, res, display)

    def __manage_reg(self, loops, div_size, i, res):
        start_div, end_div = 0, div_size
        # divide the output into mini-graphs for better representation
        for j in range(loops):
            self.create_alleles_hist_with_threshold(res[ALLELE_INDEX][start_div:end_div],
                                                    res[TRAIN_INDEX][start_div:end_div],
                                                    res[TEST_INDEX][start_div:end_div], i, j)
            start_div += div_size
            end_div += div_size
            if end_div > len(res[ALLELE_INDEX]):
                end_div = len(res[ALLELE_INDEX])

    def __manage_multiple_display(self, loops, div_size, i, res, display_factor):
        loops = int(loops / display_factor) + 1
        start_div, end_div = 0, div_size
        # divide the output into mini-graphs for better representation
        for j in range(loops):
            data = []
            for k in range(display_factor):
                data.append([res[ALLELE_INDEX][start_div:end_div],
                             res[TRAIN_INDEX][start_div:end_div],
                             res[TEST_INDEX][start_div:end_div]])
                start_div += div_size
                end_div += div_size
                if end_div > len(res[ALLELE_INDEX]):
                    end_div = len(res[ALLELE_INDEX])
                    break
            if display_factor == 4:
                self.create_alleles_hist_of_four(data, i, j)
            elif display_factor == 2:
                self.create_alleles_hist_of_two(data, i, j)

    def create_alleles_hist_with_threshold(self, alleles, trains, tests, i, fig_num):
        f, ax = plt.subplots(1, 1)
        x = np.arange(len(alleles))
        y1, y2 = trains, tests
        width = 0.25
        ax.bar(x, y1, width, color='pink')
        ax.bar(x + width, y2, width, color='navy')
        ax.set_xticks(x + width)
        ax.set_xticklabels(alleles)
        ax.legend(['Train', 'Test'])
        fig_num = str(fig_num)
        plt.savefig(self.__labels[i] + "_allele_histogram_results_single" + fig_num + ".png")

    def create_alleles_hist_of_four(self, data, i, fig_num):
        fig, axes = plt.subplots(ncols=2, nrows=2)
        ax1, ax2, ax3, ax4 = axes.ravel()
        axes = [ax1, ax2, ax3, ax4]
        x = np.arange(len(data[0][ALLELE_INDEX]))
        width = 0.25

        for j, ax in enumerate(axes):
            if j <= len(data) - 1:  # valid index
                print("axe", j)
                y1, y2 = data[j][TRAIN_INDEX], data[j][TEST_INDEX]
                ax.bar(x, y1, width, color='pink')
                ax.bar(x + width, y2, width, color='navy')
                ax.set_xticks(x + width)
                ax.set_xticklabels(x)
                ax.legend(['Train', 'Test'])

        fig_num = str(fig_num)
        plt.savefig(self.__labels[i] + "_allele_histogram_results_4display_" + fig_num + ".png")

    def create_alleles_hist_of_two(self, data, i, fig_num):
        fig, axes = plt.subplots(ncols=1, nrows=2)
        ax1, ax2 = axes.ravel()
        axes = [ax1, ax2]
        x = np.arange(len(data[0][ALLELE_INDEX]))
        width = 0.25

        for j, ax in enumerate(axes):
            if j <= len(data) - 1:  # valid index
                y1, y2 = data[j][TRAIN_INDEX], data[j][TEST_INDEX]
                ax.bar(x, y1, width, color='pink')
                ax.bar(x + width, y2, width, color='navy')
                ax.set_xticks(x + width)
                ax.set_xticklabels(data[j][ALLELE_INDEX])
                ax.legend(['Train', 'Test'])

        fig_num = str(fig_num)
        plt.savefig(self.__labels[i] + "_allele_histogram_results_2display_" + fig_num + ".png")


def run_calculator(labels, results_label):
    basic_path = 'C:\\Users\\Yasmin\\PycharmProjects\\MHCBinding\\' + results_label + '_results_'
    paths = [basic_path + name + '.csv' for name in labels]
    calculator = Error_Calc(paths, labels)
    error_types = ['Pearson', 'MSE']
    for t in error_types:
        calculator.loss(t)


# data
full_ls = ['FULL', 'HIGHEST10', 'HIGHEST20', 'YEAR1995', 'YEAR2000', 'EC50_FULL', 'MAJORITY_VOTE',
           'YEAR', 'EC50']
labels = full_ls
RUN = False
if RUN:
    print("TRAIN results")
    run_calculator(labels, "train")
    print("TEST results")
    run_calculator(labels, "test")

HIST = False
EXTRACT_CSV = False
THRESHOLD = 0.05
DISPLAY = 2
hist_labels = ['FULL']
if HIST:
    allele_errors = AlleleSpecificError(hist_labels)
    allele_errors.create_errors(EXTRACT_CSV)
    allele_errors.create_alleles_hist(THRESHOLD, DISPLAY)

allele_labels = ['FULL']
ALLELES_NUM = True
if ALLELES_NUM:
    for t in ["train", "test"]:
        print(t, "results")
        basic_path = 'C:\\Users\\Yasmin\\PycharmProjects\\MHCBinding\\' \
                     + "FINAL_" + t + "_40_allele__results_"
        paths = [basic_path + name + '.csv' for name in allele_labels]
        calculator = Error_Calc(paths, labels)
        for l in allele_labels:
            calculator.multiple_greatest_n_loss(l)
# AUC
# predicts, true = np.array(data[PREDICT_COL]), \
#                  np.array(data[TRUE_COL])
# # predicts, true = np.array(self.__normalize_values(data[PREDICT_COL])), \
# #                  np.array(self.__normalize_values(data[TRUE_COL]))
# # predicts[predicts > 1] = 1
# # predicts[predicts < 0] = 0
# # true[true > 1] = 1
# # true[true < 0] = 0
# # fpr, tpr, thresholds = metrics.roc_curve(true, predicts, pos_label=2)
# # return metrics.auc(fpr, tpr)
# return stats.pearsonr(predicts, true)

# def __data_loss(self, data):
#     predicts, true = np.array(data[PREDICT_COL]), np.array(data[TRUE_COL])
#     fpr, tpr, thresholds = metrics.roc_curve(true, predicts, pos_label=2)
#     return metrics.auc(fpr, tpr)
    # todo old
    # true[true == 0] = MINIMAL_VALUE
    # predicts[predicts == 0] = MINIMAL_VALUE
    # normalized_predicts = predicts / true
    # normalized_predicts = normalized_predicts
    # normalized_true = np.full((1, np.alen(true)), 1.0, dtype=float)
    # return np.sum(np.power(normalized_true - normalized_predicts, 2)) / (np.alen(normalized_true) - 1)
    # # todo: create data sets with the inequalities
# for elemnt in predicts:
#     if elemnt == 0:
#         print(elemnt)
#     elif elemnt == MINIMAL_VALUE:
#         print((elemnt))
# for elemnt in true:
#     if elemnt == 0:
#         print(elemnt)
#     elif elemnt == MINIMAL_VALUE:
#         print((elemnt))
