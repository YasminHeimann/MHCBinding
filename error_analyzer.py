import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np, pandas as pd
from collections import Counter, defaultdict

# whether to show data plot
SHOW_PLOT = False
# whether to create the data's errors' plots
DATA_PLOTS = False
# whether to extract a file with the highest errors
TOP_RESULTS_EXCEL = False
# check the john case
JOHN_ARTICLES = False
WITH_NAN = True

##
# values for the log ratio
##
NORMALIZATION = np.log(50000)
EPSILON = 1E-15
LOG_EPSILON = np.log(EPSILON)

ERROR_INDEX = 0
GROUP_INDEX = 1
VAL1 = 2
VAL2 = 3
TOP_ERRORS_NUM = 100


COL_ORDER = ['MHC ligand ID', 'Reference ID', 'Epitope ID', 'Description', 'Allele Name',
             'Units', 'Qualitative Measure', 'Measurement Inequality', 'Quantitative measurement',
             'PubMed ID', 'Date', 'Assay Comments', 'Type', 'Authors', 'Journal', 'Title', 'Submission ID',
             'Object Type', 'Starting Position', 'Ending Position', 'Non-peptidic epitope ChEBI ID',
             'Antigen Name', 'Parent Protein', 'Parent Protein Accession', 'Organism Name', 'Parent Species',
             'Parent Species Accession', 'Epitope Comments', 'Epitope Relationship', 'Object Type.1',
             'Description.1', 'Starting Position.1', 'Ending Position.1', 'Non-peptidic object Accession',
             'Synonyms', 'Antigen Name.1', 'Parent Protein.1', 'Organism Name.1', 'Parent Organism',
             'Name', 'Host ID', 'Geolocation', 'MHC Types Present', 'Process Type', 'Disease State',
             'Disease Stage', 'Processed Antigen Epitope Relation', 'Processed Antigen Object Type',
             'Processed Antigen Description', 'Processed Antigen Starting Position',
             'Processed Antigen Ending Position', 'Non-peptidic Processed Antigen ChEBI ID',
             'Processed Antigen Source Molecule Name', 'Processed Antigen protein parent Name',
             'Processed Antigen protein parent Accession', 'Processed Antigen Organism Name',
             'Processed Antigen Organism Species', 'Processed Antigen Organism Species ID',
             'In vitro administration type', 'Processed Antigen Epitope Relation.1',
             'Processed Antigen Object Type.1', 'Processed Antigen Description.1',
             'Processed Antigen Starting Position.1', 'Processed Antigen Ending Position.1',
             'Non-peptidic Processed Antigen ChEBI ID.1', 'Processed Antigen Source Molecule Name.1',
             'Protein Parent Name', 'Protein Parent Accession', 'Processed Antigen Organism Name.1',
             'Immunogen Organism Species', 'Immunogen Organism Species ID', 'Processed Antigen Comments',
             'Location of assay data in the manuscript', 'Method/Technique', 'Assay Group',
             'Number of Subjects Tested', 'Number of Subjects Responded', 'Response Frequency',
             'PDB ID', 'Cell Tissue Type', 'Cell Type', 'Cell Culture Conditions', 'Allele Evidence Code',
             'MHC allele class']
ERROR_COLS = ['Description', 'Allele Name', 'Measurement Inequality', 'Quantitative measurement',
              'PubMed ID', 'Date', 'Title', 'Method/Technique']
AFFINITY_COL = 'Quantitative measurement'
METHOD_COL = 'Method/Technique'

METHOD = 0
AFFINITIES = 1
M1 = 0
M2 = 1
M_ERR = 2
SINGLE_METHOD = False
DOUBLE_METHOD = True


def create_john_histogram(errors):
    """
    Saves a figure with a histogram indicating the distribution of errors that were caused due to
    an assay including John Sidney.
    Also plots an excel of the assays.
    :param params: includes a list of pairs (as list) of the assays' group and the measured error
    """
    f, ax = plt.subplots(1, 1)
    ax.hist(errors, bins=20)
    mean = '%.5f' % np.mean(errors)
    mean_msg = "mean of error: " + str(mean) + ". #assys: " + str(len(errors))
    anchored_text = AnchoredText(mean_msg, loc=2)
    ax.add_artist(anchored_text)
    plt.xlabel('amount of errors')
    plt.ylabel('error values')
    plt.title('A Histogram of the identical assays\' affinity values difference(John)')
    plt.savefig("errors_john_blanks.png", dpi=600)


def create_errors_histogram(errors, mean):
    f, ax = plt.subplots(1, 1)
    ax.hist(errors, bins=20)
    mean = '%.5f' % mean
    mean_msg = "mean of error: " + str(mean)
    anchored_text = AnchoredText(mean_msg, loc=2)
    ax.add_artist(anchored_text)
    plt.xlabel('amount of errors')
    plt.ylabel('error values')
    plt.title('A Histogram of the identical assays\' affinity values difference')
    if WITH_NAN:
        plt.savefig("errors_blanks.png", dpi=600)
    else:
        plt.savefig("errors_eq.png", dpi=600)
    if SHOW_PLOT:
        plt.show()


def create_pairs_affinity_plot(x_axis, y_axis):
    plt.xlabel('Affinity Values')
    plt.ylabel('Paired Values of Identical Assay')
    plt.title('Pairs of identical assays - showing the affinity values')
    plt.subplot(2,1,1)
    plt.plot(x_axis, y_axis,'.b')
    plt.plot(y_axis, x_axis, '.b')
    plt.subplot(2, 1, 2)
    plt.hist2d(x_axis+y_axis,y_axis+x_axis,bins=20, cmap='Blues')
    plt.legend()
    plt.show()
    if WITH_NAN:
        plt.savefig("affinities_blanks_hist.png", dpi=600)
    else:
        plt.savefig("affinities_eq.png",dpi=600)
    if SHOW_PLOT:
        plt.show()


def create_greatest_errors_csv(errors):
    if WITH_NAN:
        file_name = "greatest" + str(TOP_ERRORS_NUM) + "_errors_blanks.xlsx"
    else:
        file_name = "greatest" + str(TOP_ERRORS_NUM) + "_errors_eq.xlsx"
    row = 0
    writer = pd.ExcelWriter(file_name, engine="xlsxwriter")
    count = 1
    for err in errors:
        #if err is not None: -for sample
        text1 = "Error: " + str(count)
        text2 = str(err[ERROR_INDEX])
        text3 = " with values: "
        text4 = str(err[VAL1]) + " ," + str(err[VAL2])

        df = err[GROUP_INDEX]
        df = df[COL_ORDER]
        df.to_excel(writer, startrow=(row + 2), startcol=0, sheet_name='Sheet1')
        worksheet = writer.sheets['Sheet1']
        worksheet.write(row, 0, text1)
        worksheet.write(row, 1, text2)
        worksheet.write(row, 2, text3)
        worksheet.write(row, 3, text4)

        row = row + 4 + len(err[GROUP_INDEX])
        count += 1
    writer.save()


def create_percentages(errors):
    err = np.array(errors)
    print("over 0.2 errors: ", len(err[err >= 0.2]))
    print("over 0.4 errors: ", len(err[err >= 0.4]))
    print("over 0.6 errors: ", len(err[err >= 0.6]))


def analyze_error(groups):  # groupBy object
    # https://www.tutorialspoint.com/python_pandas/python_pandas_groupby.htm
    pairs_errors, pairs_names, pairs_x_vals, pairs_y_vals = [], [], [], []

    greatest_ten = [None]*TOP_ERRORS_NUM
    greatest_ten_errors = np.zeros(TOP_ERRORS_NUM)
    argmin_error = 0
    val_min_error = 1E-15
    count = 0

    for name, g in groups:
        # generates a list with the affinity values
        affinity_vals = g[AFFINITY_COL].tolist()
        affinity_vals.sort(reverse=True)  # descending order
        size = len(affinity_vals)
        # replaces 0 values and apply log scale
        for k in range(size):
            affinity_vals[k] = np.log(affinity_vals[k] + 1) / NORMALIZATION
        # go through all the pairs in the group, and calculate the error
        for i in range(size):
            val1 = affinity_vals[i]
            for j in range(i+1, size):
                val2 = affinity_vals[j]
                error = abs(val1 - val2)
                if error != -float("Inf"):
                    pairs_errors.append(error)
                else:
                    pairs_errors.append(LOG_EPSILON)
                pairs_names.append(name)
                pairs_x_vals.append(val1)
                pairs_y_vals.append(val2)
                # update the greatest errors
                if count < TOP_ERRORS_NUM:
                    greatest_ten[count] = (error, g.reset_index(), val1, val2)
                    greatest_ten_errors[count] = error
                    count += 1
                elif (error > val_min_error):
                    greatest_ten[argmin_error] = (error, g.reset_index(), val1, val2)
                    greatest_ten_errors[argmin_error] = error
                    argmin_error = np.argmin(greatest_ten_errors)
                    val_min_error = greatest_ten_errors[argmin_error]

    create_percentages(pairs_errors)
    if DATA_PLOTS:
        create_pairs_affinity_plot(pairs_x_vals, pairs_y_vals)
        create_errors_histogram(pairs_errors, np.mean(pairs_errors))
    #print(greatest_ten)
    if TOP_RESULTS_EXCEL:
        create_greatest_errors_csv(greatest_ten)


def create_method_errors_histogram(method_name, errors, mean):
    f, ax = plt.subplots(1, 1)
    ax.hist(errors, bins=20)
    mean = '%.5f' % mean
    mean_msg = "mean of error: " + str(mean) + ". number of errors: " + str(len(errors))
    anchored_text = AnchoredText(mean_msg, loc=2)
    ax.add_artist(anchored_text)
    plt.ylabel('amount of errors')
    plt.xlabel('error values')
    plt.title(method_name)
    # plt.show()
    method_name = method_name.replace("/", "_")
    method_name = method_name.replace(" ", "_")
    method_name = method_name.replace(",", "_")
    plt.savefig(method_name, dpi=600)


def normalize_values(l1):
    return list(map(lambda val: (np.log(val + 1) / NORMALIZATION), l1))


def calc_single_error(val1, val2):
    error = abs(val1 - val2)
    return error if error != -float("Inf") else LOG_EPSILON


def calc_two_methods_error(l1, l2):
    error_list = []
    for i in range(len(l1)):
        for j in range(len(l2)):
            error_list.append(calc_single_error(l1[i], l2[j]))
    return max(error_list)


def calc_dict_errors(dups_dict):
    loops = len(dups_dict) - 1
    all_errors = []
    for l in range(loops):
        cur_method_pair = dups_dict.popitem()  # notice: error if empty
        for temp_method, temp_affinities in dups_dict.items():
            all_errors.append((cur_method_pair[METHOD], temp_method, calc_two_methods_error(
                normalize_values(cur_method_pair[AFFINITIES]),
                normalize_values(temp_affinities))))
    return all_errors


def analyze_method_error(groups):
    all_methods_single = defaultdict(list)
    all_methods_double = MethodHandler()
    for name, g in groups:
        methods = g[METHOD_COL].tolist()
        affinities = g[AFFINITY_COL].tolist()
        # create a dictionary of the methods and the affinities
        dups_dict = defaultdict(list)
        for m, a in zip(methods, affinities):
            dups_dict[m].append(a)
        if len(dups_dict) > 1:
            # list of each methods' pairs and their highest error
            errors = calc_dict_errors(dups_dict)
            for i in range(len(errors)):
                if SINGLE_METHOD:
                    all_methods_single[errors[i][M1]].append(errors[i][M_ERR])
                    all_methods_single[errors[i][M2]].append(errors[i][M_ERR])
                if DOUBLE_METHOD:
                    all_methods_double.add_method_error(errors[i][M1], errors[i][M2], errors[i][M_ERR])

    # extract the errors into a graph
    if SINGLE_METHOD:
        for method, err_list in all_methods_single.items():
                create_method_errors_histogram(method, err_list, np.mean(err_list))
    if DOUBLE_METHOD:
        all_methods_double.get_final_results()


class MethodHandler:
    def __init__(self):
        self.__all_method_pairs = []

    def add_method_error(self, m1, m2, err):
        exist = False
        for pair in self.__all_method_pairs:
            if pair.equal(m1, m2):
                pair.add_error(err)
                exist = True
                break
        if not exist:
            new_pair = self.MethodsPair(m1, m2)
            new_pair.add_error(err)
            self.__all_method_pairs.append(new_pair)

    def get_final_results(self):
        for pair in self.__all_method_pairs:
            res = pair.get_single_pair_results()
            title = res[0] + ", " + res[1]
            create_method_errors_histogram(title, res[2], np.mean(res[2]))

    class MethodsPair():
        def __init__(self, m1, m2):
            self.__m1, self.__m2, self.__err_list = m1, m2, []

        def add_error(self, err):
            self.__err_list.append(err)

        def equal(self, other_m1, other_m2):
            if (self.__m1 == other_m1 and self.__m2 == other_m2) or (self.__m1 == other_m2 and self.__m2 == other_m1):
                return True
            return False

        def get_single_pair_results(self):
            return self.__m1, self.__m2, self.__err_list
