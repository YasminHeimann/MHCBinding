from mhcflurry import Class1AffinityPredictor, custom_loss
import pandas as pd, numpy as np
from sklearn import metrics
from modelError import run_calculator
from collections import Counter
from scipy import stats

EPITOPE_INDEX = 0
MHC_INDEX = 1
AFFINITY_INDEX = 2
TRAIN = 0
TEST = 1
MHC_COL = "allele"
EPITOPE_COL = "peptide"
AFFINITY_COL = "measurement_value"
INEQ_COL = 'measurement_inequality'
PREDICTION_INDEX = 0
NORMALIZATION = np.log(50000)
# indicates whether this is a test run (reads from excel and not csv)
SAMPLE = False

OUTPUT_COLS = ["allele","peptide","measurement_value","measurement_inequality",
               "measurement_type","measurement_source","original_allele"]


def create_alleles_histogram(train, test):
    train_distinct_mhcs = Counter(train[MHC_COL].to_list())
    test_mhcs = dict(Counter(test[MHC_COL].to_list()))  # TODO CHECK
    results = []
    col_names = ["allele", "peptides in train", "train error", "peptides in test", "test error"]
    for mhc, count in train_distinct_mhcs:
        train_allele_error, test_allele_error = 0, 0
        results.append([mhc, count, train_allele_error, test_mhcs[mhc], test_allele_error])
    # convert the list of lists into a df with the columns above
    # extract to csv
    # create histogram


def normalization(l1):
    return list(map(lambda val: (np.log(val + 1) / NORMALIZATION), l1))


def calc_test_mse(true, predicted, label):
    print(label, "MSE without normalization: ", metrics.mean_squared_error(true, predicted))
    n_t, n_p = normalization(true), normalization(predicted)
    print(label, "MSE with normalization: ", metrics.mean_squared_error(n_t,n_p))
    print(label, "Pearson with normalization: ", stats.pearsonr(n_p, n_t))


def calculate_model_error(predictor, test, label):
    print(label, " error results:")
    # get the predictions on the test set
    alleles, predicts, real_affinities = [], [], []
    valid_alleles = predictor.supported_alleles
    # iterate the rows in the test data frame to predict valid alleles
    for i, assay in test.iterrows():
        allele = assay[MHC_COL]
        if allele not in valid_alleles:
            continue
        else:
            predicts.append(predictor.predict([assay[EPITOPE_COL]], allele=allele)[0])
            real_affinities.append(assay[AFFINITY_COL])
            alleles.append(allele)
    predicts = np.array(predicts)
    real_affinities = np.array(real_affinities)
    alleles = np.array(alleles)

    # extract the test output to an excel
    df = pd.DataFrame({'mhc allele': alleles,'predicted labels': predicts, 'true labels': real_affinities})
    df.to_csv(label + "_results_" + LABEL + ".csv", index=False)
    # calculate the error
    calc_test_mse(real_affinities, predicts, label)
    return df


def run_network(allele_train_data, predictor_path, distinct_mhcs, n_model):
    new_predictor = Class1AffinityPredictor()
    # train all the alleles
    for i, allele in enumerate(distinct_mhcs):
        single_allele_train_data = allele_train_data.loc[allele_train_data.allele == allele]
        single_allele_train_data = single_allele_train_data.loc[(single_allele_train_data.peptide.str.len() >= 8) &
                                                                (single_allele_train_data.peptide.str.len() <= 15)]
        try:
            new_predictor.fit_allele_specific_predictors(
                n_models=n_model, architecture_hyperparameters_list=[{"layer_sizes": [16],
                                                                "max_epochs": 5, "random_negative_constant": 5}],
                peptides=single_allele_train_data.peptide.values,
                affinities=single_allele_train_data.measurement_value.values,
                allele=allele)  # , models_dir_for_save=predictor_path)
        except Exception as exp:
            print(exp, ". in allele: ", allele)

    # save the predictor
    # new_predictor.save(models_dir=predictor_path)
    return new_predictor


def run_all(test_path, train_path):
    # create df from the data sets paths
    if not SAMPLE:
        train = pd.read_csv(train_path)
        test = pd.read_csv(test_path)
    else:
        train = pd.read_excel(train_path)
        test = pd.read_excel(test_path)
    distinct_mhcs = list(set(train[MHC_COL].tolist()))

    predictor_path = "predictors/"
    predictor = run_network(train, predictor_path, distinct_mhcs, 1)
    # if FIRST_RUN:
    #     predictor = run_network(train, predictor_path, distinct_mhcs)
    # else:
    #     predictor = Class1AffinityPredictor.load(predictor_path)

    calculate_model_error(predictor, test, "test")
    calculate_model_error(predictor, train, "train")
    #create_alleles_histogram(final_train_df, final_test_df)


def run_most_common_alleles(test_path, train_path, allele_nums, n_model):
    # create df from the data sets paths
    train = pd.read_csv(train_path)
    test = pd.read_csv(test_path)
    for n in allele_nums:
        all_mhcs = train[MHC_COL].to_list()
        most_common_alleles = [mhc for mhc, count in Counter(all_mhcs).most_common(n)]
        filtered_train_data = train[train[MHC_COL].isin(most_common_alleles)]
        filtered_test_data = test[test[MHC_COL].isin(most_common_alleles)]

        predictor_path = "predictors/"
        predictor = run_network(filtered_train_data, predictor_path, most_common_alleles, n_model)

        test_name = "FINAL_test_" + str(n) + "_allele_"
        train_name = "FINAL_train_" + str(n) + "_allele_"
        calculate_model_error(predictor, filtered_test_data, test_name)
        calculate_model_error(predictor, filtered_train_data, train_name)


def test():
    label = "sample"
    test_path = f"test_train/test_iedb_data_{label}.xlsx"
    train_path = f"test_train/train_iedb_data_{label}.xlsx"
    print("Running MHCflurry network for label: ", label)
    run_all(test_path, train_path)


def main(filters, mode="most_common"):
    common_alleles_max = [40]
    model_num = 20
    for ft in filters:
        print("learning label ", l, "with", model_num, " models")
        test_path = f"test_train/test_iedb_data_{ft}.csv"
        train_path = f"test_train/train_iedb_data_{ft}.csv"
        if mode == "most_common":
            run_most_common_alleles(test_path, train_path, common_alleles_max, model_num)
        print("Running MHCflurry network for label: ", LABEL)
        run_all(test_path, train_path)


if __name__ == "___main__":
    if SAMPLE:
        test()
    else:
        types = ["FULL", "HIGHEST10", "HIGHEST20", "EC50_FULL", "YEAR2000", "YEAR1995", "MAJORITY_VOTE", "EC50", "YEAR"]
        # can send a list with only a single filter type
        main(types)