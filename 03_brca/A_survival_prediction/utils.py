'''
Contains the utility functions
'''

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

import sys
sys.path.append('../')

from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import brier_score
from sklearn.model_selection import GridSearchCV, KFold, PredefinedSplit
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage

from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter


def get_structured_array(data_bool, data_value):
    all_bools = data_bool
    all_values = data_value

    new_list = []
    for idx in range(len(all_bools)):
        new_list.append(tuple((all_bools[idx], all_values[idx])))
    return np.array(new_list, dtype='bool, i8')


def run_coxph(x_train, y_train, x_valid, y_valid,
        x_test, y_test):
    coxph = CoxPHFitter()

    train_df = pd.DataFrame(x_train)
    return

def run_coxnet(l1_ratio, n_alphas,
        x_train, y_train, 
        x_valid, y_valid,
        x_test, y_test, plot=False):

    alphas = 10. ** np.linspace(-3.5, 6, n_alphas)

    coxnet = make_pipeline(StandardScaler(),
            CoxnetSurvivalAnalysis(l1_ratio=l1_ratio,
                fit_baseline_model=True))
    
    if isinstance(x_train, np.ndarray) or isinstance(x_train, np.array):
        X_comb = np.concatenate((x_train, x_valid))
    else:
        X_comb = x_train + x_valid

    Y_comb = np.concatenate((y_train, y_valid))

    train_ind = np.arange(len(x_train))
    valid_ind = len(x_train) + np.arange(len(x_valid))

    my_cv = [[np.array(train_ind), np.array(valid_ind)]]

    gcv = GridSearchCV(
        make_pipeline(StandardScaler(),
            CoxnetSurvivalAnalysis(
                l1_ratio=l1_ratio, fit_baseline_model=True)),
        param_grid={
            "coxnetsurvivalanalysis__alphas": [[v] for v in alphas]},
        cv = my_cv,
        error_score=0.5,
        n_jobs=4).fit(X_comb, Y_comb)

    cv_results = pd.DataFrame(gcv.cv_results_)

    best_model = gcv.best_estimator_
    outputs = best_model.predict(x_test)

    c_idx_score = best_model.score(x_test, y_test)
    print("**************************")
    print(c_idx_score)
    
    return outputs, c_idx_score, best_model


def run_clustering(x_train, x_test, n_clusters=4):
    scaler = StandardScaler()
    scaler.fit(x_train)
    x_test_transformed = scaler.transform(x_test)
    linked = linkage(x_test_transformed, 'ward')

    return scaler, linked, x_test_transformed

def generate_km_curve(durations, events):
    kmf = KaplanMeierFitter()
    kmf.fit(durations, events)

    return kmf
