from utils import *
from config import save_loc
from dataloader import Dataset
import numpy as np
import matplotlib
matplotlib.use('agg')
# from matplotlib import pyplot as plt
# from math import isnan
# from scipy.cluster.hierarchy import fcluster
# from sklearn.cluster import AgglomerativeClustering
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.preprocessing import OneHotEncoder
# import seaborn as sns
# from collections import Counter

from lifelines import CoxPHFitter
coxph = CoxPHFitter(penalizer=0.1)

import numpy as np
import pandas as pd
# from matplotlib import pyplot as plt

import sys
sys.path.append('../')

# from sksurv.linear_model import CoxnetSurvivalAnalysis
# from sksurv.metrics import brier_score
# from sklearn.model_selection import GridSearchCV, KFold
# from sklearn.pipeline import make_pipeline
# from sklearn.preprocessing import StandardScaler


def coxph_settings(fold=0, mode='feature', method='scca',
                   deflation='opd', suffix ='log_histogram', 
                   prior=0, num_genes=1000):
    
    brca = Dataset(num_genes, fold, suffix)
    features, labels = brca.feature_generator(method, prior, deflation)
    
    events_train = [x for (x,y) in labels['train']]
    times_train = [y for (x,y) in labels['train']]

    events_test = [x for (x,y) in labels['test']]
    times_test = [y for (x,y) in labels['test']]

    if mode == 'feature':
        train_df = pd.DataFrame(features["train"][2])
        test_df = pd.DataFrame(features["test"][2])
        
    elif mode == 'PCAgenomics':
        train_df = pd.DataFrame(brca.PCA_genomics_train)
        test_df = pd.DataFrame(brca.PCA_genomics_test)
        
    elif mode == 'genomics':
        train_df = pd.DataFrame(features["train"][0])
        test_df = pd.DataFrame(features["test"][0])
        
    elif mode == 'PCAimaging':
        train_df = pd.DataFrame(brca.PCA_imaging_train)
        test_df = pd.DataFrame(brca.PCA_imaging_test)
        
    elif mode == 'imaging':
        train_df = pd.DataFrame(features["train"][1])
        test_df = pd.DataFrame(features["test"][1])

    elif mode == 'PCAconcat':
        train_df = pd.DataFrame(brca.PCA_concat_train)
        test_df = pd.DataFrame(brca.PCA_concat_test)
        
    else:
        raise(NotImplementedError)
              
    if train_df.isnull().values.any() or train_df.isnull().values.any():
        return None, None
    
    train_df['Events'] = pd.Series(np.array(events_train, dtype=int))
    train_df['Times'] = pd.Series(np.array(times_train, dtype=int))

    test_df['Events'] = pd.Series(np.array(events_test, dtype=int))
    test_df['Times'] = pd.Series(np.array(times_test, dtype=int))

    print(mode, train_df.shape, test_df.shape)
    coxph.fit(train_df, duration_col='Times', event_col='Events')    
    test_score = coxph.score(test_df, 'concordance_index')
    return test_score, coxph


def run_all_fold(setting, fold_range):
    
    
    setting_scores = []
    all_models = []
    for fold_num in fold_range:
        curr_run = coxph_settings(mode=setting['mode'], fold=fold_num, 
                                             suffix=setting['suffix'], 
                                             method=setting['method'],
                                             deflation=setting['deflation'])
        setting_scores.append(curr_run[0])
        all_models.append(curr_run[1])
    return setting_scores, all_models

if __name__=='__main__':
    
    fold_range = list(range(5))
    
    deflation = 'opd'
    suffix = 'log_histogram'
    PCAgenomics_setting = {'mode':'PCAgenomics', 'method':'scca', 'suffix':suffix, 'deflation':deflation}
    PCAgenomics_score, PCAgenomics_models = run_all_fold(PCAgenomics_setting, fold_range)
    np.savez(save_loc + 'PCAgenomics', PCAgenomics_score, PCAgenomics_models)
    print('PCAgenomics', PCAgenomics_score)

    suffix = 'log_histogram'
    PCAimaging_setting = {'mode':'PCAimaging', 'method':'scca', 'suffix':suffix, 'deflation':deflation}
    PCAimaging_score, PCAimaging_models = run_all_fold(PCAimaging_setting, fold_range)
    np.savez(save_loc + 'PCAimaging', PCAimaging_score, PCAimaging_models)
    print('PCAimaging', PCAimaging_score)

    deflation = 'opd'
    suffix = 'log_histogram'
    PCAconcat_setting = {'mode':'PCAconcat', 'method':'scca', 'suffix':suffix, 'deflation':deflation}
    PCAconcat_score, PCAconcat_models = run_all_fold(PCAconcat_setting, fold_range)
    np.savez(save_loc + 'PCAconcat', PCAconcat_score, PCAconcat_models)
    print('PCAconcat', PCAconcat_score)
    
    for suffix in ['log_histogram']:
        for method in ['scca', 'gnscca']: #, 'gcca']:
            for deflation in ['hd', 'pd', 'opd']:
                genomics_setting = {'mode':'genomics', 'method':method, 'suffix':suffix, 'deflation':deflation}
                genomics_score, genomics_models = run_all_fold(genomics_setting, fold_range)
                np.savez(save_loc + 'feature_' + method + '_' + 
                         suffix + '_' + deflation + '_genomics', genomics_score, genomics_models)
                print('genomics', genomics_score)
                suffix = 'log_histogram'
                imaging_setting = {'mode':'imaging', 'method':method, 'suffix':suffix, 'deflation':deflation}
                imaging_score, imaging_models = run_all_fold(imaging_setting, fold_range)
                np.savez(save_loc + 'feature_' + method + '_' + 
                         suffix + '_' + deflation + '_imaging', imaging_score, imaging_models)
                print('imaging', imaging_score)
              
                feature_setting_scca = {'mode':'feature', 'method':method, 'suffix':suffix, 'deflation':deflation}
                print(feature_setting_scca)
                feature_scca_score, scca_models = run_all_fold(feature_setting_scca, fold_range)
                np.savez(save_loc + 'feature_' + method + '_' + 
                         suffix + '_' + deflation, feature_scca_score, scca_models)
                print(suffix, method, deflation)
                print(feature_scca_score)
               