import pandas as pd
import numpy as np
from scipy.io import loadmat
from math import isnan

from config import data_loc, results_loc
from utils import get_structured_array
from sklearn.preprocessing import OneHotEncoder


class Dataset:
    '''
    This class is to load the dataset pertaining to a particular split
    for running survival prediction
    '''
    def __init__(self, num_genes, fold_num, suffix, data_loc=data_loc, mode='encoded'):
        '''
        Save the configuration of the dataset
        '''
        self.num_genes = num_genes
        self.fold_num = fold_num
        self.data_loc = data_loc
        self.mode = mode
        
        if len(suffix) != 0:
            suffix = '_' + suffix
        self.suffix = suffix
        
        self.train_loc = data_loc + 'train' + self.suffix + '/' + str(num_genes) + "/" + str(fold_num)
        self.valid_loc = data_loc + 'valid' + self.suffix + '/' + str(num_genes) + "/" + str(fold_num)
        self.test_loc = data_loc + 'test' + self.suffix + '/' + str(num_genes) + "/" + str(fold_num)
        
        self.train_pca_loc =  data_loc + 'train_pca/100/' + str(fold_num)
        self.valid_pca_loc =  data_loc + 'valid_pca/100/' + str(fold_num)
        self.test_pca_loc =  data_loc + 'test_pca/100/' + str(fold_num)
        
        self.load_patients_data()
        self.load_genomics_data()
        self.load_PCA_genomics_data()
        self.load_imaging_data()
        self.load_PCA_imaging_data()
        self.load_PCA_concat_data()
#         self.load_concatenated_data()
#         self.load_clinical_data()
#         self.load_clinical_data_enc()
        
        return
    
#     def load_clinical_data(self):
# #         import ipdb
# #         ipdb.set_trace()
        
#         clinical_df = pd.read_csv(data_loc + 'clinical_patient_brca_subset.csv', sep='\t')
#         clinical_df = clinical_df.replace("[Not Available]", "0")
#         clinical_df = clinical_df.replace("not amplified", "0")
#         clinical_df = clinical_df.replace("<4.0", "4.0")
#         clinical_df = clinical_df.replace(">6", "6.0")
#         clinical_df = clinical_df.replace(">6.0", "6.0")
#         clinical_df = clinical_df.replace("polisomy", "0")
#         new_pids = []
        
#         for x in clinical_df['bcr_patient_barcode'].values:
#             new_pids.append('-'.join(x.split('.')))
            
# #         ipdb.set_trace()
#         clinical_df['PatientID'] = new_pids
#         clinical_df = clinical_df.set_index('PatientID')

#         columns_clinical = clinical_df.columns.values
# #         print(self.clinical_dict.head())
        
#         skip_columns = ['bcr_patient_barcode', 'PatientID']
#         numeric_columns = ['birth_days_to', 'age_at_diagnosis', 'her2_copy_number', 'cent17_copy_number',
#                            'her2_and_cent17_cells_count', 'her2_cent17_ratio']
        
#         cat_clinical_df = pd.DataFrame()
#         cat_clinical_df['PatientID'] = clinical_df.index.values
#         cat_clinical_df = cat_clinical_df.set_index('PatientID')
        
# #         ipdb.set_trace()
#         for column_name in columns_clinical:
            
#             if column_name in skip_columns:
#                 continue
#             elif column_name in numeric_columns:
#                 cat_clinical_df[column_name] = clinical_df[column_name]
#                 continue
                
#             clinical_df[column_name] = clinical_df[column_name].astype('category')
#             temp_column = clinical_df[column_name].cat.codes
#             if temp_column.std() < 1:
#                 continue
#             cat_clinical_df[column_name + '_cat'] = clinical_df[column_name].cat.codes
            
#         self.clinical_train = []
#         for x in self.train_patients:
#             self.clinical_train.append(cat_clinical_df.loc[x])
#         self.clinical_train = np.array(self.clinical_train)
        
#         self.clinical_valid = []
#         for x in self.valid_patients:
#             self.clinical_valid.append(cat_clinical_df.loc[x])
#         self.clinical_valid = np.array(self.clinical_valid)
            
#         self.clinical_test = []
#         for x in self.test_patients:
#             self.clinical_test.append(cat_clinical_df.loc[x])
#         self.clinical_test = np.array(self.clinical_test)
        
#     def load_clinical_data_enc(self):
        
#         self.clinical_dict = pd.read_csv(data_loc + 'clinical_subset.csv')
        
#         new_pids = []
#         for x in self.clinical_dict['PatientID'].values:
#             new_pids.append('-'.join(x.split('.')))
#         self.clinical_dict['PatientID'] = new_pids

#         if self.mode == 'encoded':
#             ERvalue_map = {'Positive': np.array([0, 1], dtype=int),
#                            'Negative': np.array([1, 0], dtype=int),
#                            'Indeterminate':np.array([1,1], dtype=int),
#                            'N/A': np.array([1,1], dtype=int)}
#         else:
#             ERvalue_map = {'Positive': [3], 'Negative':[-3], 'Indeterminate': [0], 'N/A': [0]}

#         ER_dict = self.clinical_dict[['PatientID', 'ER_STATUS_BY_IHC']].copy()
#         ER_dict = ER_dict.set_index('PatientID')
#         self.ER_dict = {}
#         for ind in ER_dict.index:
#             self.ER_dict[ind] = ER_dict.loc[ind][0]
#             if type(ER_dict.loc[ind][0]) == float and isnan(ER_dict.loc[ind][0]):
#                 self.ER_dict[ind] = "N/A"
#             else:
#                 self.ER_dict[ind] = self.ER_dict[ind]
                


#         self.ER_train = []
#         for x in self.train_patients:
#             if x + '-01' in self.ER_dict.keys():
# #                 self.ER_train.append(self.ER_dict[x + '-01'])
#                 self.ER_train.append(ERvalue_map[self.ER_dict[x + '-01']])
#             else:
# #                 self.ER_train.append('N/A')   
#                 self.ER_train.append(ERvalue_map['N/A'])   

#         self.ER_valid = []
#         for x in self.valid_patients:
#             if x + '-01' in self.ER_dict.keys():
# #                 self.ER_valid.append(self.ER_dict[x + '-01'])
#                 self.ER_valid.append(ERvalue_map[self.ER_dict[x + '-01']])
#             else:
# #                 self.ER_valid.append('N/A')   
#                 self.ER_valid.append(ERvalue_map['N/A'])  

#         self.ER_test = []
#         for x in self.test_patients:
#             if x + '-01' in self.ER_dict.keys():
# #                 self.ER_test.append(self.ER_dict[x + '-01'])
#                 self.ER_test.append(ERvalue_map[self.ER_dict[x + '-01']])
#             else:
# #                 self.ER_test.append('N/A')   
#                 self.ER_test.append(ERvalue_map['N/A'])  

#         PR_dict = self.clinical_dict[['PatientID', 'PR_STATUS_BY_IHC']].copy()
#         PR_dict = PR_dict.set_index('PatientID')

#         if self.mode == 'encoded':
#             PRvalue_map = {'Positive': np.array([0, 1], dtype=int),
#                            'Negative': np.array([1, 0], dtype=int),
#                            'Indeterminate':np.array([1,1], dtype=int),
#                            'N/A': np.array([1,1], dtype=int)}
#         else:
#             PRvalue_map = {'Positive': [3], 'Negative':[-3], 'Indeterminate': [0], 'N/A': [0]}

#         self.PR_dict = {}
#         for ind in PR_dict.index:
#             self.PR_dict[ind] = PR_dict.loc[ind][0]
#             if type(PR_dict.loc[ind][0]) == float and isnan(PR_dict.loc[ind][0]):
#                 self.PR_dict[ind]  = "N/A"
#             else:
#                 self.PR_dict[ind] = self.PR_dict[ind]

#         self.PR_train = []
#         for x in self.train_patients:
#             if x + '-01' in self.PR_dict.keys():
# #                 self.PR_train.append(self.PR_dict[x + '-01'])
#                 self.PR_train.append(PRvalue_map[self.PR_dict[x + '-01']])
#             else:
# #                 self.PR_train.append('N/A')
#                 self.PR_train.append(PRvalue_map['N/A'])   

#         self.PR_valid = []
#         for x in self.valid_patients:
#             if x + '-01' in self.PR_dict.keys():
# #                 self.PR_valid.append(self.PR_dict[x + '-01'])
#                 self.PR_valid.append(PRvalue_map[self.PR_dict[x + '-01']])
#             else:
# #                 self.PR_valid.append('N/A')
#                 self.PR_valid.append(PRvalue_map['N/A'])  

#         self.PR_test = []
#         for x in self.test_patients:
#             if x + '-01' in self.PR_dict.keys():
# #                 self.PR_test.append(self.PR_dict[x + '-01'])
#                 self.PR_test.append(PRvalue_map[self.PR_dict[x + '-01']])
#             else:
# #                 self.PR_test.append('N/A')
#                 self.PR_test.append(PRvalue_map['N/A'])   

#         HER2_dict = self.clinical_dict[['PatientID', 'HER2_FISH_STATUS', 'IHC_HER2']].copy()
#         HER2_dict = HER2_dict.set_index('PatientID')

#         if self.mode == 'encoded':
#             HER2value_map = {'Positive': np.array([0, 1], dtype=int),
#                        'Negative': np.array([1, 0], dtype=int),
#                        'N/A': np.array([1,1], dtype=int)}
#         else:
#             HER2value_map = {'Positive': [3], 'Negative':[-3],  'N/A': [0]}

#         self.HER2_dict = {}
#         for ind in HER2_dict.index:
#             value = HER2_dict.loc[ind]
#             if type(value[0]) == float and isnan(value[0]):
#                 new_val = value[1]
#             else:
#                 new_val = value[0]
#             if type(new_val) == str:
#                 if new_val in ['Equivocal', 'Indeterminate']:
#                     new_val = "N/A"
#             if type(new_val) == float and isnan(new_val):
#                 new_val  = "N/A"
#             self.HER2_dict[ind] = new_val

#         self.HER2_train = []
#         for x in self.train_patients:
#             if x + '-01' in self.HER2_dict.keys():
# #                 self.HER2_train.append(self.HER2_dict[x + '-01'])
#                 self.HER2_train.append(HER2value_map[self.HER2_dict[x + '-01']])
#             else:
# #                 self.HER2_train.append('N/A') 
#                 self.HER2_train.append(HER2value_map['N/A']) 

#         self.HER2_valid = []
#         for x in self.valid_patients:
#             if x + '-01' in self.HER2_dict.keys():
# #                 self.HER2_valid.append(self.HER2_dict[x + '-01'])
#                 self.HER2_valid.append(HER2value_map[self.HER2_dict[x + '-01']])
#             else:
# #                 self.HER2_valid.append('N/A')
#                 self.HER2_valid.append(HER2value_map['N/A'])

#         self.HER2_test = []
#         for x in self.test_patients:
#             if x + '-01' in self.HER2_dict.keys():
# #                 self.HER2_test.append(self.HER2_dict[x + '-01'])
#                 self.HER2_test.append(HER2value_map[self.HER2_dict[x + '-01']])
#             else:
# #                 self.HER2_test.append('N/A') 
#                 self.HER2_test.append(HER2value_map['N/A']) 

#         if self.mode == 'encoded':
#             PAM50value_map = {'Normal': np.array([0, 0, 1], dtype=int),
#                               'LumA': np.array([0, 1, 0], dtype=int),
#                               'LumB': np.array([1, 1, 0], dtype=int),
#                               'Basal': np.array([0, 1, 1], dtype=int),
#                               'Her2': np.array([1, 1, 1], dtype=int),
#                               'N/A': np.array([0, 0, 0], dtype=int)}
#         else:
#             PAM50value_map = {'Normal': [0], 'LumA': [1], 'LumB': [2],
#                               'Basal': [3], 'Her2':[4], 'N/A': [5]}


#         PAM50_dict = pd.read_csv(data_loc + 'PAM50.csv').to_dict()['x']
#         self.PAM50_dict = {}
#         for ind in PAM50_dict.keys():
#             new_ind = ind.strip() + '-01'
#             self.PAM50_dict[new_ind] = PAM50_dict[ind].strip()

#         self.PAM50_train = []
#         for x in self.train_patients:
#             if x + '-01' in self.PAM50_dict.keys():
# #                 self.PAM50_train.append(self.PAM50_dict[x + '-01'])
#                 self.PAM50_train.append(PAM50value_map[self.PAM50_dict[x + '-01']])
#             else:
# #                 self.PAM50_train.append('N/A')
#                 self.PAM50_train.append(PAM50value_map['N/A'])

#         self.PAM50_valid = []
#         for x in self.valid_patients:
#             if x + '-01' in self.PAM50_dict.keys():
# #                 self.PAM50_valid.append(self.PAM50_dict[x + '-01'])
#                 self.PAM50_valid.append(PAM50value_map[self.PAM50_dict[x + '-01']])
#             else:
# #                 self.PAM50_valid.append('N/A')
#                 self.PAM50_valid.append(PAM50value_map['N/A'])

#         self.PAM50_test = []
#         for x in self.test_patients:
#             if x + '-01' in self.PAM50_dict.keys():
# #                 self.PAM50_test.append(self.PAM50_dict[x + '-01'])
#                 self.PAM50_test.append(PAM50value_map[self.PAM50_dict[x + '-01']])
#             else:
# #                 self.PAM50_test.append('N/A')
#                 self.PAM50_test.append(PAM50value_map['N/A'])

#         self.clinical_train = []
#         for x in range(len(self.train_patients)):
#             self.clinical_train.append([list(self.ER_train[x]) +
#                                    list(self.PR_train[x]) +
#                                    list(self.HER2_train[x]) +
#                                    list(self.PAM50_train[x])])
#         self.clinical_train = np.array(self.clinical_train)

#         self.clinical_valid = []
#         for x in range(len(self.valid_patients)):
#             self.clinical_valid.append([list(self.ER_valid[x]) +
#                                    list(self.PR_valid[x]) +
#                                    list(self.HER2_valid[x]) + 
#                                    list(self.PAM50_valid[x])])
#         self.clinical_valid = np.array(self.clinical_valid)

#         self.clinical_test = []
#         for x in range(len(self.test_patients)):
#             self.clinical_test.append([list(self.ER_test[x]) +
#                                    list(self.PR_test[x]) +
#                                    list(self.HER2_test[x]) +
#                                    list(self.PAM50_test[x])])
#         self.clinical_test = np.array(self.clinical_test)

#         return
    
    def load_patients_data(self):
        '''
        Load the splits and imaging-genomics features within the
        class object
        '''
        
        # Load the patients data
        self.train_patients = list(pd.read_csv(
            self.train_loc + '/patients.csv', header=None)[0])
#         self.train_patients = ['-'.join(p.split('-')[:-1]) for p in self.train_patients]

        self.valid_patients = list(pd.read_csv(
            self.valid_loc + '/patients.csv', header=None)[0])
#         self.valid_patients = ['-'.join(p.split('-')[:-1]) for p in self.valid_patients]

        self.test_patients = list(pd.read_csv(
            self.test_loc + '/patients.csv', header=None)[0])
#         self.test_patients = ['-'.join(p.split('-')[:-1]) for p in self.test_patients]

        # Load survival times
        self.train_times = list(pd.read_csv(
            self.train_loc + '/times.csv', header=None)[0])
        self.valid_times = list(pd.read_csv(
            self.valid_loc + '/times.csv', header=None)[0])
        self.test_times = list(pd.read_csv(
            self.test_loc + '/times.csv', header=None)[0])

        # Load survival events
        self.train_events = list(pd.read_csv(
            self.train_loc + '/events.csv', header=None)[0])
        self.valid_events = list(pd.read_csv(
            self.valid_loc + '/events.csv', header=None)[0])
        self.test_events = list(pd.read_csv(
            self.test_loc + '/events.csv', header=None)[0])

    # train_times.extend(valid_times)
    # train_events.extend(valid_events)

        print("Length of Train", len(self.train_events))
        print("Length of Valid", len(self.valid_events))
        print("Length of Test",  len(self.test_events))

        print("Number of events in Train", sum(self.train_events))
        print("Number of events in Valid", sum(self.valid_events))
        print("Number of events in Test", sum(self.test_events))

        self.y_train = get_structured_array(
                self.train_events, self.train_times)
        self.y_valid = get_structured_array(
                self.valid_events, self.valid_times)
        self.y_test  = get_structured_array(
                self.test_events, self.test_times)


        return

    def load_genomics_data(self):
        '''
        Load genomics and imaging features
        '''
        
        self.genomics_train = np.array(pd.read_csv(
            self.train_loc + '/genes.csv', header=None));
        self.genomics_train = np.asarray(
            self.genomics_train, np.float64)

        self.genomics_valid = np.array(pd.read_csv(
            self.valid_loc + '/genes.csv', header=None));
        self.genomics_valid = np.asarray(
            self.genomics_valid, np.float64)

        # genomics_train = np.concatenate((genomics_train, genomics_valid), axis=0)

        self.genomics_test = np.array(pd.read_csv(
            self.test_loc + '/genes.csv', header=None));
        self.genomics_test = np.asarray(
            self.genomics_test, np.float64)

        return
        
    def load_PCA_genomics_data(self):
        '''
        Load genomics and imaging features
        '''
        
        self.PCA_genomics_train = np.array(pd.read_csv(
            self.train_pca_loc + '/genes.csv', header=None));
        self.PCA_genomics_train = np.asarray(
            self.PCA_genomics_train, np.float64)

        self.PCA_genomics_valid = np.array(pd.read_csv(
            self.valid_pca_loc + '/genes.csv', header=None));
        self.PCA_genomics_valid = np.asarray(
            self.PCA_genomics_valid, np.float64)

        # genomics_train = np.concatenate((genomics_train, genomics_valid), axis=0)

        self.PCA_genomics_test = np.array(pd.read_csv(
            self.test_pca_loc + '/genes.csv', header=None));
        self.PCA_genomics_test = np.asarray(
            self.PCA_genomics_test, np.float64)

        return

    def load_imaging_data(self):
        '''
        Load imaging features
        '''
        
        self.imaging_train = {}
        self.imaging_valid = {}
        self.imaging_test  = {}
                
        self.imaging_train = np.array(pd.read_csv(
            self.train_loc +'/image.csv', header=None));
        self.imaging_train = np.asarray(self.imaging_train, np.float64)

        self.imaging_valid = np.array(pd.read_csv(
            self.valid_loc + '/image.csv', header=None));
        self.imaging_valid = np.asarray(self.imaging_valid, np.float64)

        # imaging_train = np.concatenate((imaging_train, imaging_valid), axis=0)
        self.imaging_test = np.array(pd.read_csv(
            self.test_loc + '/image.csv', header=None));
        self.imaging_test = np.asarray(self.imaging_test, np.float64)

        return


    def load_PCA_imaging_data(self):
        '''
        Load imaging features
        '''
        
        self.PCA_imaging_train = {}
        self.PCA_imaging_valid = {}
        self.PCA_imaging_test  = {}
                
        self.PCA_imaging_train = np.array(pd.read_csv(
            self.train_loc + '/PCA_image.csv', header=None));
        self.PCA_imaging_train = np.asarray(self.PCA_imaging_train, np.float64)

        self.PCA_imaging_valid = np.array(pd.read_csv(
            self.valid_loc + '/PCA_image.csv', header=None));
        self.PCA_imaging_valid = np.asarray(self.PCA_imaging_valid, np.float64)

        # imaging_train = np.concatenate((imaging_train, imaging_valid), axis=0)
        self.PCA_imaging_test = np.array(pd.read_csv(
            self.test_loc + '/PCA_image.csv', header=None));
        self.PCA_imaging_test = np.asarray(self.PCA_imaging_test, np.float64)

        return



    def load_PCA_concat_data(self):
        '''
        Load concatenated PCA features
        '''

        self.PCA_concat_train = np.concatenate(
                        (self.PCA_genomics_train, self.PCA_imaging_train), axis=1)
        self.PCA_concat_valid = np.concatenate(
                        (self.PCA_genomics_valid, self.PCA_imaging_valid), axis=1)
        self.PCA_concat_test  = np.concatenate(
                        (self.PCA_genomics_test, self.PCA_imaging_test), axis=1)

        return

#     def load_concatenated_data(self):
#         '''
#         Concatenate the imaging and genomics data for the
#         early fusion model
#         '''
#         self.concat_train = {}
#         self.concat_valid = {}
#         self.concat_test  = {}

#         for subtype in ['Normal', 'Her2', 'LumA', 'LumB', 'Basal']:
#             self.concat_train[subtype] = np.concatenate(
#                     (self.genomics_train[subtype], self.imaging_train[subtype]),
#                      axis=1)
#             self.concat_valid[subtype] = np.concatenate(
#                     (self.genomics_valid[subtype], self.imaging_valid[subtype]),
#                      axis=1)
#             self.concat_test[subtype]  = np.concatenate(
#                     (self.genomics_test[subtype], self.imaging_test[subtype]),
#                      axis=1)

#         return

    def feature_generator(self, mode_flag, prior_flag, deflation_flag):
        '''
        The codes to generate features depending on the setting
        are all in this function. 
        '''
       
        mode_list = ['genomics', 'imaging', 'early-fusion', 'cca',
                     'scca', 'gnscca', 'gcca']
        
#         mode_list = ['genomics', 'imaging', 'early-fusion', 'clinical', 
#                      'late-fusion', 'scca', 'gn-scca', 'gcca']
       
        print(mode_flag, prior_flag, deflation_flag)
        if mode_flag not in mode_list:
            raise(NotImplementedError)
            
        labels = {}
        labels['train'] = self.y_train
        labels['valid'] = self.y_valid
        labels['test'] = self.y_test

        features = {}
        if mode_flag == 'genomics':
            features['train'] = self.genomics_train
            features['valid'] = self.genomics_valid
            features['test']  = self.genomics_test
            return features, labels

        if mode_flag == 'imaging':
            features['train'] = self.imaging_train
            features['valid'] = self.imaging_valid
            features['test']  = self.imaging_test
            return features, labels

        if mode_flag == 'early-fusion':
            features['train'] = self.concat_train
            features['valid'] = self.concat_valid
            features['test']  = self.concat_test
            return features, labels
        
#         if mode_flag == 'clinical':
#             features['train'] = self.clinical_train
#             features['valid'] = self.clinical_valid
#             features['test'] = self.clinical_test
#             return features, labels
        
        if mode_flag == 'late-fusion':
            raise(NotImplementedError)

        if mode_flag == 'scca':

            if prior_flag not in [0, 1]:
                raise(NotImplementedError)
            
            if prior_flag == 1:
                raise(NotImplementedError)

            if prior_flag == 0: 

                if deflation_flag is None:
                    raise(NotImplementedError)

                values = loadmat(results_loc +
                        'sparse_' + deflation_flag + '_' +
                        str(self.num_genes) + 
                        '_' + str(self.fold_num) + self.suffix + '_best.mat')

                U = (np.squeeze(values['u']))
                V = (np.squeeze(values['v']))

                Xproject_train = np.dot(
                        self.genomics_train,U).astype(np.float64)
                Yproject_train = np.dot(
                        self.imaging_train,V).astype(np.float64)

                Xproject_valid = np.dot(
                        self.genomics_valid,U).astype(np.float64)
                Yproject_valid = np.dot(
                        self.imaging_valid,V).astype(np.float64)
                
                Xproject_test = np.dot(
                        self.genomics_test,U).astype(np.float64)
                Yproject_test = np.dot(
                        self.imaging_test,V).astype(np.float64)

                E_train = np.concatenate(
                        (Xproject_train, Yproject_train), axis=1)
                E_valid = np.concatenate(
                        (Xproject_valid, Yproject_valid), axis=1)
                E_test  = np.concatenate(
                        (Xproject_test, Yproject_test), axis=1)

                features['train'] = [Xproject_train, Yproject_train, E_train]
                features['valid'] = [Xproject_valid, Yproject_valid, E_valid]
                features['test']  = [Xproject_test, Yproject_test, E_test]
                return features, labels

        if mode_flag in ['cca', 'gnscca', 'gcca']:

            prior_flag = ''
            
            values = loadmat(results_loc +
                    mode_flag + prior_flag + '_' +
                    deflation_flag + '_' +
                    str(self.num_genes) + 
                    '_' + str(self.fold_num) + self.suffix +'_best.mat')

            U = (np.squeeze(values['u']))
            V = (np.squeeze(values['v']))

            Xproject_train = np.dot(
                    self.genomics_train,U).astype(np.float64)
            Yproject_train = np.dot(
                    self.imaging_train,V).astype(np.float64)

            Xproject_valid = np.dot(
                    self.genomics_valid,U).astype(np.float64)
            Yproject_valid = np.dot(
                    self.imaging_valid,V).astype(np.float64)
            
            Xproject_test = np.dot(
                    self.genomics_test,U).astype(np.float64)
            Yproject_test = np.dot(
                    self.imaging_test,V).astype(np.float64)

            E_train = np.concatenate(
                    (Xproject_train, Yproject_train), axis=1)
            E_valid = np.concatenate(
                    (Xproject_valid, Yproject_valid), axis=1)
            E_test  = np.concatenate(
                    (Xproject_test, Yproject_test), axis=1)

            features['train'] = [Xproject_train, Yproject_train, E_train]
            features['valid'] = [Xproject_valid, Yproject_valid, E_valid]
            features['test']  = [Xproject_test, Yproject_test, E_test]
            return features, labels

        return
