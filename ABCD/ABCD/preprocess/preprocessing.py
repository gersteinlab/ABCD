import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.impute import SimpleImputer
import os
from tqdm import tqdm
import re
import warnings
from sktime.transformations.series.impute import Imputer
from ABCD.preprocess.utils import processing_genomics
from ABCD.preprocess.utils.features import *
from ABCD.preprocess.utils.subsetting_days import *
from ABCD.preprocess.utils.normalize import normalize
warnings.filterwarnings('default')

def process_data(config):
    '''
    This function preprocess the raw wearables data into a normalized numpy dataframe.
    return -> numpy [N, C, L]
    N for # subjects
    C for # features
    L for # measurements
    '''
    if config['add_genomics']:
        genomics = processing_genomics()

    group = config['group']
    label_path = config['label_path']
    raw_data_path = config['raw_path']

    df_labels = pd.read_csv(label_path)
    path_of_data, labels = parse_labels(df_labels, raw_data_path)
    # path_of_data, labels = class_balance(np.array(path_of_data), np.array(labels))
    path_of_data = list(path_of_data)
    labels = list(labels)
    list_of_x = []
    list_of_y = []
    fileNames = []
    subjectList = []
    print('Begin Preprocessing the Files')
    
    cov = pd.read_csv(config['cov_path'])
    fitted = False
    add_genomics = config['add_genomics']
    filter_threshold = config['filter_threshold']
    # set logger to tqdm handler

    counter = 0
    for index, file in tqdm(enumerate(path_of_data), total = len(path_of_data), desc = "Processing"):
        
        subjectId = parse_file_name(file)
        if config['add_genomics']:
            if len(genomics.loc[genomics['subjectkey_intersection'] == subjectId]) == 0:
                print('No Genomics Info Found')
                continue
        try:
            output = SubsettingDay(file, day = config['subseted_days'], filter_threshold = filter_threshold)
        except AssertionError:
            print('Subsetting failure, value error')
            continue
        
        if output is None:
            continue
        
        output = add_cov(cov, output, subjectId, cov = names_)
        
        if add_genomics:
            output = add_cov(genomics, output, subjectId, 
                            cov = ['PRS_x', 'PRS_y'], cov_key = 'subjectkey_intersection')
        try:
            output = postprocessing(output, False)
        except:
            continue
        # except TypeError:
        #     print('TypeError in Postprocessing')
        #     continue
        if fitted == False:
            fitted = True
        if not checkNa(output):
            print('Error! Final output still containes NaNs')
            continue 
        print('Success')
        list_of_x.append(output)
        list_of_y.append(labels[index])
        fileNames.append(file)
        subjectList.append(subjectId)

    X = np.stack(list_of_x, axis = 0)
    Y = np.array(list_of_y)
    subjectList = np.array(subjectList)
    X = X.transpose([0, 2, 1])
    X_out = normalize(X, add_genomics = add_genomics)
    
    if add_genomics:
        np.save(os.path.join(config['outpath'],f'new_genomic_all_{filter_threshold}_{group}_non_clinical_2_3_X.npy'), X_out)
        np.save(os.path.join(config['outpath'],f'new_genomic_all_{filter_threshold}_{group}_non_clinical_2_3_subject.npy'), subjectList)
        np.save(os.path.join(config['outpath'],f'new_genomic_all_{filter_threshold}_{group}_non_clinical_2_3_Y.npy'), Y)
    else:
        np.save(os.path.join(config['outpath'],f'all_{filter_threshold}_{group}_non_clinical_2_3_X.npy'), X_out)
        np.save(os.path.join(config['outpath'],f'all_{filter_threshold}_{group}_non_clinical_2_3_subject.npy'), subjectList)
        np.save(os.path.join(config['outpath'],f'all_{filter_threshold}_{group}_non_clinical_2_3_Y.npy'), Y)
    
    print(f'Final data of size {X_out.shape[0]}')
    return X_out

