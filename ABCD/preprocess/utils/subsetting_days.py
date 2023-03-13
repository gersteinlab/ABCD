from sklearn.preprocessing import LabelEncoder
import numpy as np
import os
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
try:
    from sktime.transformations.series.impute import Imputer
except:
    print('No Sktime Package found, please install it to process the data.')
from .features import *
warnings.filterwarnings('default')

le = LabelEncoder()
les = [LabelEncoder()]

def class_balance(X, Y):
    '''
    Takes X,Y as input, takes the X,Y balanced as the output
    '''
    num_positive_class = np.sum(Y == 1)
    num_negative_class = np.sum(Y == 0)
    min_class = min(num_positive_class, num_negative_class)
    min_flag = 1 if min_class == num_positive_class else 0
    max_flag = 1 - min_flag
    frac = min_class /(num_positive_class+num_negative_class)
    if frac <= .4:
        idx= np.where(Y == max_flag)[0]
        idx_idx = np.arange(len(idx))
        np.random.shuffle(idx_idx)
        subsetted_shuffled_idx_idx = idx_idx[:min_class]
        idx_subsetted = idx[subsetted_shuffled_idx_idx]
        final_idx = np.concatenate([idx_subsetted, np.where(Y == min_flag)[0]])
        return X[final_idx], Y[final_idx]
    else:
        return X,Y


def get_path_from_subjectID(ids, data_path):
    return [os.path.join(data_path, id + '_combined.csv') for id in ids]
    
    
def parse_labels(df, data_path):
    df['label'] = le.fit_transform(df['label'])
    return get_path_from_subjectID(df['subjectID'].to_list(), data_path), df['label'].to_list()    
    
    
def SubsettingDay(fileName, day = [1,2], filter_threshold = 11000): 
    '''
    This function takes a csv file as an input, and return => a dataframe with desirable measurements.
    The desirable measurements are a subset of the original raw meansurements which contains:
    @param day: a list or int, specify which day/days (continuous days in a week) of the measurements to contain
    @param filter_threshold: acts as a filter => Each of the day has to have # of measurements greater than some number.
    '''
    day = set(day) if isinstance(day, list) else {day}
    try:
        a = pd.read_csv(fileName)
    except FileNotFoundError:
        print('Subsetting Failure. Original Data Not Found')
        return None
    except Warning:
        print(fileName)
        return None
    try:
        a['Wear_Time'] = pd.to_datetime(a['Wear_Time'])
    except KeyError:
        a['Wear_Time'] = a['Time']
        a['Wear_Time'] = pd.to_datetime(a['Wear_Time'])
    
    a = a.sort_values(by = 'Wear_Time').reset_index().drop(['index'], axis = 1)
    a['DayOfWeek'] = a['Wear_Time'].dt.dayofweek
    a['NumOfWeek'] = a['Wear_Time'].dt.isocalendar().week
    try:
        df = a.groupby(['DayOfWeek', 'NumOfWeek']).filter(lambda x: x.notna().sum().sum() > filter_threshold). \
        groupby(['DayOfWeek', 'NumOfWeek']).apply(lambda x: x.notna().sum().sum()).reset_index(). \
        rename(columns={0:'count'})
    except ValueError:
        print('Subsetting Failure. Failed at subsetting function. The failed filename is ' + os.path.basename(fileName))
        return None
    df = df.loc[df['DayOfWeek'].isin(day)]
    df = df[['NumOfWeek','count']].groupby('NumOfWeek').filter(lambda sb: len(sb) == len(day))
    if df.empty:    
        print('Subsetting Failure. No satisfiable day frames founded. The failed filename is ' + os.path.basename(fileName))
        return None
    targetingNumOfWeek = df.sort_values(by = ['count']).iloc[0,:]['NumOfWeek']
    ret =  a.loc[(a['DayOfWeek'].isin(day)) & (a['NumOfWeek'] == targetingNumOfWeek)]
    assert len(ret) == len(day) * 1440
    return ret.iloc[:,:-2]


def postprocessing(df, fitted = False, categorical = ['Level']):
    from sktime.datasets import load_airline
    from sktime.forecasting.arima import ARIMA

    df = df.drop(['SleepStage', 'ShortWakes'], axis = 1)
    imputation_dic = {'Intensity': Imputer(method="constant", missing_values = 0), 
                      'SleepLevel': Imputer(method="constant", missing_values = 0),
                     'Steps': Imputer(method="constant", missing_values = 0),
                      'Calories': Imputer(method="constant", missing_values = df['Calories'].min()),
                     }                  
    int_mask = df.isna().iloc[:, 1:8]
    if not fitted: # sklearn fitting
        for index, cat in enumerate(categorical):
            df[cat] = les[index].fit_transform(df[cat]) 
    else:
        for index, (cat) in enumerate(categorical):
            df[cat] = les[index].transform(df[cat])
    # imp_mean = Imputer(method="forecaster", forcaster = ARIMA(order=(1, 1, 0), seasonal_order=(0, 1, 0, 12), suppress_warnings=True))
    imp_mean = Imputer(method ='drift')
    # df_fitted = imp_mean.fit_transform(df)
    # df_fitted = df
    
    if not df_fitted.isna().any().any():
        imp_mean = Imputer(method="linear")
        df_fitted = imp_mean.fit_transform(df_fitted)
    if not df_fitted.isna().any().any():
        imp_mean = Imputer(method="bfill")
        df_fitted = imp_mean.fit_transform(df_fitted)
        imp_mean = Imputer(method="ffill")
        df_fitted = imp_mean.fit_transform(df_fitted)
    out = np.concatenate([df_fitted.iloc[:, 1:], int_mask.to_numpy()], axis = 1) 
    return out

def checkNa(df):
    return np.isnan(df).sum() == 0


def add_cov(cov_df, data_df, subjectkey, cov = names_ + target_1, cov_key = 'subjectkey'):
    val = list(cov_df.loc[cov_df[cov_key] == subjectkey][cov].values)[0]
        
    for col_name, value in zip(cov, val):
        data_df[col_name] = value
    return data_df




def get_target(cov_df, subjectkey, cov = target_1 + target_2):
    return np.array(list(cov_df.loc[cov_df['subjectkey'] == subjectkey][cov].values)[0])

get_subject_id_from_file_name = lambda file: re.search('(NDAR_\w*)_combined\.csv', file).group(1)
get_subject_id_from_file_name = lambda file: re.search('(NDAR\w*)_combined\.csv', file).group(1)


def parse_file_name(file):
    file = re.search('(NDAR\w*)_combined\.csv', file).group(1)
    if file[4] != '_':
        file = file[:4] + '_' + file[4:]
    return file


def find_index(query, key):
    key_value_dic = {k: i for i, k in enumerate(key)}
    return np.array([key_value_dic[q] for q in query])

def ablation_ts(mode):
    ts = ['Intensity', 'SleepLevel','Steps','Calories','HeartRate','SleepValue','METs']
    covs =['parent_div_cat_1.0', 'parent_div_cat_2.0', 'parent_div_cat_3.0',
       'parent_div_cat_4.0', 'parent_div_cat_5.0', 'parent_div_cat_6.0',
       'parent_div_cat_nan', 'race_eth_cat_1.0', 'race_eth_cat_2.0',
       'race_eth_cat_3.0', 'race_eth_cat_4.0', 'race_eth_cat_5.0',
       'race_eth_cat_nan', 'sex_F', 'sex_M', 'family_income',
       'interview_age', ]
    cog_task = ['stopsigreactiontime','stopsignalgort', 'standevgort','picvocab','flanker','processspeed',
 'picmemory','readingscore',]   
    family_his = ['famhis_bip', 'famhis_schiz', 'famhis_antisocial', 'famhis_nerves',
       'famhis_treatment', 'famhis_hospital', 'famhis_suicide',
       'parent_grade', 'adopted',]
    score = ['cbcl_externalizing',
 'cbcl_internalizing']
    genomic = ['PRS_x', 'PRS_y']
    ts_ind = [t + '_ind' for t in ts]
    grouping = {'W': ts + ts_ind, 'C1': ['parent_div_cat_1.0', 'parent_div_cat_2.0', 'parent_div_cat_3.0',
       'parent_div_cat_4.0', 'parent_div_cat_5.0', 'parent_div_cat_6.0', 'race_eth_cat_1.0', 'race_eth_cat_2.0',
       'race_eth_cat_3.0', 'race_eth_cat_4.0', 'race_eth_cat_5.0', 'sex_F', 'family_income',
       'interview_age', 'parent_grade'] , 'C2': ['famhis_bip', 'famhis_schiz', 'famhis_antisocial', 'famhis_nerves',
       'famhis_treatment', 'famhis_hospital', 'famhis_suicide','adopted'], 'C3': score, 'G' : genomic}
    query = []
    key = ts + covs + cog_task + family_his + score + ts_ind
    assert X.shape[1] == len(key)
    for m in mode:
        query += grouping[m]
    index = find_index(query, key)
    return index, query