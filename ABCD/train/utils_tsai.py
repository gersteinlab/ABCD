import numpy as np
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
    'parent_div_cat_4.0', 'parent_div_cat_5.0', 'parent_div_cat_6.0', 'parent_div_cat_nan', 'race_eth_cat_1.0', 'race_eth_cat_2.0',
    'race_eth_cat_3.0', 'race_eth_cat_4.0', 'race_eth_cat_5.0', 'race_eth_cat_nan', 'sex_F', 'family_income',
    'interview_age', 'parent_grade'] , 'C2': ['famhis_bip', 'famhis_schiz', 'famhis_antisocial', 'famhis_nerves',
    'famhis_treatment', 'famhis_hospital', 'famhis_suicide','adopted'], 'cog': cog_task, 'C3': score, 'G' : genomic}

def class_balance(X, Y, subjectList = None, seed = 42):
    '''
    Takes X,Y as input, takes the X,Y balanced as the output
    '''
    if seed:
        np.random.seed(seed)
    else:
        pass
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
        if subjectList is not None:
            return X[final_idx], Y[final_idx], subjectList[final_idx]
        return X[final_idx], Y[final_idx]
    else:
        return X,Y

def find_index(query, key):
    key_value_dic = {k: i for i, k in enumerate(key)}
    return np.array([key_value_dic[q] for q in query])


def get_feature_names(mapped_from = ['ts', 'covs', 'cog_task', 'family_his', 'score', 'genomic', 'ts_ind']):
    key = []
    for k in mapped_from:
        key += eval(k)
    return key




def ablation_ts(X, mapped_to, mapped_from = ['ts', 'covs', 'cog_task', 'family_his', 'score', 'genomic', 'ts_ind']):
    '''
    Mapping function to subset time series.
    @params X: data
    @params mapped_to: a list of values containing following options
    W -> Wearables
    C1 -> Covariates set one mentioned in the paper
    C2 -> covariates set two mentioned in the paper
    C3 -> cognitive scores
    G -> Genomics
    @params mapped from: composition of the original time-series
    can be either:
    'ts', 'covs', 'cog_task', 'family_his', 'score', 'genomic', 'ts_ind'
    or
    'ts', 'covs', 'cog_task', 'family_his', 'score', 'ts_ind'
    '''


    query = []
    key = get_feature_names(mapped_from=mapped_from)
    assert X.shape[1] == len(key), "Mapped from items do not match X size"
    for m in mapped_to:
        query += grouping[m]
    index = find_index(query, key)
    return index, query