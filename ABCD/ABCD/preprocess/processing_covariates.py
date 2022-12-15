import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import StandardScaler
def processing_covariates():
    df = pd.read_csv('/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/metadata/ABCD_metadata_06082022.tsv', sep = '\t')
    df_2 = pd.read_csv('/home/yl2428/ABCD/ABCD_metadata_03152022.csv')[['subjectkey', 'interview_age','sex']]
    df = df.merge(df_2, left_on = 'src_subject_id', right_on = 'subjectkey')
    column_trans = ColumnTransformer(
        [('categories', OneHotEncoder(dtype='int'), cat),
        ('quant', StandardScaler(), quan)],
        remainder='passthrough', verbose_feature_names_out=False)
    column_trans.fit(df)
    names = column_trans.get_feature_names_out()
    arr = column_trans.transform(df)
    df = pd.DataFrame(arr, columns=[column_trans.get_feature_names_out()])
    df = df.drop(['src_subject_id'], axis = 1)
    df.to_csv('covariates_processed_reg.csv')
    return df