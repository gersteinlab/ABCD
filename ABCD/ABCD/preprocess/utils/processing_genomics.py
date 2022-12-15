import pandas as pd
from ABCD.preprocess.utils import *
def processing_genomics(genomics_phenotypes_path_1 = '/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/bb926/PRSs/adhd2019/mj.pipeline/ldpred2/RData/ldpred2_prs_auto.out.tsv',
    genomics_phenotypes_path_2 = '/gpfs/gibbs/pi/gerstein/jjl86/data/ABCD/bb926/PRSs/adhd2019/mj.pipeline/ldpred2/RData/ldpred2_prs_infinitesimal.out.tsv'):
    genomics_phenotypes_1 = pd.read_csv(genomics_phenotypes_path_1, sep = '\t')
    genomics_phenotypes_2 = pd.read_csv(genomics_phenotypes_path_2, sep = '\t')
    genomics_phenotypes_1 = genomics_phenotypes_1.reset_index()
    genomics_phenotypes_1['subjectkey_intersection'] = genomics_phenotypes_1['index'].str.extract(r'\w*_(NDAR.*)')
    genomics_phenotypes_2 = genomics_phenotypes_2.reset_index()
    genomics_phenotypes_2['subjectkey_intersection'] = genomics_phenotypes_2['index'].str.extract(r'\w*_(NDAR.*)')
    genomics = genomics_phenotypes_1.merge(genomics_phenotypes_2, on = 'subjectkey_intersection')[['subjectkey_intersection','PRS_x','PRS_y']]
    return genomics