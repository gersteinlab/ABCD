cat = ['parent_div_cat', 'race_eth_cat','sex']
encoded_cats = [ 'famhis_bip',
 'famhis_schiz',
 'famhis_antisocial',
 'famhis_nerves',
 'famhis_treatment',
 'famhis_hospital',
 'famhis_suicide',
 'parent_div_cat',
 'parent_grade',
 'adopted']

quan = ['family_income', 'interview_age']
target_1 = ['cbcl_externalizing',
 'cbcl_internalizing']

target_2 =['stopsigreactiontime',
 'stopsignalgort',
 'standevgort',
 'picvocab',
 'flanker',
 'processspeed',
 'picmemory',
 'readingscore']

names_ = ['parent_div_cat_1.0', 'parent_div_cat_2.0', 'parent_div_cat_3.0',
       'parent_div_cat_4.0', 'parent_div_cat_5.0', 'parent_div_cat_6.0',
       'parent_div_cat_nan', 'race_eth_cat_1.0', 'race_eth_cat_2.0',
       'race_eth_cat_3.0', 'race_eth_cat_4.0', 'race_eth_cat_5.0',
       'race_eth_cat_nan', 'sex_F', 'sex_M', 'family_income',
       'interview_age', 'stopsigreactiontime',
 'stopsignalgort',
 'standevgort',
 'picvocab',
 'flanker',
 'processspeed',
 'picmemory',
 'readingscore','famhis_bip', 'famhis_schiz', 'famhis_antisocial', 'famhis_nerves',
       'famhis_treatment', 'famhis_hospital', 'famhis_suicide',
       'parent_grade', 'adopted'] + target_1

names_reg = ['parent_div_cat_1.0', 'parent_div_cat_2.0', 'parent_div_cat_3.0',
       'parent_div_cat_4.0', 'parent_div_cat_5.0', 'parent_div_cat_6.0',
       'parent_div_cat_nan', 'race_eth_cat_1.0', 'race_eth_cat_2.0',
       'race_eth_cat_3.0', 'race_eth_cat_4.0', 'race_eth_cat_5.0',
       'race_eth_cat_nan', 'sex_F', 'sex_M', 'family_income',
       'interview_age', 
       'famhis_bip', 'famhis_schiz', 'famhis_antisocial', 'famhis_nerves',
       'famhis_treatment', 'famhis_hospital', 'famhis_suicide',
       'parent_grade', 'adopted']