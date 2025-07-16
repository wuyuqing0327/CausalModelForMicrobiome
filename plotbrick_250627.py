import pandas as pd
import numpy as np

path = '/Users/chelseawu/UChicago/Part-time/BSD/result/20250627/'
data_path = '/Users/chelseawu/UChicago/Part-time/BSD/Datasets/'

#
# pm25 = pd.read_csv(path + 'ancombc_result_PM25_prv01.csv')
# ghrelin = pd.read_csv(path + 'ancombc_result_Ghrelin_prv01.csv')
# resistin = pd.read_csv(path + 'ancombc_result_Resistin_prv01.csv')
# insluin = pd.read_csv(path + 'ancombc_result_Insulin_prv01.csv')
#
#
# a1 = pm25[pm25['q_M12_v2']<0.05][['taxon', 'lfc_M12_v2', 'q_M12_v2']]
# a2 = ghrelin[ghrelin['q_Ghrelin_v3']<0.05][['taxon', 'lfc_Ghrelin_v3', 'q_Ghrelin_v3']]
# a3 = resistin[resistin['q_Resistin_v3']<0.05][['taxon', 'lfc_Resistin_v3', 'q_Resistin_v3']]
# a4 = insluin[insluin['q_Insulin_v3']<0.05][['taxon', 'lfc_Insulin_v3', 'q_Insulin_v3']]
#
# sig_micro = list( set(set(a1['taxon']) & set(a2['taxon'])) | set(set(a1['taxon']) & set(a3['taxon']))
#       | set(set(a1['taxon']) & set(a4['taxon'])) )
# len(sig_micro)
#
#
# result = pm25[np.in1d(pm25['taxon'], sig_micro)][['taxon', 'lfc_M12_v2', 'q_M12_v2']].\
#     merge(ghrelin[np.in1d(ghrelin['taxon'], sig_micro)][['taxon', 'lfc_Ghrelin_v3', 'q_Ghrelin_v3']], on=['taxon'], how='left').\
#     merge(resistin[np.in1d(resistin['taxon'], sig_micro)][['taxon', 'lfc_Resistin_v3', 'q_Resistin_v3']], on=['taxon'], how='left').\
#     merge(insluin[np.in1d(insluin['taxon'], sig_micro)][['taxon', 'lfc_Insulin_v3', 'q_Insulin_v3']], how='left')
#
#
# result.to_csv(path + 'sig_micro.csv', index=False)
#

data = pd.read_csv(path + 'sig_micro.csv')
taxonIDname = pd.read_excel(data_path + 'taxonIDmappingdict.xlsx')

taxonIDname['Genus'] = taxonIDname['attribute'].apply(lambda x: x.split(',')[-1])
data.merge(taxonIDname, left_on='taxon', right_on='Genus')







