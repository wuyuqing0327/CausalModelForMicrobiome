import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import t

# Open the text file for metabolitic
df1 = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\Data\Breast_Cancer\GCMS_feature_04252022.csv')
print("df1 shape:", df1.shape)
df1['metaid'] = 'v' + df1['metaid'].astype(str)




df2 = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\Data\Breast_Cancer\vip_matrix.csv')
print("df2 shape:", df2.shape)
df2.rename(columns={'MetabolomicsID':'metaid'}, inplace=True)


for name in ['White_alone_not_Hispanic_or_Latino_rate_binary',
       'high_school_and_lessthen_rate_binary',
       'poverty_below_100_rate_binary']:
    print("name:", name)
    data = df1[['metaid', 'Mass', 'RT']].merge(df2[['metaid', name]], on=['metaid'])
    print("datashape:", data.shape)

    data['p-value'] = np.where(data[name]>=2, 0.03, 2)
    data['statistic'] = data['p-value'].apply(lambda p: t.ppf(1 - p / 2, 416))
    data['statistic'] = np.where(np.isinf(data['statistic']), 1000, data['statistic'])
    data.rename(columns={'Mass':'m/z', 'RT': 'retention time'}, inplace=True)
    # print("datacols:", data.columns)
    data = data.drop_duplicates(['m/z', 'retention time', 'p-value', 'statistic'])

    data = data[['m/z', 'retention time', 'p-value', 'statistic']]
    data['m/z'] = pd.to_numeric(data['m/z'], errors='coerce')
    data['retention time'] = pd.to_numeric(data['retention time'], errors='coerce')
    data['p-value'] = pd.to_numeric(data['p-value'], errors='coerce')
    data['statistic'] = pd.to_numeric(data['statistic'], errors='coerce')

    data.columns = ['mz', 'rtime', 'p-value', 't-score']

    print("datashape final:", data.shape)
    data.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Breast_Cancer\mummichog_%s.txt' %(name), sep='\t', index=False)




