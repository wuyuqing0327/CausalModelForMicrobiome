import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import t

# Open the text file for metabolitic
df1 = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\MSC1445_list_filterID.csv')
print("df1 shape:", df1.shape)
df1 = df1.reset_index()
df1.rename(columns={'index':'MetaID'}, inplace=True)
# df1[['MetaID', 'Mass', 'Retention Time']].shape
df1['MetaID'] = 'Meta_' + df1['MetaID'].astype(str)



df2 = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\vip_scores_filterID.csv')
print("df2 shape:", df2.shape)
df2.rename(columns={'Metabolite':'MetaID'}, inplace=True)

for name in [f"{i:03}" for i in range(1, 486)]:
    print("name:", name)
    data = df1[['MetaID', 'Mass', 'Retention Time']].merge(df2[['MetaID', 'Comp1_UPDASV%s' %name]], on=['MetaID'])
    print("datashape:", data.shape)

    data['p-value'] = np.where(data['Comp1_UPDASV%s' %name]>2, 0.03, 2)
    data['statistic'] = data['p-value'].apply(lambda p: t.ppf(1 - p / 2, 416))
    data['statistic'] = np.where(np.isinf(data['statistic']), 1000, data['statistic'])
    data.rename(columns={'Mass':'m/z', 'Retention Time': 'retention time'}, inplace=True)
    # print("datacols:", data.columns)
    data = data.drop_duplicates(['m/z', 'retention time', 'p-value', 'statistic'])

    data = data[['m/z', 'retention time', 'p-value', 'statistic']]
    data['m/z'] = pd.to_numeric(data['m/z'], errors='coerce')
    data['retention time'] = pd.to_numeric(data['retention time'], errors='coerce')
    data['p-value'] = pd.to_numeric(data['p-value'], errors='coerce')
    data['statistic'] = pd.to_numeric(data['statistic'], errors='coerce')

    data.columns = ['mz', 'rtime', 'p-value', 't-score']

    print("datashape final:", data.shape)
    data.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\mummichog\microbiome_%s.txt' %(name), sep='\t', index=False)

#
# import mummichog.functional_analysis as fa
#
# from mummichog1.functional_analysis import PathwayAnalysis
#
# analysis = fa.PathwayAnalysis(data['m/z'], data['p-value'], 'positive')
# # Run the analysis
# results = analysis.run()
#
# # Output or further process the results
# print(results)