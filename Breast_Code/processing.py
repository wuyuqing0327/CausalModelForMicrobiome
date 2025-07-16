import pandas as pd
import numpy as np

df1 = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\Data\Breast_Cancer\data_combined_12152022.csv')
print(df1.columns)
print(df1.shape)
set(df1['zipcode'])
df1['zipcode'] = np.where(df1['zipcode']==60680, 60606, df1['zipcode'])

df2 = pd.read_excel(r'C:\Users\yuqingw1\Workfolder\Data\Breast_Cancer\Geographicdata.xlsx')
print(df2.shape)
print(df2.columns)
df2['White_alone_not_Hispanic_or_Latino_rate'] = df2['White_alone_not_Hispanic_or_Latino']/df2['total_population']
df2['high_school_and_lessthen_rate'] = df2['high_school_and_lessthen']/df2['total_population']
# df2['salary_less_75k_rate'] = df2['salary_less_75k']/df2['total_population']
df2['poverty_below_100_rate'] = df2['poverty_below_100']/df2['total_population']


data = df1.merge(df2[['zipcode',  'White_alone_not_Hispanic_or_Latino_rate',
       'high_school_and_lessthen_rate',
       'poverty_below_100_rate']], on=['zipcode'], how='left')

selectcols = ['StudyID', 'zipcode','Age_serum', 'marital_status', 'tobacco_history', 'alcohol_history', 'T2D', 'adi_natrank',
              'White_alone_not_Hispanic_or_Latino_rate',
            'high_school_and_lessthen_rate', 'poverty_below_100_rate']
data[selectcols + list(df1.columns[67:])].to_csv(r'C:\Users\yuqingw1\Workfolder\Data\Breast_Cancer\breastcancer_pls_data.csv', index=False)











