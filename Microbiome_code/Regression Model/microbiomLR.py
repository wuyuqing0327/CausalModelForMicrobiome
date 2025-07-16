import pandas as pd
import numpy as np
path=r'C:\Users\yuqingw1\Workfolder'

### read the datasets
df = pd.read_csv(path + '/Data/compass_deidentified.csv')
print(df.shape)
print(sum(pd.isnull(df['specimen_id'])))

### remove the nan specimen_id
df = df[~pd.isnull(df['specimen_id'])]
df = df.reset_index()
df['specimen_id'] = df['specimen_id'].astype(int).astype(str)
print(df.shape)
print(sum(pd.isnull(df['specimen_id'])))


### check systolic >> diastolic

df[df['systolic']<=df['diastolic']][['specimen_id','systolic', 'diastolic']]
 # .to_csv('incorrectblood.csv', index=False))

# ### remove the incorrect blood datasets
# print("the incorrect blood datasets shape:", df[(abs(df['systolic'] - df['diastolic'])>15)].shape)
# df = df[(abs(df['systolic'] - df['diastolic'])>15)]
# df = df.reset_index()
### exchange systolic and diastolic
mask = df['systolic']<=df['diastolic']
df.loc[mask, ['systolic', 'diastolic']] = df.loc[mask, ['diastolic', 'systolic']].values
print(df[df['systolic']<=df['diastolic']][['specimen_id','systolic', 'diastolic']])

### processing age gender bmi missing value
df[pd.isnull(df['age'])][['specimen_id','age']].shape
print(df['age_binned'].value_counts())


# Bin numeric data, treating NaN as a separate category
cutPoint = list(np.unique(np.percentile(df[~pd.isnull(df['age'])]['age'], range(10, 100, 10))))
df['age_binned2'] = pd.cut(df['age'], bins=[-np.inf] + cutPoint+ [np.inf], include_lowest=True)
df['age_binned2'] = np.where(pd.isnull(df['age_binned2']), 'missing', df['age_binned2'])

cutPoint = list(np.unique(np.percentile(df[~pd.isnull(df['bmi'])]['bmi'], range(10, 100, 10))))
df['bmi_binned2'] = pd.cut(df['bmi'], bins=[-np.inf] + cutPoint+ [np.inf], include_lowest=True)
df['bmi_binned2'] = np.where(pd.isnull(df['bmi_binned2']), 'missing', df['bmi_binned2'])

df['gender'] = df['gender'].astype(str)
df['gender'] = np.where(pd.isnull(df['gender'])|(df['gender']=='nan'), 'missing', df['gender'])
df['gender'] = np.where(df['gender']=='I prefer not to answer', 'missing', df['gender'])




###### merge data
df['specimen_id'] = df['specimen_id'].astype('int').astype('str')
finaldf_trans['sampleID'] = finaldf_trans['sampleID'].astype('int').astype('str')


datamodel = finaldf_trans.merge(df[['specimen_id', 'systolic', 'diastolic', 'gender', 'bmi', 'bmi_binned2', 'age', 'age_binned2']], left_on=['sampleID'], right_on=['specimen_id'])

print(datamodel.shape)


### check the repeat ID
from collections import Counter
item_counts = Counter(list(datamodel['sampleID']))
# Finding items with more than one occurrence
repeat_items = [item for item, count in item_counts.items() if count > 1]
print("data408_416 Repeated items:", repeat_items)

item_counts = Counter(list(datamodel['sampleID']))
repeat_items = [item for item, count in item_counts.items() if count > 1]
print("data426_428 Repeated items:", repeat_items)

item_counts = Counter(list(datamodel['sampleID']))
repeat_items = [item for item, count in item_counts.items() if count > 1]
print("datamicrbio Repeated items:", repeat_items)

df_repeat = datamodel[(datamodel['sampleID']=='104726') | (datamodel['sampleID']=='106503') | (datamodel['sampleID']=='106615') | (datamodel['sampleID']=='108359')]
df_repeat[['sampleID', 'systolic', 'diastolic', 'gender', 'bmi', 'age']]





########################################################################################################
### processing
dataproc408_416 = pd.get_dummies(datamodel, columns=['gender', 'bmi_binned2', 'age_binned2'])


### regression model

del_cols = ['specimen_id', 'sampleID', 'systolic', 'diastolic', 'gender', 'bmi', 'bmi_binned2', 'age', 'age_binned2']
target1 = 'systolic'
target2 = 'diastolic'
# base_cols= list(set(dataproc408_416.columns) - set(data408_416.columns))
base_cols = ['age_binned2_(-inf, 38.19]',
 'age_binned2_(38.19, 43.61]',
 'age_binned2_(43.61, 48.2]',
 'age_binned2_(48.2, 51.328]',
 'age_binned2_(51.328, 54.1]',
 'age_binned2_(54.1, 56.69]',
 'age_binned2_(56.69, 59.36]',
 'age_binned2_(59.36, 62.54]',
 'age_binned2_(62.54, 67.41]',
 'bmi_binned2_(-inf, 20.884]',
 'bmi_binned2_(20.884, 22.864]',
 'bmi_binned2_(22.864, 24.531]',
 'bmi_binned2_(24.531, 26.232]',
 'bmi_binned2_(26.232, 27.891]',
 'bmi_binned2_(27.891, 29.662]',
 'bmi_binned2_(29.662, 31.75]',
 'bmi_binned2_(31.75, 34.497]',
 'bmi_binned2_(34.497, 38.965]',
 'bmi_binned2_(38.965, inf]',
 'gender_Male',
 'gender_Female']

base_cols2 = ['age','bmi',
    'gender_Male',
    'gender_Female']

features = list(set(dataproc408_416.columns) - set([target1, target2]) - set(del_cols) - set(base_cols))
print(len(features))

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

def Regressionforblood(data, target, base_cols, features, method):

    # X_train, X_test, y_train, y_test = train_test_split(data[features+base_cols], data[target], test_size=1,
    #                                                     random_state=2024)
    X_train, X_test, y_train, y_test = data[features+base_cols], data[features+base_cols], data[target], data[target]
    select_var=[]
    select_coef=[]
    pvalue=[]
    for i in features:
        cols = [i] + base_cols
        X = sm.add_constant(X_train[cols])
        model = sm.OLS(y_train, X).fit()
        f_pvalue = model.f_pvalue
        if f_pvalue > 0 and f_pvalue < 0.05:
            # print("feature:", i)
            select_var.append(i)
            select_coef.append(model.params[i])
            pvalue.append(f_pvalue)
            # print("pvalue:", f_pvalue)

    # Apply FDR correction using the Benjamini-Hochberg method
    rejected, p_values_corrected, _, _ = multipletests(pvalue, alpha=0.05, method=method, is_sorted=False)
    print("resulttttt:", sum(rejected))
    condition = p_values_corrected < 0.05
    select_var = list(np.array(select_var)[condition])
    select_coef = list(np.array(select_coef)[condition])
    pvalue = list(np.array(p_values_corrected)[condition])
    print("lenpvalue:", len(pvalue))

    return select_var,select_coef,pvalue



varSBD, coefSBD, pvalueSBD = Regressionforblood(dataproc408_416, target=target1, base_cols=base_cols,  features=features, method='fdr_bh') ###'fdr_bh'

varDBD, coefDBD, pvalueDBD = Regressionforblood(dataproc408_416, target=target2, base_cols=base_cols,  features=features, method='fdr_bh')






print(len(varSBD))
print(len(varDBD))
