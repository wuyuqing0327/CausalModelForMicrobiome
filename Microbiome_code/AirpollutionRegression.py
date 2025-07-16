import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.model_selection import cross_val_score
from scipy import stats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import re

path = 'Datasets/'
save_path = 'result/202506new/microbiome_filtermissingover50/'
trans = pd.read_csv(path + 'microbiome_allvariable_806_transformation.csv')
# shannon = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\processing_4_microbiome_new\microbiome_combine_806row_485col_shannon.csv')


print(trans['smoking_status'].value_counts())

print(trans['diabetes2'].value_counts())



trans = trans[~pd.isnull(trans['sampleID'])]
trans['sampleID'] = trans['sampleID'].astype(int).astype(str)
print(trans.shape)

# shannon = shannon[~pd.isnull(shannon['sampleID'])]
# shannon['sampleID'] = shannon['sampleID'].astype(int).astype(str)
# print(shannon.shape)


data = trans.copy()
# data = trans.merge(shannon, on=['sampleID'])

data = data[~ (pd.isnull(data['M12']) | pd.isnull(data['M18']) | pd.isnull(data['M6']))]
print(data.shape)



cutPoint = [18.5, 25, 30]
data['bmi_binned2'] = pd.cut(data[~pd.isnull(data['bmi'])]['bmi'],bins=[-np.inf] + cutPoint + [np.inf]).astype(str)
data['bmi_binned2'] = np.where(pd.isnull(data['bmi_binned2']), 'missing', data['bmi_binned2'])
print(data['bmi_binned2'].value_counts())
data['bmi_binned2_30.0_inf'] = np.where(data['bmi_binned2']=='(30.0, inf]', 1, 0)
data['bmi_binned2_25.0_30.0'] = np.where(data['bmi_binned2']=='(25.0, 30.0]', 1, 0)
data['bmi_binned2_inf_18.5'] = np.where(data['bmi_binned2']=='(-inf, 18.5]', 1, 0)
data['bmi_binned2_missing'] = np.where(data['bmi_binned2']=='missing', 1, 0)

print(data.shape)
data['gender_Male']=np.where(data['gender']=='Male', 1, 0)
data['smoking_status_Current']=np.where(data['smoking_status']=='Current', 1, 0)
data['smoking_status_Former']=np.where(data['smoking_status']=='Former', 1, 0)
data['diabetes2_Yes'] = np.where(data['diabetes2']=='Yes', 1, 0)
print(data.shape)

data['gender'] = np.where(pd.isnull(data['gender']), 'missing', data['gender'])
data['race'] = np.where(pd.isnull(data['race']), 'missing', data['race'])
data['education'] = np.where(pd.isnull(data['education']), 'missing', data['education'])
data['household_income'] = np.where(pd.isnull(data['household_income']), 'missing', data['household_income'])

print(data[['age', 'bmi']].describe())
for i in ['bmi_binned2', 'gender', 'race', 'household_income', 'education']:
    print(data[i].value_counts())

cols = ['sampleID', 'age', 'gender_Male',
        'bmi_binned2_30.0_inf','bmi_binned2_25.0_30.0','bmi_binned2_inf_18.5',
        'smoking_status_Current','smoking_status_Former',
        'diabetes2_Yes',
         'M6', 'M12', 'M18', 'M24']

# Microcols = [f"UPDASV{i:03}" for i in range(1, 486)]
Microcols = ['UPDASV030', 'UPDASV051', 'UPDASV053', 'UPDASV065', 'UPDASV073', 'UPDASV139', 'UPDASV169', 'UPDASV184', 'UPDASV186', 'UPDASV187', 'UPDASV190', 'UPDASV197', 'UPDASV224', 'UPDASV236', 'UPDASV250', 'UPDASV264', 'UPDASV270', 'UPDASV285', 'UPDASV291', 'UPDASV302', 'UPDASV310', 'UPDASV345', 'UPDASV358', 'UPDASV359', 'UPDASV372', 'UPDASV374', 'UPDASV385', 'UPDASV387', 'UPDASV394', 'UPDASV402', 'UPDASV413', 'UPDASV424', 'UPDASV429', 'UPDASV430', 'UPDASV434', 'UPDASV442', 'UPDASV467', 'UPDASV480']
# Microcols_filtermissing75 = ['UPDASV024', 'UPDASV030', 'UPDASV039', 'UPDASV051', 'UPDASV053', 'UPDASV060', 'UPDASV061', 'UPDASV065', 'UPDASV073', 'UPDASV139', 'UPDASV150', 'UPDASV153', 'UPDASV155', 'UPDASV169', 'UPDASV176', 'UPDASV184', 'UPDASV186', 'UPDASV187', 'UPDASV190', 'UPDASV197', 'UPDASV198', 'UPDASV207', 'UPDASV224', 'UPDASV234', 'UPDASV236', 'UPDASV243', 'UPDASV250', 'UPDASV264', 'UPDASV268', 'UPDASV270', 'UPDASV273', 'UPDASV285', 'UPDASV291', 'UPDASV296', 'UPDASV302', 'UPDASV310', 'UPDASV345', 'UPDASV350', 'UPDASV358', 'UPDASV359', 'UPDASV361', 'UPDASV366', 'UPDASV372', 'UPDASV374', 'UPDASV376', 'UPDASV385', 'UPDASV387', 'UPDASV388', 'UPDASV391', 'UPDASV393', 'UPDASV394', 'UPDASV398', 'UPDASV401', 'UPDASV402', 'UPDASV413', 'UPDASV424', 'UPDASV429', 'UPDASV430', 'UPDASV434', 'UPDASV442', 'UPDASV444', 'UPDASV447', 'UPDASV450', 'UPDASV467', 'UPDASV480', 'UPDASV485']

Microcols = list(set(Microcols) - set(['UPDASV186'])) # remove 'UPDASV186' because (具体不太清楚，可能是这个microbiome与某个microbiome重复了）
data = data[cols+Microcols]
# data = data[cols+['Shannon_Diversity']]
print("final data shape:", data.shape)


ylist = Microcols
# ylist2=['Shannon_Diversity']
air_col = ['M6', 'M12', 'M18', 'M24']
restcols=list(set(data.columns) - set(ylist) - set(['sampleID']) - set(air_col))
# restcols=list(set(data.columns) - set(ylist2) - set(['sampleID']) - set(air_col))
print(len(restcols))
featuresall = []
for i in ['M12']:
    features_tmp = restcols+ [i]
    featuresall.append(features_tmp)

traindata = data[air_col+ylist+restcols]

print("traindata shape:", traindata.shape)



#### processing + Regression
def airpollutionreg(data, features, ylist, microcols=Microcols):
    # X_train, X_test, y_train_all, y_test_all = train_test_split(data[features], data[ylist], test_size=0.1, random_state=42)
    X_train = data[features]
    X_test = data[features]
    y_train_all = data[ylist]
    y_test_all = data[ylist]

    X_train['age'] = np.where((X_train['age'] == 0) | (pd.isnull(X_train['age'])), X_train['age'].mean(),
                              X_train['age'])
    X_test['age'] = np.where((X_test['age'] == 0) | (pd.isnull(X_test['age'])), X_train['age'].mean(), X_test['age'])

    restcols = list(set(X_train.columns) - set(microcols))
    print(" microcols len:", len(microcols))
    print(" restcols len:", len(restcols))

    dropcols = []
    for i in X_train.columns:
        if len(np.unique(X_train[i])) == 1:
            dropcols.append(i)
    print("dropcols:", dropcols)
    X_train.drop(dropcols, axis=1, inplace=True)
    X_test.drop(dropcols, axis=1, inplace=True)

    model_coef = {}
    model_pvalue = {}
    model_error = {}
    for yi in ylist:
        # print("current target:", yi)
        y_train = y_train_all[yi]
        y_test = y_test_all[yi]


        # Initialize models
        X_train = sm.add_constant(X_train)
        model = sm.OLS(y_train, X_train).fit()

        # Display results
        model_coef[yi] = list(model.params)
        model_pvalue[yi] = list(model.pvalues)
        model_error[yi] = list(model.bse)

    pattern = r'^M\d+$'
    replacement = 'PM25'
    orgname = list(X_train.columns)
    updname = [replacement if re.match(pattern, item) else item for item in orgname]

    model_coef_df = pd.DataFrame(model_coef).transpose()
    model_coef_df.columns = updname
    model_pvalue_df = pd.DataFrame(model_pvalue).transpose()
    model_pvalue_df.columns = updname
    model_error_df = pd.DataFrame(model_error).transpose()
    model_error_df.columns = updname

    # adjust p_value utilizing
    model_fdrpvalue = []
    pval = model_pvalue_df
    pval = pval.fillna(0)
    # print("pval:", pval)
    for i in list(model_pvalue_df.columns):
        eachpval = list(pval[i])
        # print("eachpval:", eachpval)
        _, each_pvalues_adjusted, _, _ = multipletests(eachpval, alpha=0.05, method='fdr_bh')
        model_fdrpvalue.append(each_pvalues_adjusted)
    model_fdrpvalue_df = pd.DataFrame(model_fdrpvalue)
    model_fdrpvalue_df.columns = model_pvalue_df.index
    model_fdrpvalue_df.index = model_pvalue_df.columns
    return model_coef_df, model_pvalue_df, model_error_df, model_fdrpvalue_df.transpose()


tmp=pd.DataFrame()
for features in featuresall:
    model_coef, model_pvalue, model_error, adjpvalue = airpollutionreg(traindata, features, ylist)

    modelcoefdf = pd.DataFrame(model_coef.reset_index())
    modelpvaluedf = pd.DataFrame(model_pvalue.reset_index())
    modelerrordf = pd.DataFrame(model_error.reset_index())
    adjpvaluedf = pd.DataFrame(adjpvalue.reset_index())
    # print("modelcoefdf:", modelcoefdf.columns)

    modelcoefdf.to_csv(save_path + 'PM25modelcoef_%s.csv' %features[-1], index=False)
    modelpvaluedf.to_csv(save_path + 'PM25modelpvalue_%s.csv' %features[-1], index=False)
    modelerrordf.to_csv(save_path + 'PM25modelerror_%s.csv' %features[-1], index=False)
    adjpvaluedf.to_csv(save_path + 'PM25adjustpvalue_%s.csv' %features[-1], index=False)

    analysis = modelcoefdf[['index', 'PM25']].merge(modelpvaluedf[['index', 'PM25']],on=['index'])\
                                                  .merge(adjpvaluedf[['index', 'PM25']], on=['index'])
    analysis.columns = ['index', '%s_coeff' % features[-1], '%s_pvalue' % features[-1], '%s_adjpvalue' % features[-1]]

    analysis = analysis.merge(modelerrordf[['index', 'PM25']], on=['index'])
    analysis.rename(columns={'PM25': '%s_stderror' % features[-1]}, inplace=True)

    analysis.sort_values(['%s_adjpvalue' % features[-1]], ascending=True).to_csv(
        save_path + 'PM25_%sresult.csv' % features[-1], index=False)

    condition = set(analysis[analysis['%s_adjpvalue' % features[-1]] < 0.05]['index'])
    print("condition:", condition, len(condition))


#     analysis.columns = ['index', 'Shannon_coeff', 'Shannon_pvalue', 'Shannon_adjpvalue']
#     analysis = analysis.merge(modelerrordf[['index', 'PM25']], on=['index'])
#     analysis.rename(columns={'PM25': 'Shannon_error'}, inplace=True)
#     tmp = pd.concat([tmp, analysis], axis=0)
# tmp['index'] = air_col
# tmp.sort_values(['Shannon_adjpvalue'], ascending=True).to_csv('result/PM25_result_shannon.csv', index=False)


