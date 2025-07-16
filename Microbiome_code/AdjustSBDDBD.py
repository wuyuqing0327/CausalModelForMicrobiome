import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, LogisticRegression
import statsmodels.api as sm
from sklearn.model_selection import cross_val_score
from scipy import stats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import re


trans = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\microbiome_allvariable_806_transformation.csv')
shannon = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\processing_4_microbiome_new\microbiome_combine_806row_485col_shannon.csv')


print(trans['smoking_status'].value_counts())

print(trans['diabetes2'].value_counts())



trans = trans[~pd.isnull(trans['sampleID'])]
trans['sampleID'] = trans['sampleID'].astype(int).astype(str)
print(trans.shape)

shannon = shannon[~pd.isnull(shannon['sampleID'])]
shannon['sampleID'] = shannon['sampleID'].astype(int).astype(str)
print(shannon.shape)


trans = trans[~ (pd.isnull(trans['systolic']) & pd.isnull(trans['diastolic'])) ]
trans = trans.reset_index(drop=True)
print(trans.shape)

"""
Regualations for adjusting systolic and diastolic
Regarding systolic and diastolic
1) LOOK INTO systolic_1, systolic_2 & diastolic_1, diastolic_2
do turn over its value if it is needed
2) if all the value in the same value then delete
3) 

"""

data = trans.copy()
# data = trans.merge(shannon, on=['sampleID'])
### exchange systolic and diastolic
equalcon = (data['systolic']==data['diastolic'])
print(sum(equalcon))
mask1 = data['systolic'] <= data['diastolic']
print( " incorrrect SBD AND DBD 1:", sum(~equalcon & mask1))
mask3 = ( (data['systolic'] - data['diastolic']) <=15 )
print( " incorrrect SBD AND DBD 2:", sum( ~mask1  & mask3 ))


data = data[~equalcon]
# data = data[~np.in1d(data['sampleID'], ['103358', '104610', '106851', '106433', '106411' ])]
data = data[~ (pd.isnull(data['systolic']) | pd.isnull(data['diastolic']))]
print(data.shape)


cutPoint = [18.5, 25, 30]
data['bmi_binned2'] = pd.cut(data[~pd.isnull(data['bmi'])]['bmi'],bins=[-np.inf] + cutPoint + [np.inf]).astype(str)
data['bmi_binned2'] = np.where(pd.isnull(data['bmi_binned2']), 'missing', data['bmi_binned2'])
print(data['bmi_binned2'].value_counts())
data['bmi_binned2_30.0_inf'] = np.where(data['bmi_binned2']=='(30.0, inf]', 1, 0)
data['bmi_binned2_25.0_30.0'] = np.where(data['bmi_binned2']=='(25.0, 30.0]', 1, 0)
data['bmi_binned2_18.5_25.0'] = np.where(data['bmi_binned2']=='(18.5, 25.0]', 1, 0)
data['bmi_binned2_missing'] = np.where(data['bmi_binned2']=='missing', 1, 0)

print(data.shape)
data['gender_Male']=np.where(data['gender']=='Male', 1, 0)
data['smoking_status_Current']=np.where(data['smoking_status']=='Current', 1, 0)
data['smoking_status_Former']=np.where(data['smoking_status']=='Former', 1, 0)
data['smoking_status_Missing']=np.where(pd.isnull(data['smoking_status']), 1, 0)
data['diabetes2_Yes'] = np.where(data['diabetes2']=='Yes', 1, 0)
data['diabetes2_No'] = np.where(data['diabetes2']=='No', 1, 0)

cols = ['sampleID', 'age', 'gender_Male',
        'bmi_binned2_30.0_inf','bmi_binned2_25.0_30.0','bmi_binned2_18.5_25.0','bmi_binned2_missing',
        'smoking_status_Current','smoking_status_Former','smoking_status_Missing',
        'diabetes2_Yes','diabetes2_No',
        'systolic', 'diastolic', 'HBP']

microcols = [f"UPDASV{i:03}" for i in range(1, 494)]
data = data.drop_duplicates(cols)
data = data[cols+microcols]
# data = data[cols + ['Shannon_Diversity']]
data = data.reset_index(drop=True)
print(data.shape)


target = ['systolic', 'diastolic']
""" 
HIGH BLOOD: 
SBP > 130 / DBP>80 then 1, else 0
SBP >= 140 / DBP>=90 then 1, else 0
"""
target2=['HBP']

features = list(set(data.columns) - set(target)-set(target2) - set(['sampleID']))
print(len(features))



########################################################################################################
### regression model

def Regressionforblood(data, target, features, microcols = [f"UPDASV{i:03}" for i in range(1, 494)]):

    # X_train, X_test, y_train_all, y_test_all = train_test_split(data[features], data[ylist], test_size=0.1, random_state=42)
    X_train_all = data[features]
    X_test_all = data[features]
    y_train_all = data[target]
    y_test_all = data[target]

    X_train_all['age'] = np.where((X_train_all['age'] == 0) | (pd.isnull(X_train_all['age'])), X_train_all['age'].mean(),
                              X_train_all['age'])
    X_test_all['age'] = np.where((X_test_all['age'] == 0) | (pd.isnull(X_test_all['age'])), X_train_all['age'].mean(), X_test_all['age'])

    restcols = list(set(X_train_all.columns) - set(microcols))
    print(" microcols len:", len(microcols))
    print(" restcols len:", len(restcols))

    result = []
    for yi in target:
        # print("current target:", yi)
        y_train = y_train_all[yi]
        y_test = y_test_all[yi]

        model_coef = {}
        model_pvalue = {}
        model_error = {}
        for fi in microcols:
            X_train = X_train_all[restcols + [fi]]
            X_test = X_test_all[restcols + [fi]]

            dropcols = []
            for i in X_train.columns:
                if len(np.unique(X_train[i])) == 1:
                    dropcols.append(i)
            print("dropcols:", dropcols)
            X_train.drop(dropcols, axis=1, inplace=True)
            X_test.drop(dropcols, axis=1, inplace=True)
            # Initialize models
            X_train = sm.add_constant(X_train)
            model = sm.OLS(y_train, X_train).fit()
            # Display results
            model_coef[fi] = list(model.params)
            model_pvalue[fi] = list(model.pvalues)
            model_error[fi] = list(model.bse)

        pattern = r'^UPDASV\d+$'
        replacement = 'microbiome'
        orgname = list(X_train.columns)
        updname = [replacement if re.match(pattern, item) else item for item in orgname]

        model_coef_df = pd.DataFrame(model_coef).transpose()
        model_coef_df.columns = updname
        model_pvalue_df = pd.DataFrame(model_pvalue).transpose()
        model_pvalue_df.columns = updname
        model_error_df = pd.DataFrame(model_error).transpose()
        model_error_df.columns = updname

        model_coef_df.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\linear_moel_coef_%s_shannon.csv' %yi)
        model_pvalue_df.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\linear_moel_pvalue_%s_shannon.csv' %yi)
        model_error_df.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\linear_moel_error_%s_shannon.csv' %yi)

        # adjust p_value utilizing
        pval = list(model_pvalue_df['microbiome'])
        print("length pval:", len(pval))
        pval = [0 if (x is None) | (np.isnan(x)) else x for x in pval]
        _, model_fdrpvalue, _, _ = multipletests(pval, alpha=0.05, method='fdr_bh')
        print("adj pvalue length:", len(model_fdrpvalue))
        # pd.DataFrame(model_fdrpvalue, columns=[restcols + ['microbiome']]).to_csv('result/linear_moel_adjpvalue_%s.csv' %yi, index=False)

        analysis = model_coef_df['microbiome'].reset_index().merge(model_pvalue_df['microbiome'].reset_index(), on=['index']).merge(
                   model_error_df['microbiome'].reset_index(), on=['index'])
        analysis['adjust_Pvalue'] = model_fdrpvalue
        analysis.columns = ['index', 'coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
        # print(analysis)
        result.append(analysis)
    return result


modelresult = Regressionforblood(data, target, features)


modelresult[0].sort_values(['adjust_Pvalue'], ascending=True).to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\systolicresult.csv', index=False)
modelresult[1].sort_values(['adjust_Pvalue'], ascending=True).to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\diastolicresult.csv', index=False)


########################################### Binary Classification ########################################################################



#### Logistic Regression
def LRforblood(data, target, features, is_categoryprocess=True, microcols = [f"UPDASV{i:03}" for i in range(1, 494)]):

    # X_train, X_test, y_train_all, y_test_all = train_test_split(data[features], data[ylist], test_size=0.1, random_state=42)
    X_train_all = data[features]
    X_test_all = data[features]
    y_train_all = data[target]
    y_test_all = data[target]

    X_train_all['age'] = np.where((X_train_all['age'] == 0) | (pd.isnull(X_train_all['age'])), X_train_all['age'].mean(),
                              X_train_all['age'])
    X_test_all['age'] = np.where((X_test_all['age'] == 0) | (pd.isnull(X_test_all['age'])), X_train_all['age'].mean(), X_test_all['age'])

    restcols = list(set(X_train_all.columns) - set(microcols))
    # print(" microcols len:", len(microcols))
    # print(" restcols len:", len(restcols))


    if is_categoryprocess:
        ncols = list(set(restcols) - set(['age', 'Shannon_Diversity']))
        X_train_all[ncols] = X_train_all[ncols].apply(lambda x: (x - x.mean()) / x.std(), axis=0)
        X_test_all[ncols] = (X_test_all[ncols] - X_train_all[ncols].mean()) / X_train_all[ncols].std()

    result = []
    for yi in target:
        # print("current target:", yi)
        y_train = y_train_all[yi]
        y_test = y_test_all[yi]

        model_coef = {}
        model_pvalue = {}
        model_error = {}
        for fi in microcols:
            print("feature:", fi)
            X_train = X_train_all[restcols + [fi]]
            X_test = X_test_all[restcols + [fi]]

            dropcols = []
            for i in X_train.columns:
                if len(np.unique(X_train[i])) == 1:
                    dropcols.append(i)
            print("dropcols:", dropcols)
            X_train.drop(dropcols, axis=1, inplace=True)
            X_test.drop(dropcols, axis=1, inplace=True)

            X_train = sm.add_constant(X_train)
            smml = sm.Logit(np.array(y_train), X_train).fit()
            # Display results
            model_coef[fi] = list(smml.params) #model.coef_[0]
            model_pvalue[fi] = list(smml.pvalues) #p_values
            model_error[fi] = list(smml.bse) #SE

        pattern = r'^UPDASV\d+$'
        replacement = 'microbiome'
        orgname = list(X_train.columns)
        updname = [replacement if re.match(pattern, item) else item for item in orgname]

        model_coef_df = pd.DataFrame(model_coef).transpose()
        model_coef_df.columns = updname
        model_pvalue_df = pd.DataFrame(model_pvalue).transpose()
        model_pvalue_df.columns = updname
        model_error_df = pd.DataFrame(model_error).transpose()
        model_error_df.columns = updname

        model_coef_df.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\lr_moel_coef_%s_shannon.csv' %yi)
        model_pvalue_df.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\lr_moel_pvalue_%s_shannon.csv' %yi)
        model_error_df.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\lr_moel_error_%s_shannon.csv' %yi)

#        # adjust p_value utilizing
#         pval = list(model_pvalue_df['microbiome'])
#         print("length pval:", len(pval))
#         pval = [0 if (x is None) | (np.isnan(x)) else x for x in pval]
#         _, model_fdrpvalue, _, _ = multipletests(pval, alpha=0.05, method='fdr_bh')
#         print("adj pvalue length:", len(model_fdrpvalue))
#
#         analysis = model_coef_df['microbiome'].reset_index().merge(model_pvalue_df['microbiome'].reset_index(), on=['index']).merge(
#                    model_error_df['microbiome'].reset_index(), on=['index'])
#         analysis['adjust_Pvalue'] = model_fdrpvalue
#         analysis.columns = ['index','coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
#         # print(analysis)
#         result.append(analysis)
    return result


modelresult = LRforblood(data, target2, features, is_categoryprocess=True, microcols=['Shannon_Diversity'])


modelresult[0].sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv(r'C:\Users\yuqingw1\Workfolder\result\Hypertension\HighBP140_90.csv', index=False)





#
# sum(data['NEWASV141']!=0)
#
# X_train = pd.read_csv('result/test.csv')
# y_train = data['HBP1']
# X_train2 = sm.add_constant(X_train)
# smml = sm.Logit(y_train, X_train).fit()
#
# variances = X_train2.var()
# zero_variance_cols = variances[variances == 0].index.tolist()
# print("Columns with zero variance:", zero_variance_cols)
# corr_matrix = X_train2.corr().abs()
# high_corr_var=np.where(corr_matrix>0.8)
# high_corr_var=[(corr_matrix.columns[x],corr_matrix.columns[y]) for x,y in zip(*high_corr_var) if x!=y and x<y]
# print("Highly correlated variables:", high_corr_var)