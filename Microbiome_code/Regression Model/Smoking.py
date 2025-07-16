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


compass_deidentified = pd.read_csv('./Datasets/compass_deidentified.csv')
print(compass_deidentified.shape)
compass_deidentified = compass_deidentified[~pd.isnull(compass_deidentified['specimen_id'])]
compass_deidentified['specimen_id'] = compass_deidentified['specimen_id'].astype(int).astype(str)
print(compass_deidentified.shape)
compass_deidentified = compass_deidentified.drop_duplicates(['specimen_id'], keep='first')
compass_deidentified.loc[compass_deidentified['specimen_id']=='106503', 'age'] = (53.36 + 56.77)/2
print(compass_deidentified.shape)
compass_deidentified = compass_deidentified[~ (pd.isnull(compass_deidentified['systolic']) & pd.isnull(compass_deidentified['diastolic'])) ]
compass_deidentified = compass_deidentified.reset_index(drop=True)
print(compass_deidentified.shape)


micro = pd.read_csv('Datasets/processing_4_microbiome/microbiome_combine_trans.csv')
micro['sampleID'] = micro['sampleID'].astype(str)
print(micro.shape)


data = micro.merge(compass_deidentified, left_on=['sampleID'], right_on=['specimen_id'])
print("data.shape:", data.shape)


ylist = [f"NEWASV{i:03}" for i in range(1, 481)]
features = ['age', 'bmi', 'gender', 'smoking_status']
print("length features:", len(features))
print("length target:", len(ylist))
data['gender'] = np.where(data['gender'] == 'Male', 1, 0)
print(data['gender'].value_counts())

print(sum(pd.isnull(data['smoking_status'])))


def Smokinregression(data, target, features, is_categoryprocess=True):

    # X_train, X_test, y_train_all, y_test_all = train_test_split(data[features], data[ylist], test_size=0.1, random_state=42)
    X_train = data[features]
    X_test = data[features]
    y_train_all = data[target].apply(lambda x: np.log(x))
    y_test_all = data[target].apply(lambda x: np.log(x))

    if is_categoryprocess:
        # Bin numeric data, treating NaN as a separate category

        # cutPoint = list(np.unique(np.percentile(X_train[~pd.isnull(X_train['age'])]['age'], range(10, 100, 10))))
        # cutPoint = [35, 50, 65]
        # X_train_all['age_binned2'] = pd.cut(X_train_all['age'], bins=[-np.inf] + cutPoint + [np.inf], include_lowest=True)
        # X_train_all['age_binned2'] = np.where(pd.isnull(X_train_all['age_binned2']), 'missing', X_train_all['age_binned2'])
        # X_test_all['age_binned2'] = pd.cut(X_test_all[~pd.isnull(X_test_all['age'])]['age'], bins=[-np.inf] + cutPoint + [np.inf],
        #                                include_lowest=True)
        # X_test_all['age_binned2'] = np.where(pd.isnull(X_test_all['age_binned2']), 'missing', X_test_all['age_binned2'])

        cutPoint = [18.5, 25, 30]
        X_train['bmi_binned2'] = pd.cut(X_train[~pd.isnull(X_train['bmi'])]['bmi'], bins=[-np.inf] + cutPoint + [np.inf])
        X_train['bmi_binned2'] = np.where(pd.isnull(X_train['bmi_binned2']), 'missing', X_train['bmi_binned2'])
        X_test['bmi_binned2'] = pd.cut(X_test[~pd.isnull(X_test['bmi'])]['bmi'], bins=[-np.inf] + cutPoint + [np.inf])
        X_test['bmi_binned2'] = np.where(pd.isnull(X_test['bmi_binned2']), 'missing', X_test['bmi_binned2'])
        print(X_train['bmi_binned2'].value_counts())

        dummies = pd.get_dummies(X_train, columns=['bmi_binned2'], prefix=["bmi"], drop_first=True)
        X_train = dummies.drop(['bmi'], axis=1)
        print(X_train.shape)
        print(X_train.columns)

        X_test_dummies = pd.get_dummies(X_test, columns=['bmi_binned2'])
        for column in X_train.columns:
            if column not in X_test_dummies.columns:
                X_test_dummies[column] = 0

        X_test = X_test_dummies[X_train.columns]
        print(X_test.shape)
        print(X_test.columns)

    else:
        X_train['bmi'] = np.where(pd.isnull(X_train['bmi']), 0, X_train['bmi'])
        X_test['bmi'] = np.where(pd.isnull(X_test['bmi']), 0, X_test['bmi'])

    X_train['age'] = np.where((X_train['age'] == 0) | (pd.isnull(X_train['age'])), X_train['age'].mean(),
                              X_train['age'])
    X_test['age'] = np.where((X_test['age'] == 0) | (pd.isnull(X_test['age'])), X_train['age'].mean(), X_test['age'])


    X_train['smoking_status_current'] = np.where(X_train['smoking_status'] == 'Current', 1, 0)
    X_test['smoking_status_current'] = np.where(X_test['smoking_status'] == 'Current', 1, 0)
    X_train['smoking_status_former'] = np.where(X_train['smoking_status'] == 'Former', 1, 0)
    X_test['smoking_status_former'] = np.where(X_test['smoking_status'] == 'Former', 1, 0)
    X_train['smoking_status_missing'] = np.where(pd.isnull(X_train['smoking_status']), 1, 0)
    X_test['smoking_status_missing'] = np.where(pd.isnull(X_test['smoking_status']), 1, 0)

    X_train.drop(['smoking_status'], axis=1, inplace=True)
    X_test.drop(['smoking_status'], axis=1, inplace=True)
    print("X_train columns:", X_train.columns)
    print("X_test columns:", X_test.columns)

    microcols = [f"NEWASV{i:03}" for i in range(1, 481)]
    restcols = list(set(X_train.columns) - set([f"NEWASV{i:03}" for i in range(1, 481)]))
    print(" microcols len:", len(microcols))
    print(" restcols len:", len(restcols))

    y_train_all[microcols] = y_train_all[microcols].fillna(0)
    y_test_all[microcols] = y_test_all[microcols].fillna(0)


    result = {}
    model_coef = {}
    model_pvalue = {}
    model_error = {}

    for yi in target:
        # print("current target:", yi)
        y_train = y_train_all[yi]
        y_test = y_test_all[yi]

        # Initialize models
        X_train = sm.add_constant(X_train)
        # print("X_train2 columns:", X_train.columns)
        model = sm.OLS(y_train, X_train).fit()

        # Display results
        model_coef[yi] = list(model.params)
        model_pvalue[yi] = list(model.pvalues)
        model_error[yi] = list(model.bse)
        # print("Coefficients:", model.params)
        # print("Standard Errors:", model.bse)  # exclude intercept
        # print("p-values:", model.pvalues)


    updname = list(X_train.columns)

    model_coef_df = pd.DataFrame(model_coef).transpose()
    model_coef_df.columns = updname
    model_pvalue_df = pd.DataFrame(model_pvalue).transpose()
    model_pvalue_df.columns = updname
    model_error_df = pd.DataFrame(model_error).transpose()
    model_error_df.columns = updname
    # print("model_coef_df:", model_coef_df)

    # adjust p_value utilizing
    pval = model_pvalue_df
    pval = pval.fillna(0)
    adjpvalue={}
    for i in list(pval.columns):
        _, model_fdrpvalue, _, _ = multipletests(pval[i], alpha=0.05, method='fdr_bh')
        adjpvalue[i] = model_fdrpvalue
    adjpvalue_df = pd.DataFrame(adjpvalue)
    adjpvalue_df['index'] = model_coef_df.index
    print("adj pvalue length:", adjpvalue_df.head())

    return model_coef_df, model_pvalue_df, model_error_df, adjpvalue_df


model_coef_df,model_pvalue_df,model_error_df, model_adjpvalue_df = Smokinregression(data, ylist, features, is_categoryprocess=True)

model_coef_df.to_csv('result/model_result/smoking_on_microbiome_coef.csv')
model_pvalue_df.to_csv('result/model_result/smoking_on_microbiome_pvalue.csv')
model_error_df.to_csv('result/model_result/smoking_on_microbiome_error.csv')
model_adjpvalue_df.to_csv('result/model_result/smoking_on_microbiome_adjpvalue.csv', index=False)



analysis_former = model_coef_df['smoking_status_former'].reset_index().merge(model_pvalue_df['smoking_status_former'].reset_index(),\
                on=['index']).merge(model_error_df['smoking_status_former'].reset_index(), on=['index'])\
                .merge(model_adjpvalue_df[['index','smoking_status_former']], on=['index'])
analysis_former.columns = ['index', 'coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
analysis_former.sort_values(['adjust_Pvalue',  'pvalue']).to_csv('result/smoking_former_onmicro.csv', index=False)


analysis_current = model_coef_df['smoking_status_current'].reset_index().merge(model_pvalue_df['smoking_status_current'].reset_index(),\
                on=['index']).merge(model_error_df['smoking_status_current'].reset_index(), on=['index'])\
                .merge(model_adjpvalue_df[['index','smoking_status_current']], on=['index'])
analysis_current.columns = ['index', 'coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
analysis_current.sort_values(['adjust_Pvalue',  'pvalue']).to_csv('result/smoking_current_onmicro.csv', index=False)


analysis_missing = model_coef_df['smoking_status_missing'].reset_index().merge(model_pvalue_df['smoking_status_missing'].reset_index(),\
                on=['index']).merge(model_error_df['smoking_status_missing'].reset_index(), on=['index'])\
                .merge(model_adjpvalue_df[['index','smoking_status_missing']], on=['index'])
analysis_missing.columns = ['index', 'coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
analysis_missing.sort_values(['adjust_Pvalue',  'pvalue']).to_csv('result/smoking_missing_onmicro.csv', index=False)


#
# test = pd.read_csv('/Users/yuqingwu/UChicago_MSDS/Part-time/BSD/Datasets/processing_4_microbiome/microbiome_combine_790row_480col.csv')
# test.shape
#
# cols = [
# 'NEWASV188',
# 'NEWASV193',
# 'NEWASV262',
# 'NEWASV281',
# 'NEWASV320',
# 'NEWASV334',
# 'NEWASV357',
# 'NEWASV358',
# 'NEWASV361',
# 'NEWASV369',
# 'NEWASV371',
# 'NEWASV373',
# 'NEWASV382',
# 'NEWASV385',
# 'NEWASV387',
# 'NEWASV395',
# 'NEWASV397',
# 'NEWASV398',
# 'NEWASV409',
# 'NEWASV410',
# 'NEWASV412',
# 'NEWASV413',
# 'NEWASV416',
# 'NEWASV421',
# 'NEWASV423',
# 'NEWASV426',
# 'NEWASV430',
# 'NEWASV431',
# 'NEWASV434',
# 'NEWASV435',
# 'NEWASV438',
# 'NEWASV446',
# 'NEWASV449',
# 'NEWASV450',
# 'NEWASV458',
# 'NEWASV460',
# 'NEWASV462',
# 'NEWASV463',
# 'NEWASV466',
# 'NEWASV470',
# 'NEWASV471',
# 'NEWASV472',
# 'NEWASV475',
# ]
#
# test[cols] = test[cols].fillna(0)
# test['sampleID'] = test['sampleID'].astype(str)
# for i in cols:
#     print("num of non-zero:", sum((test[i]!=0.0)))
#     print(data[data['sampleID']== list(test[test[i]!=0]['sampleID'])[0]]['smoking_status'])
