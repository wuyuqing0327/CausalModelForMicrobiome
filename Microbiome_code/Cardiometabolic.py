import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.model_selection import cross_val_score
from scipy import stats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
from sklearn.model_selection import cross_val_score
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import re


path = 'Datasets/'
save_path = 'result/202506new/microbiome_filtermissingover50/'

trans = pd.read_csv( path + 'microbiome_allvariable_806_transformation.csv')
# shannon = pd.read_csv( path + 'microbiome_combine_806row_485col_shannon.csv')

trans = trans[~pd.isnull(trans['sampleID'])]
trans['sampleID'] = trans['sampleID'].astype(int).astype(str)
# shannon['sampleID'] = shannon['sampleID'].astype(int).astype(str)
print(trans.shape)



data = trans.copy()
# data = trans.merge(shannon, on=['sampleID'])
print(data.shape)
data = data[~ (pd.isnull(data['Ghrelin']) | pd.isnull(data['Resistin']) | pd.isnull(data['Insulin']))]
# data = data[~ (pd.isnull(data['MCP1']) | pd.isnull(data['Cpeptide']) | pd.isnull(data['CKMB']))]
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

cols = ['age', 'gender_Male',
        'bmi_binned2_30.0_inf','bmi_binned2_25.0_30.0','bmi_binned2_inf_18.5',
        'smoking_status_Current','smoking_status_Former',
        'diabetes2_Yes']


# Micro = [f"UPDASV{i:03}" for i in range(1, 486)]
Micro = ['UPDASV030', 'UPDASV051', 'UPDASV053', 'UPDASV065', 'UPDASV073', 'UPDASV139', 'UPDASV169', 'UPDASV184', 'UPDASV186', 'UPDASV187', 'UPDASV190', 'UPDASV197', 'UPDASV224', 'UPDASV236', 'UPDASV250', 'UPDASV264', 'UPDASV270', 'UPDASV285', 'UPDASV291', 'UPDASV302', 'UPDASV310', 'UPDASV345', 'UPDASV358', 'UPDASV359', 'UPDASV372', 'UPDASV374', 'UPDASV385', 'UPDASV387', 'UPDASV394', 'UPDASV402', 'UPDASV413', 'UPDASV424', 'UPDASV429', 'UPDASV430', 'UPDASV434', 'UPDASV442', 'UPDASV467', 'UPDASV480']
#Micro_filtermissing75 = ['UPDASV024', 'UPDASV030', 'UPDASV039', 'UPDASV051', 'UPDASV053', 'UPDASV060', 'UPDASV061', 'UPDASV065', 'UPDASV073', 'UPDASV139', 'UPDASV150', 'UPDASV153', 'UPDASV155', 'UPDASV169', 'UPDASV176', 'UPDASV184', 'UPDASV186', 'UPDASV187', 'UPDASV190', 'UPDASV197', 'UPDASV198', 'UPDASV207', 'UPDASV224', 'UPDASV234', 'UPDASV236', 'UPDASV243', 'UPDASV250', 'UPDASV264', 'UPDASV268', 'UPDASV270', 'UPDASV273', 'UPDASV285', 'UPDASV291', 'UPDASV296', 'UPDASV302', 'UPDASV310', 'UPDASV345', 'UPDASV350', 'UPDASV358', 'UPDASV359', 'UPDASV361', 'UPDASV366', 'UPDASV372', 'UPDASV374', 'UPDASV376', 'UPDASV385', 'UPDASV387', 'UPDASV388', 'UPDASV391', 'UPDASV393', 'UPDASV394', 'UPDASV398', 'UPDASV401', 'UPDASV402', 'UPDASV413', 'UPDASV424', 'UPDASV429', 'UPDASV430', 'UPDASV434', 'UPDASV442', 'UPDASV444', 'UPDASV447', 'UPDASV450', 'UPDASV467', 'UPDASV480', 'UPDASV485']


Micro = list(set(Micro) - set(['UPDASV186'])) # remove 'UPDASV186' because (具体不太清楚，可能是这个microbiome与某个microbiome重复了）
target1 = ['Ghrelin', 'Resistin', 'Insulin'] # 'Leptin',
# target2 = ['TNFalpha_gp', 'Troponin_gp', 'IL6_gp']
# target3 = ['MCP1', 'Cpeptide', 'CKMB']

for i in target1:
    print("sum nan in label_%s:"%i, sum(pd.isnull(data[i])))
    # print("sum nan in label_%s:" % i, sum(pd.isnull(diversity[i])))

features = list(cols+Micro)
print("features:", len(features))


def RegressionforCardio(data, target, features, microcols=Micro):

    restcols = list(set(features) - set(microcols))
    print(" microcols len:", len(microcols))
    print(" restcols len:", restcols)


    result = []
    for yi in target:
        data = data[~pd.isnull(data[yi])].reset_index(drop=True)
        data[yi] = data[yi].apply(lambda x: np.sign(x) * np.log(np.abs(x)))
        data[yi] = (data[yi] - data[yi].mean()) / data[yi].std()

        y_train = data[yi]
        y_test = data[yi]

        model_coef = {}
        model_pvalue = {}
        model_error = {}
        if microcols:
            for fi in microcols:
                X_train = data[restcols + [fi]]
                X_test = data[restcols + [fi]]

                dropcols=[]
                for i in restcols:
                    if len(np.unique(X_train[i]))==1:
                        dropcols.append(i)
                print("dropcols:", dropcols)
                X_train.drop(dropcols, axis=1, inplace=True)
                X_test.drop(dropcols, axis=1, inplace=True)

                # Initialize models
                X_train = sm.add_constant(X_train)
                # print("X_train2 columns:", X_train.columns)
                model = sm.OLS(y_train, X_train).fit()

                # Display results
                model_coef[fi] = list(model.params)
                model_pvalue[fi] = list(model.pvalues)
                model_error[fi] = list(model.bse)
                # print("Coefficients:", model.params)
                # print("Standard Errors:", model.bse)  # exclude intercept
                # print("p-values:", model.pvalues)

            pattern = r'^UPDASV\d+'
            replacement = 'microbiome'
            orgname = list(X_train.columns)
            updname = [replacement if re.match(pattern, item) else item for item in orgname]

            model_coef_df = pd.DataFrame(model_coef).transpose()
            model_coef_df.columns = updname
            model_pvalue_df = pd.DataFrame(model_pvalue).transpose()
            model_pvalue_df.columns = updname
            model_error_df = pd.DataFrame(model_error).transpose()
            model_error_df.columns = updname
            # print("model_coef_df:", model_coef_df)

            model_coef_df.to_csv(save_path + 'biomarker_linear_coef_%s.csv' %yi)
            model_pvalue_df.to_csv(save_path + 'biomarker_linear_pvalue_%s.csv' %yi)
            model_error_df.to_csv(save_path + 'biomarker_linear_error_%s.csv' %yi)

            # adjust p_value utilizing
            pval = list(model_pvalue_df['microbiome'])
            print("length pval:", len(pval))
            pval = [0 if (x is None) | (np.isnan(x)) else x for x in pval]
            _, model_fdrpvalue, _, _ = multipletests(pval, alpha=0.05, method='fdr_bh')
            print("adj pvalue length:", len(model_fdrpvalue))

            analysis = model_coef_df['microbiome'].reset_index().merge(model_pvalue_df['microbiome'].reset_index(), on=['index']).merge(
                       model_error_df['microbiome'].reset_index(), on=['index'])
            analysis['adjust_Pvalue'] = model_fdrpvalue
            analysis.columns = ['index','coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
            # print(analysis)
            result.append(analysis)
        else:
            print("Shannon feature length:", len(restcols), restcols)
            X_train = data[restcols]
            X_test = data[restcols]
            dropcols = []
            for i in restcols:
                if len(np.unique(X_train[i])) == 1:
                    dropcols.append(i)
            print("dropcols:", dropcols)
            X_train.drop(dropcols, axis=1, inplace=True)
            X_test.drop(dropcols, axis=1, inplace=True)

            # models
            X_train = sm.add_constant(X_train)
            model = sm.OLS(y_train, X_train).fit()

            model_coef[yi] = list(model.params)
            model_pvalue[yi] = list(model.pvalues)
            model_error[yi] = list(model.bse)

            orgname = list(X_train.columns)
            model_coef_df = pd.DataFrame(model.params).transpose()
            model_coef_df.columns = orgname
            model_pvalue_df = pd.DataFrame(model.pvalues).transpose()
            model_pvalue_df.columns = orgname
            model_error_df = pd.DataFrame(model.bse).transpose()
            model_error_df.columns = orgname
            # print("model_coef_df:", model_coef_df)

            model_coef_df.to_csv(save_path + 'biomarker_linear_coef_shannon_%s.csv' %yi)
            model_pvalue_df.to_csv(save_path + 'biomarker_linear_pvalue_shannon_%s.csv' %yi)
            model_error_df.to_csv(save_path + 'biomarker_linear_error_shannon_%s.csv' %yi)

            # adjust p_value utilizing
            pval = list(model_pvalue_df['Shannon_Diversity'])
            print("length pval:", len(pval))
            pval = [0 if (x is None) | (np.isnan(x)) else x for x in pval]
            _, model_fdrpvalue, _, _ = multipletests(pval, alpha=0.05, method='fdr_bh')
            print("adj pvalue length:", len(model_fdrpvalue))

            analysis = model_coef_df['Shannon_Diversity'].reset_index().merge(model_pvalue_df['Shannon_Diversity'].reset_index(), on=['index']).merge(
                       model_error_df['Shannon_Diversity'].reset_index(), on=['index'])
            analysis['adjust_Pvalue'] = model_fdrpvalue
            analysis.columns = ['index','coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
            # print(analysis)
            result.append(analysis)

    return result


modelresult = RegressionforCardio(data, target1, features)
for i in range(len(modelresult)):
    name = target1[i]
    modelresult[i].sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv(save_path + 'Biomarker_%s.csv' %name, index=False)


# modelresult = RegressionforCardio(data, target3, features)
# for i in range(len(modelresult)):
#     name = target3[i]
#     modelresult[i].sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv(r'C:\Users\yuqingw1\Workfolder\result\Biomarker_%s.csv' %name, index=False)
#

# modelresult = RegressionforCardio(data, target1, features,  microcols=[])
# tmp=pd.DataFrame()
# for i in range(len(modelresult)):
#     tmp = pd.concat([tmp, pd.DataFrame(modelresult[i])], axis=0)
# tmp['index']=target1
# tmp.sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv(r'C:\Users\yuqingw1\Workfolder\result\Biomarker_shannon.csv', index=False)
#
#
# features = list(cols + ['Shannon_Diversity'])
# print(len(features))
# modelresult = RegressionforCardio(data, target1+target3, features, microcols=[])
# tmp=pd.DataFrame()
# for i in range(len(modelresult)):
#     tmp = pd.concat([tmp, pd.DataFrame(modelresult[i])], axis=0)
# tmp['index']=target1+target3
# tmp.sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv(r'C:\Users\yuqingw1\Workfolder\result\Biomarker_shannon_Regression.csv', index=False)
#

#################################### Logistic Regression ############################################
#
# def varsele_findcorrelation(corrMatrix, thresh=0.9):
#     """
#     Given a 2D numpy array, this will find highly correlated features,
#     and return a list of features to remove(position of each variable)
#     params:
#     - df : numpy matrix
#     - thresh : correlation threshold, will remove one of pairs of features with
#                a correlation greater than this value
#     """
#
#     if corrMatrix.shape[0] != corrMatrix.shape[1]:
#         raise Exception('corralation matrix should be square matrix')
#
#     n = corrMatrix.shape[0]
#     removeList = []
#
#     for i in range(n - 1):
#         colindex = np.where(np.abs(corrMatrix[:, i]) > thresh)[0]
#         colindex = colindex[colindex > i]
#         removeList = removeList + list(colindex)
#     return list(set(removeList))
#
#
# def varsele_findlinearcombos(mat, thresh=1e-10):
#     """
#     Given a numpy matrix(row is record and column is feature), this will find collinear features,
#     and return a list of features to remove(position of each variable)
#     params:
#     - df : numpy array
#     - thresh : correlation threshold, will remove one of pairs of features with
#                a correlation greater than this value
#     """
#     Q, R = np.linalg.qr(mat)
#     removeList = np.where(np.abs(R.diagonal()) <= thresh)[0]
#     return list(removeList)
#
#
# def LRforCardio(data, target, features, is_categoryprocess=True, microcols = [f"NEWASV{i:03}" for i in range(1, 481)]):
#
#     restcols = list(set(features) - set(microcols))
#     print(" microcols len:", len(microcols))
#     print(" restcols:", len(restcols), restcols)
#     # print("check if have nan:", sum(pd.isnull(data[microcols[0]])))
#
#     if is_categoryprocess:
#         ncols = list(set(restcols) - set(['age', 'Shannon_Diversity']))
#         data[ncols] = data[ncols].apply(lambda x: (x - x.mean()) / x.std(), axis=0)
#
#
#     result = []
#     for yi in target:
#         data = data[~pd.isnull(data[yi])].reset_index(drop=True)
#
#         y_train = data[yi]
#         y_test = data[yi]
#
#         model_coef = {}
#         model_pvalue = {}
#         model_error = {}
#         if microcols:
#             for fi in microcols:
#                 print("feature:", fi)
#                 X_train = data[restcols + [fi]]
#                 X_test = data[restcols + [fi]]
#                 dropcols=[]
#                 for i in restcols:
#                     if len(np.unique(X_train[i]))==1:
#                         dropcols.append(i)
#                 print("dropcols:", dropcols)
#                 X_train.drop(dropcols, axis=1, inplace=True)
#                 X_test.drop(dropcols, axis=1, inplace=True)
#
#                 ### remove colinear variable
#                 varsele_corr_linear_remain = np.array(list(set(restcols) - set(dropcols)))
#                 res = varsele_findcorrelation(np.corrcoef(X_train[X_train.columns], rowvar=False))
#                 varsele_corr_linear = set(varsele_corr_linear_remain[res])
#                 print("remove colinear var:", varsele_corr_linear)
#                 restcols = list(set(varsele_corr_linear_remain) - varsele_corr_linear)
#                 # # ### remove vif
#                 vif_df = pd.DataFrame(
#                     [{"variable": restcols[i], "vif": variance_inflation_factor(X_train[restcols], i)} for i in
#                      range(len(restcols))])
#                 print("vif_df:", vif_df[vif_df['vif']>8])
#                 restcols = list(set(restcols) - set(vif_df[vif_df['vif'] > 8]['variable']))
#                 X_train = data[restcols + [fi]]
#                 X_test = data[restcols + [fi]]
#
#                 # print("X_train.columns:", X_train.columns)
#                 # # Initialize models
#                 # model = LogisticRegression(penalty='l2', C=0.01, solver='lbfgs', max_iter=1000, random_state=42, tol=0.0001)
#                 # # print("************* trainning moddel ***************")
#                 # model.fit(X_train, y_train)
#                 # # Predictions probabilities
#                 # y_pred_prob = model.predict_proba(X_test)[:, 1]
#                 #
#                 # # Add intercept to X_test
#                 # X_with_intercept = np.c_[np.ones((X_test.shape[0], 1)), X_test]
#                 # # W matrix
#                 # W = np.diag(y_pred_prob * (1 - y_pred_prob))
#                 #
#                 #
#                 # # Variance-covariance matrix of the coefficients
#                 # V = np.linalg.inv(X_with_intercept.T @ W @ X_with_intercept)
#                 # # Standard errors of the coefficients
#                 # SE = np.sqrt(np.diag(V))
#                 # # Coefficients from the model
#                 # beta = np.concatenate(([model.intercept_[0]], model.coef_[0]))
#                 # # t-statistics
#                 # t_stats = beta / SE
#                 # # Degrees of freedom
#                 # df = len(y_test) - len(beta)
#                 # # p-values
#                 # p_values = [2 * (1 - stats.t.cdf(np.abs(t_stat), df)) for t_stat in t_stats]
#                 # # print("Standard Errors:", SE)
#                 # # print("p-Values:", p_values)
#                 # # print("coefficients:", model.coef_)
#                 # # print("Intercept:", model.intercept_)
#
#                 X_train = sm.add_constant(X_train)
#                 smml = sm.Logit(np.array(y_train), X_train).fit()
#                 # Coefficients
#                 # print("Coefficients:", smml.params)
#                 # # Standard errors
#                 # print("Standard Errors:", smml.bse)
#                 # P-values
#                 # print("P-values:", smml.pvalues)
#
#                 # Display results
#                 model_coef[fi] = list(smml.params) #model.coef_[0]
#                 model_pvalue[fi] = list(smml.pvalues) #p_values
#                 model_error[fi] = list(smml.bse) #SE
#
#             pattern = r'^NEWASV\d+$'
#             replacement = 'microbiome'
#             orgname = list(X_train.columns)
#             updname = [replacement if re.match(pattern, item) else item for item in orgname]
#
#             model_coef_df = pd.DataFrame(model_coef).transpose()
#             model_coef_df.columns = updname
#             model_pvalue_df = pd.DataFrame(model_pvalue).transpose()
#             model_pvalue_df.columns = updname
#             model_error_df = pd.DataFrame(model_error).transpose()
#             model_error_df.columns = updname
#
#             model_coef_df.to_csv('result/model_result/cardiometob_lr_coef_%s.csv' %yi)
#             model_pvalue_df.to_csv('result/model_result/cardiometob_lr_pvalue_%s.csv' %yi)
#             model_error_df.to_csv('result/model_result/cardiometob_lr_error_%s.csv' %yi)
#
#             # adjust p_value utilizing
#             pval = list(model_pvalue_df['microbiome'])
#             print("length pval:", len(pval))
#             pval = [0 if (x is None) | (np.isnan(x)) else x for x in pval]
#             _, model_fdrpvalue, _, _ = multipletests(pval, alpha=0.05, method='fdr_bh')
#             print("adj pvalue length:", len(model_fdrpvalue))
#
#             analysis = model_coef_df['microbiome'].reset_index().merge(model_pvalue_df['microbiome'].reset_index(), on=['index']).merge(
#                        model_error_df['microbiome'].reset_index(), on=['index'])
#             analysis['adjust_Pvalue'] = model_fdrpvalue
#             analysis.columns = ['index','coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
#             # print(analysis)
#             result.append(analysis)
#         else:
#             X_train = data[restcols]
#             X_test = data[restcols]
#             dropcols = []
#             for i in restcols:
#                 if len(np.unique(X_train[i])) == 1:
#                     dropcols.append(i)
#             print("dropcols:", dropcols)
#             X_train.drop(dropcols, axis=1, inplace=True)
#             X_test.drop(dropcols, axis=1, inplace=True)
#
#             X_train = sm.add_constant(X_train)
#             smml = sm.Logit(np.array(y_train), X_train).fit()
#
#             model_coef[yi] = list(smml.params)  # model.coef_[0]
#             model_pvalue[yi] = list(smml.pvalues)  # p_values
#             model_error[yi] = list(smml.bse)  # SE
#
#             orgname = list(X_train.columns)
#             model_coef_df = pd.DataFrame(smml.params).transpose()
#             model_coef_df.columns = orgname
#             model_pvalue_df = pd.DataFrame(smml.pvalues).transpose()
#             model_pvalue_df.columns = orgname
#             model_error_df = pd.DataFrame(smml.bse).transpose()
#             model_error_df.columns = orgname
#             # print("model_coef_df:", model_coef_df)
#
#             model_coef_df.to_csv('result/model_result/cardiometob_lr_coef_shannon_%s.csv' %yi)
#             model_pvalue_df.to_csv('result/model_result/cardiometob_lr_pvalue_shannon_%s.csv' %yi)
#             model_error_df.to_csv('result/model_result/cardiometob_lr_error_shannon_%s.csv' %yi)
#
#             # adjust p_value utilizing
#             pval = list(model_pvalue_df['Shannon_Diversity'])
#             print("length pval:", len(pval))
#             pval = [0 if (x is None) | (np.isnan(x)) else x for x in pval]
#             _, model_fdrpvalue, _, _ = multipletests(pval, alpha=0.05, method='fdr_bh')
#             print("adj pvalue length:", len(model_fdrpvalue))
#
#             analysis = model_coef_df['Shannon_Diversity'].reset_index().merge(model_pvalue_df['Shannon_Diversity'].reset_index(), on=['index']).merge(
#                        model_error_df['Shannon_Diversity'].reset_index(), on=['index'])
#             analysis['adjust_Pvalue'] = model_fdrpvalue
#             analysis.columns = ['index','coeff', 'pvalue', 'std_error', 'adjust_Pvalue']
#             # print(analysis)
#             result.append(analysis)
#
#     return result
#
#
# # modelresult = LRforCardio(data, target2, features, is_categoryprocess=True)
# # print(len(modelresult))
# # for i in range(len(modelresult)):
# #     name = target2[i]
# #     modelresult[i].sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv('result/Cardiometobolic_%s.csv' %name, index=False)
# #
#
# modelresult = LRforCardio(data, target2, features, is_categoryprocess=True, microcols=[])
# tmp=pd.DataFrame()
# for i in range(len(modelresult)):
#     tmp = pd.concat([tmp, pd.DataFrame(modelresult[i])], axis=0)
# tmp['index']=target2
# tmp.sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv('result/Cardiometobolic_shannon2.csv', index=False)
