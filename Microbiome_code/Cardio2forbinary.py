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
from statsmodels.stats.outliers_influence import variance_inflation_factor
import statsmodels.api as sm
import re
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
save_path = 'result/202506new/'

trans = pd.read_csv(path + 'microbiome_allvariable_806_transformation.csv')
# shannon = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\processing_4_microbiome_new\microbiome_combine_806row_485col_shannon.csv')
# beta = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\processing_4_microbiome_new\microbiome_combine_806row_bray_curtis.csv')
# shannon['sampleID'] = shannon['sampleID'].astype(int).astype(str)
# beta['sampleID'] = beta['sampleID'].astype(int).astype(str)

trans = trans[~pd.isnull(trans['sampleID'])]
trans['sampleID'] = trans['sampleID'].astype(int).astype(str)
print(trans.shape)



data = trans.copy()
# # data = trans.merge(shannon, on=['sampleID'])
# data = trans.merge(beta, on=['sampleID'])
print(data.shape)
data = data[~ (pd.isnull(data['TNFalpha_gp']) | pd.isnull(data['Troponin_gp']) | pd.isnull(data['IL6_gp']))]

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


Micro = [f"UPDASV{i:03}" for i in range(1, 486)]
# target1 = ['Ghrelin', 'Resistin', 'Leptin', 'Insulin']
target2 = ['TNFalpha_gp', 'Troponin_gp', 'IL6_gp']

for i in target2:
    print("sum nan in label_%s:"%i, sum(pd.isnull(data[i])))
    # print("sum nan in label_%s:" % i, sum(pd.isnull(diversity[i])))


features = list(cols+Micro)
print("features:", len(features))

features = list(cols+['Shannon_Diversity'])
print("features:", len(features))


features = list( set(cols) | set(data['sampleID']))
print("features:", len(features))


################################### Logistic Regression ############################################

def varsele_findcorrelation(corrMatrix, thresh=0.9):
    """
    Given a 2D numpy array, this will find highly correlated features,
    and return a list of features to remove(position of each variable)
    params:
    - df : numpy matrix
    - thresh : correlation threshold, will remove one of pairs of features with
               a correlation greater than this value
    """

    if corrMatrix.shape[0] != corrMatrix.shape[1]:
        raise Exception('corralation matrix should be square matrix')

    n = corrMatrix.shape[0]
    removeList = []

    for i in range(n - 1):
        colindex = np.where(np.abs(corrMatrix[:, i]) > thresh)[0]
        colindex = colindex[colindex > i]
        removeList = removeList + list(colindex)
    return list(set(removeList))


def varsele_findlinearcombos(mat, thresh=1e-10):
    """
    Given a numpy matrix(row is record and column is feature), this will find collinear features,
    and return a list of features to remove(position of each variable)
    params:
    - df : numpy array
    - thresh : correlation threshold, will remove one of pairs of features with
               a correlation greater than this value
    """
    Q, R = np.linalg.qr(mat)
    removeList = np.where(np.abs(R.diagonal()) <= thresh)[0]
    return list(removeList)


def LRforCardio(data, target, features, is_categoryprocess=True, microcols = Micro):

    restcols = list(set(features) - set(microcols)-set(['bmi_binned2_25.0_30.0', 'smoking_status_Former','bmi_binned2_inf_18.5']))
    print(" microcols len:", len(microcols))
    print(" restcols:", len(restcols), restcols)
    # print("check if have nan:", sum(pd.isnull(data[microcols[0]])))

    if is_categoryprocess:
        ncols = list(set(restcols) - set(['age', 'Shannon_Diversity']))
        data[ncols] = data[ncols].apply(lambda x: (x - x.mean()) / x.std(), axis=0)


    result = []
    for yi in target:
        data = data[~pd.isnull(data[yi])].reset_index(drop=True)

        y_train = data[yi]
        y_test = data[yi]

        model_coef = {}
        model_pvalue = {}
        model_error = {}
        if microcols:
            for fi in microcols:
                print("feature:", fi)
                X_train = data[restcols + [fi]]
                X_test = data[restcols + [fi]]
                dropcols=[]
                for i in restcols:
                    if len(np.unique(X_train[i]))==1:
                        dropcols.append(i)
                print("dropcols:", dropcols)
                X_train.drop(dropcols, axis=1, inplace=True)
                X_test.drop(dropcols, axis=1, inplace=True)

                ### remove colinear variable
                varsele_corr_linear_remain = np.array(list(set(restcols) - set(dropcols)))
                res = varsele_findcorrelation(np.corrcoef(X_train[X_train.columns], rowvar=False))
                varsele_corr_linear = set(varsele_corr_linear_remain[res])
                print("remove colinear var:", varsele_corr_linear)
                restcols = list(set(varsele_corr_linear_remain) - varsele_corr_linear)
                # # ### remove vif
                vif_df = pd.DataFrame(
                    [{"variable": restcols[i], "vif": variance_inflation_factor(X_train[restcols], i)} for i in
                     range(len(restcols))])
                print("vif_df:", vif_df[vif_df['vif']>8])
                restcols = list(set(restcols) - set(vif_df[vif_df['vif'] > 8]['variable']))
                print("remaining cols:", len(restcols), restcols)
                X_train = data[restcols + [fi]]
                X_test = data[restcols + [fi]]

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

            model_coef_df.to_csv(r'C:\\Users\\yuqingw1\\Workfolder\\result\\Biomarker\\biomarker_lr_coef_%s.csv' %yi)
            model_pvalue_df.to_csv(r'C:\\Users\\yuqingw1\\Workfolder\\result\\Biomarker\\biomarker_lr_pvalue_%s.csv' %yi)
            model_error_df.to_csv(r'C:\\Users\\yuqingw1\\Workfolder\\result\\Biomarker\\biomarker_lr_error_%s.csv' %yi)

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
            X_train = data[restcols]
            X_test = data[restcols]
            dropcols = []
            for i in restcols:
                if len(np.unique(X_train[i])) == 1:
                    dropcols.append(i)
            print("dropcols:", dropcols)
            X_train.drop(dropcols, axis=1, inplace=True)
            X_test.drop(dropcols, axis=1, inplace=True)

            X_train = sm.add_constant(X_train)
            smml = sm.Logit(np.array(y_train), X_train).fit()

            model_coef[yi] = list(smml.params)  # model.coef_[0]
            model_pvalue[yi] = list(smml.pvalues)  # p_values
            model_error[yi] = list(smml.bse)  # SE

            orgname = list(X_train.columns)
            model_coef_df = pd.DataFrame(smml.params).transpose()
            model_coef_df.columns = orgname
            model_pvalue_df = pd.DataFrame(smml.pvalues).transpose()
            model_pvalue_df.columns = orgname
            model_error_df = pd.DataFrame(smml.bse).transpose()
            model_error_df.columns = orgname
            # print("model_coef_df:", model_coef_df)

            model_coef_df.to_csv(r'C:\\Users\\yuqingw1\\Workfolder\\result\\Biomarker\\biomarker_lr_coef_shannon_%s.csv' %yi)
            model_pvalue_df.to_csv(r'C:\\Users\\yuqingw1\\Workfolder\\result\\Biomarker\\biomarker_lr_pvalue_shannon_%s.csv' %yi)
            model_error_df.to_csv(r'C:\\Users\\yuqingw1\\Workfolder\\result\\Biomarker\\biomarker_lr_error_shannon_%s.csv' %yi)

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


modelresult = LRforCardio(data, target2, features, is_categoryprocess=False)
print(len(modelresult))
for i in range(len(modelresult)):
    name = target2[i]
    modelresult[i].sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv(r'C:\\Users\\yuqingw1\\Workfolder\\result\\Biomarker_%s.csv' %name, index=False)





modelresult = LRforCardio(data, target2, features, is_categoryprocess=False, microcols=[])
tmp=pd.DataFrame()
for i in range(len(modelresult)):
    tmp = pd.concat([tmp, pd.DataFrame(modelresult[i])], axis=0)
tmp['index']=target2
tmp.sort_values(['adjust_Pvalue', 'pvalue'], ascending=True).to_csv(r'C:\\Users\\yuqingw1\\Workfolder\\result\\Biomarker_shannon_LR.csv', index=False)




#### beta
print(len(features))
X_train = data[list(set(features) - set(['bmi_binned2_25.0_30.0', 'smoking_status_Former','bmi_binned2_inf_18.5']))]
X_test = data[list(set(features) - set(['bmi_binned2_25.0_30.0', 'smoking_status_Former','bmi_binned2_inf_18.5']))]
X_train = sm.add_constant(X_train)
for yi in target2:
    smml = sm.Logit(np.array(data[yi]), X_train).fit()