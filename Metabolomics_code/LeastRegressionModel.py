import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

data = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\Metabolomics.csv')
data.shape

ylist = [f"UPDASV{i:03}" for i in range(1, 486)]

xlist = list(set(data.columns) - set(ylist) - set(['sampleID']))
print("len ylist:", len(ylist))
print("len xlist:", len(xlist))


coefresult={}
pvalueresult={}
for j in tqdm(ylist):
    yi = data[j]
    metacoeff=[]
    metapavalue=[]
    for i in xlist:
        xi = sm.add_constant(data[i])
        modeli = sm.OLS(yi, xi).fit()
        intercept_i, coeff_i = modeli.params
        intercept_pvaluei, pvaluei = modeli.pvalues
        metacoeff.append(coeff_i)
        metapavalue.append(pvaluei)

    coefresult[j] = metacoeff
    pvalueresult[j] = metapavalue


coefdf = pd.DataFrame(coefresult)
pvaldf = pd.DataFrame(pvalueresult)

print("coefdf shape:", coefdf.shape)
print("pvaldf shape:", pvaldf.shape)

coefdf['metabolomicID'] = xlist
pvaldf['metabolomicID'] = xlist
coefdf.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\modelresultforcoeff.csv', index=False)
pvaldf.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\modelresultforpvalue.csv', index=False)

pvaldf.head()

adjusted_p_values = []
for mi in [f"UPDASV{i:03}" for i in range(1, 486)]:
    # Example p-values
    pvi = pvaldf[mi]
    # Apply Bonferroni correction
    adjusted_p_values.append(multipletests(pvi, method='bonferroni')[1])

# Print adjusted p-values
print(pd.DataFrame(adjusted_p_values))

adjusted_p_valuesdf = pd.DataFrame(adjusted_p_values).T
adjusted_p_valuesdf.columns=[f"UPDASV{i:03}" for i in range(1, 486)]
adjusted_p_valuesdf['metabolomicID'] = xlist
adjusted_p_valuesdf.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\modelresultforadjpv.csv', index=False)

sig_micro=[]
for fi in [f"UPDASV{i:03}" for i in range(1, 486)]:
    if sum(adjusted_p_valuesdf[fi]<0.05)>=1:
        sig_micro.append(fi)

df_filter = adjusted_p_valuesdf[['metabolomicID'] + sig_micro]
print(df_filter.shape)
df_filtered = df_filter[~(df_filter == 1).all(axis=1)]
print(df_filtered.shape)

df_filtered.to_csv(r'C:\Users\yuqingw1\Workfolder\result\Metobolomics\modelresultforadjpv_forsignificant.csv', index=False)
