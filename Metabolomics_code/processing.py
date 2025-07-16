import numpy as np
import pandas as pd

namepdf = pd.read_excel(r'C:\Users\yuqingw1\Workfolder\Data\Metabolomics\MSC1445_SampleID.xlsx')
namepdf.loc[namepdf['Number']==520, 'Sample_ID'] = 'S061903715'

namelist = pd.read_csv(
    r'C:\Users\yuqingw1\Workfolder\Data\Metabolomics\COMPASS_Sample_Withdrawal_for Metabolomics_and_CVD_Assay.csv')


len(set(namepdf['Sample_ID']))
len(set(namepdf['Specimen_ID']))
sum(~pd.isnull(namepdf['Specimen_ID']))
len(set(namepdf[~pd.isnull(namepdf['Specimen_ID'])]['Specimen_ID']))


### check the repeat ID
from collections import Counter
item_counts = Counter(list(namepdf[~pd.isnull(namepdf['Specimen_ID'])]['Specimen_ID']))
# Finding items with more than one occurrence
repeat_items = [item for item, count in item_counts.items() if count > 1]
print("Repeated items:", repeat_items)



### merge
namelist.columns = ['specimenID', 'Sample_ID']
namepdf['Specimen_ID'] = np.where(pd.isnull(namepdf['Specimen_ID']), 0, namepdf['Specimen_ID'])
namepdf['Specimen_ID'] = namepdf['Specimen_ID'].astype(int).astype(str)
namepdf['Specimen_ID'] = np.where(namepdf['Specimen_ID'] == '0', '', namepdf['Specimen_ID'])

namelist['specimenID'] = namelist['specimenID'].astype(int).astype(str)

name = namepdf.merge(namelist, on=['Sample_ID'], how='left')
print("name.shape:", name.shape)
print(name.head())
# name.to_csv(r'C:\Users\yuqingw1\Workfolder\Data\Metabolomics\namemerge.csv', index=False)


# Open the text file for metabolitic
content = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\Data\Metabolomics\MSC1445_list.txt', encoding='latin1', sep='\t')
print(content.shape)
print(content.columns)
content.head()

# cols = ['1445-' + str(i) + ': Log2(normalized)' for i in range(1, 526)]
cols = ['1445-' + str(i) + '(raw)' for i in range(1, 526)]
caldata = content[cols].T.reset_index()

restcols = list(set(caldata.columns) - set(['index']))
idencol = []
for i in restcols:
    missingrate = sum(caldata[i]==1)/caldata.shape[0]
    if missingrate>0.5:
        idencol.append(i)
print("delete rows number:", len(idencol))

cols2 = ['1445-' + str(i) + ': Log2(normalized)' for i in range(1, 526)]
othcols = ['Compound','Formula','HMP ID', 'KEGG ID', 'LMP ID', 'Mass', 'Retention Time']

metabolomicdata = content.reset_index()
metabolomicdata = metabolomicdata[np.in1d(metabolomicdata['index'], list(set(restcols) - set(idencol)))][cols+cols2+othcols]
metabolomicdata['identify'] = np.where( pd.isnull(metabolomicdata['HMP ID']) & pd.isnull(metabolomicdata['KEGG ID']) & pd.isnull(metabolomicdata['LMP ID']),1,0)
print(metabolomicdata['identify'].value_counts())
metabolomicdata = metabolomicdata[metabolomicdata['identify']==0].reset_index(drop=True)
print(metabolomicdata.shape)

for i in ['HMP ID', 'KEGG ID', 'LMP ID']:
    metabolomicdata[i] = np.where(pd.isnull(metabolomicdata[i]), 'None', metabolomicdata[i])


metabolomicdata.loc[(metabolomicdata['HMP ID']=='HMDB14429') & (metabolomicdata['KEGG ID']=='C19512'), 'KEGG ID' ] ='C06802'
metabolomicdata.loc[(metabolomicdata['LMP ID']=='LMFA08020085') & (metabolomicdata['KEGG ID']=='C16952'), 'KEGG ID' ] ='C06866'
metabolomicdata.loc[(metabolomicdata['HMP ID']=='HMDB31448') & (metabolomicdata['KEGG ID']=='None'), 'KEGG ID' ] ='C01443'
metabolomicdata.loc[(metabolomicdata['HMP ID']=='HMDB40689') & (metabolomicdata['KEGG ID']=='None'), 'KEGG ID' ] ='C17140'


metabolomicdata2 = metabolomicdata[metabolomicdata['KEGG ID']!='None']
print(metabolomicdata2.shape)
newdata2 = pd.DataFrame(metabolomicdata2.groupby(['KEGG ID'], as_index=False)[cols].mean())
newdata2 = pd.merge(metabolomicdata2.groupby(['KEGG ID'], as_index=False)[['HMP ID', 'LMP ID', 'Compound', 'Formula', 'Mass', 'Retention Time']].min(), newdata2, on=['KEGG ID'])
newdata2[cols2] = newdata2[cols].apply(lambda x: np.log(x))
print(newdata2.shape)

metabolomicdata1 = metabolomicdata[metabolomicdata['KEGG ID']=='None']
metabolomicdata3 = metabolomicdata1[metabolomicdata1['HMP ID']!='None']
print(metabolomicdata3.shape)
newdata1 = pd.DataFrame(metabolomicdata3.groupby(['HMP ID'], as_index=False)[cols].mean())
newdata1 = pd.merge(metabolomicdata3.groupby(['HMP ID'], as_index=False)[['KEGG ID', 'LMP ID', 'Compound', 'Formula', 'Mass', 'Retention Time']].min(), newdata1, on=['HMP ID'])
newdata1[cols2] = newdata1[cols].apply(lambda x: np.log(x))
print(newdata1.shape)

metabolomicdata4 = metabolomicdata1[metabolomicdata1['HMP ID']=='None']
print(metabolomicdata4.shape)
newdata3 = pd.DataFrame(metabolomicdata4.groupby(['LMP ID'], as_index=False)[cols].mean())
newdata3 = pd.merge(metabolomicdata4.groupby(['LMP ID'], as_index=False)[['HMP ID', 'KEGG ID', 'Compound', 'Formula', 'Mass', 'Retention Time']].min(), newdata3, on=['LMP ID'])
newdata3[cols2] = newdata3[cols].apply(lambda x: np.log(x))
print(newdata3.shape)


newdata = pd.concat([newdata1, newdata2, newdata3], axis=0)
print(newdata.shape)

print(newdata[['HMP ID', 'KEGG ID', 'LMP ID']])
print(newdata[othcols])

### check the repeat ID
from collections import Counter
item_counts = Counter(list(newdata['HMP ID']))
# Finding items with more than one occurrence
repeat_items = [item for item, count in item_counts.items() if count > 1]
print("Repeated items:", repeat_items, len(repeat_items))
# for i in repeat_items:
#     print(metabolomicdata[metabolomicdata['HMP ID'] == i][['HMP ID', 'KEGG ID', 'LMP ID']])


newdata.to_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\MSC1445_list_filterID.csv', index=False)

newdata = newdata.reset_index(drop=True)
newdata = newdata[cols2].T
print(newdata.shape)
newdata['sampleID'] = list(name['specimenID'])
print(newdata.shape)
newdata = newdata.reset_index(drop=True)


## load microbiome dataset
microbiome = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\microbiome_allvariable_806_transformation.csv')
microcols = [f"UPDASV{i:03}" for i in range(1, 486)]
microbiome['sampleID'] = microbiome['sampleID'].astype(int).astype(str)
print("microbiome data shape:", microbiome.shape)

data = newdata.merge(microbiome[microcols+['sampleID']], on=['sampleID'])
print(data.shape)
data.columns
metabolcol = list(set(data.columns) - set(microcols) - set(['sampleID']))
print("len metaboliccols:", len(metabolcol))

mapdict = {}
for i in metabolcol:
    mapdict[i] = 'Meta_'+str(i)

data.rename(columns=mapdict, inplace=True)
print(data.columns)

data[['sampleID'] + list(mapdict.values()) + microcols].to_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\Metabolomics_filterID.csv', index=False)

# meta = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\Metabolomics_filterID.csv')








