import pandas as pd
import numpy as np

"""
todo:
<1>
matched: multi-key and change name 
nonmatched: non-key examples

replace value

<2> each microbiom value need to convert rate for its row and then do arcsin traansformation
DF2: arc-sin transformation
specimen === sampleid

LR:
systolic
diastolic


check 
(1) nan for these
(2) systolic >> diastolic 
check systolic_1, systolic_2 & diastolic_1 & diastolic_2

(3) all value is 100 
    then transform the value, for the value of bacteria of each sample(mean an individual) transform into a rate with value/sum(value for this row)
(4) each microbiome + age + BMI + gender to do linear regression (response is systolic and diastolic)
age/BMI missing & category 
age

(5) merge all the datasets and duplicates the microbiome (use a new name to represent the same microbiome)

(6) get microbiome coeff and p-value
if new p-value: 0.05/6000 microbiome and select 
return coeff, standard error, pvalue
"""

path1='./Datasets/DFI_16S_rRNA_MiSeq_Sequencing_Result_MMF.16S.404+MMF.16S.408_12+MMF.16S.415_16/'
path2='./Datasets/Microbiome_CACHET/'
path3='./Datasets/DFI_16S_rRNA_MiSeq_Sequencing_Result_MMF.16S.426_428/'



df408_416 = pd.read_csv(path1 + 'MMF.16S.404_MMF.16S.408_MMF.16S.409_MMF.16S.410_MMF.16S.411_MMF.16S.412_MMF.16S.415_MMF.16S.416_JiajunLuo.taxTable_rdp.csv')
dfmicrobiome = pd.read_csv(path2 + 'MB3_MMF.16S.109_MMF.16S.110_MMF.16S.111_SairaTasmin.taxTable_rdp.csv')
df426_428 = pd.read_csv(path3+'MMF.16S.426_MMF.16S.427_MMF.16S.428_JiajunLuo.taxTable_rdp.csv')
print(df408_416.shape)
print(dfmicrobiome.shape)
print(df426_428.shape)


### ensure the large dataset has no repeat key
print(len(set(df408_416['taxonID'].apply(lambda x: x.split(';')[0]))))
print(df408_416.shape)
print(len(set(dfmicrobiome['taxonID'].apply(lambda x: x.split(';')[0]))))
print(dfmicrobiome.shape)
print(len(set(df426_428['taxonID'].apply(lambda x: x.split(';')[0]))))
print(df426_428.shape)

####
def match_taxonID(listofdataframe, listofdfname):
    """
    :param listofdataframe: [df1,df2,df3.....]
                            requirment: each dataframe include at least 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus'
           listofdfname: the name of each dataframe
    :return: matchedresult
    """

    idlist = []
    rstlist = []
    for i in range(len(listofdataframe)):
        listofdataframe[i]['attribute'] = listofdataframe[i][['Kingdom', 'Phylum', 'Class', 'Order', 'Family']].astype(
            str).apply(lambda x: ','.join(x), axis=1)
        idlist.extend(list(listofdataframe[i]['taxonID'].apply(lambda x: x.split(';')[0]+'_%s' %listofdfname[i])))
        rstlist.extend(list(listofdataframe[i]['attribute']))
    multikeylist = []
    restkeylist = []
    for i in set(rstlist):
        matchedkey = []
        for key, value in zip(idlist, rstlist):
            if value == i:
                matchedkey.append(key)
        if len(matchedkey) >= 2:
            multikeylist.append([i, matchedkey])
        else:
            restkeylist.append([i, matchedkey])
    return multikeylist, restkeylist

result1, result2 = match_taxonID([df408_416, dfmicrobiome, df426_428], ['404and416', 'CACHET', '426and428'])
rsdf1 = pd.DataFrame(result1, columns=['attribute', 'taxonID'])
rsdf2 = pd.DataFrame(result2, columns=['attribute', 'taxonID'])
rsdf1.to_csv('result/multikeymapping_fFamily.csv', index=False)
rsdf2.to_csv('result/nonkeymapping_fFamily.csv', index=False)

# rsdf1 = pd.read_csv('result/multikeymapping_for5.csv')
# rsdf2 = pd.read_csv('result/nonkeymapping_for5.csv')



#### replace according to the same key
detail408_416 = pd.read_csv(path1 + 'MMF.16S.404_MMF.16S.408_MMF.16S.409_MMF.16S.410_MMF.16S.411_MMF.16S.412_MMF.16S.415_MMF.16S.416_JiajunLuo.otuTable_rdp.csv')
detailmicrobiome = pd.read_csv(path2 + 'MB2_MMF.16S.109_MMF.16S.110_MMF.16S.111_SairaTasmin.otuTable_rdp.csv')
detail426_428 = pd.read_csv(path3+'MMF.16S.426_MMF.16S.427_MMF.16S.428_JiajunLuo.otuTable_rdp.csv')
print(detail408_416.shape)
print(detailmicrobiome.shape)
print(detail426_428.shape)


### check sampleid if nan
print(sum(pd.isnull(detail408_416['sampleID'])))
print(sum(pd.isnull(detailmicrobiome['sampleID'])))
print(sum(pd.isnull(detail426_428['sampleID'])))

print(detail408_416['sampleID'])
print(detailmicrobiome['sampleID'])
print(detail426_428['sampleID'])

detail408_416['sampleID'] = detail408_416['sampleID'].astype(str)
detail426_428['sampleID'] = detail426_428['sampleID'].apply(lambda x: x.split('_')[-1]).astype(str)
detailmicrobiome['sampleID'] = detailmicrobiome['sampleID'].apply(lambda x: x.split('_')[-1]).astype(str)


detail408_416[detail408_416['sampleID']=='102853']['ASV46880;seqs=25;samples=1']
sum(detail426_428['sampleID']=='102853')
sum(detailmicrobiome['sampleID']=='102853')


def replace_keyid(listofdataframe, listofdfname, multimapping):

    # unique_numbers = np.random.choice(range(10000, 100000), size=multimapping.shape[0]+nonmapping.shape[0], replace=False)
    sequence = [f"{i:03}" for i in range(1, multimapping.shape[0]+1)]
    multimapping['newtaxonID'] = list(sequence)[:multimapping.shape[0]]
    multimapping['newtaxonID'] = 'FAMILY' + multimapping['newtaxonID'].astype(str)
    mappingdf = multimapping[['newtaxonID', 'taxonID', 'attribute']]
    mappingdf = mappingdf.reset_index(drop=True)
    finaldf = pd.DataFrame()
    for i in range(len(listofdataframe)):
        result = {}
        all_columns = list(set(listofdataframe[i].columns) - set(['sampleID']))
        for col in all_columns:
            tmpcol = col.split(';')[0] + '_' + listofdfname[i]
            for j in range(len(mappingdf)):
                if tmpcol in mappingdf.loc[j, 'taxonID']:
                    if mappingdf.loc[j, 'newtaxonID'] in result:
                        result[mappingdf.loc[j, 'newtaxonID']] += np.array(listofdataframe[i][col])
                    else:
                        result[mappingdf.loc[j, 'newtaxonID']] = np.array(listofdataframe[i][col])

        eachdata = pd.DataFrame(result)
        eachdata['sampleID'] = listofdataframe[i]['sampleID']
        finaldf = pd.concat([finaldf, eachdata], axis=0, join='outer', ignore_index=True)
    print("finaldf shape:", finaldf.shape)
    transcols = list(set(finaldf.columns) - set(['sampleID']))
    print("length cols:", len(transcols))
    finaltrans = finaldf.copy()
    finaltrans = finaltrans.fillna(0)
    finaltrans = finaltrans[transcols].div(finaltrans[transcols].sum(axis=1), axis=0)
    finaldf_trans = pd.concat([finaldf['sampleID'], finaltrans], axis=1)
    # from collections import Counter
    # item_counts = Counter(list(finaldf_trans['sampleID']))
    # repeat_items = [item for item, count in item_counts.items() if count > 1]
    # print("Repeated items:", repeat_items)
    finaldata = finaldf_trans.groupby(['sampleID'], as_index=False)[transcols].mean()
    print("finaldata shape:", finaldata.shape)
    finaldata = finaldata.reset_index(drop=True)
    finaldata = pd.concat([finaldata['sampleID'], np.arcsin(np.sqrt(finaldata[transcols]))], axis=1)
    finaldata[transcols] = finaldata[transcols].apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=0), axis=1)
    return finaldf, finaldata, mappingdf

fdata, fdata_trans, mappingdf = replace_keyid([detail408_416, detailmicrobiome, detail426_428], ['404and416', 'CACHET', '426and428'], rsdf1 )

print(fdata.shape)
print(fdata.columns)
print(fdata_trans.shape)
mappingdf
fdata_trans.head()

# fdata.to_csv('Datasets/microbiome_combine.csv', index=False)
fdata_trans.to_csv('Datasets/microbiome_combine_trans_forfamily.csv', index=False)
# mappingdf.to_csv('result/taxonIDmappingdict.csv', index=False)

# fdata_trans[fdata_trans['sampleID']=='102853']['NEWASV414']

### check the repeat ID
from collections import Counter
item_counts = Counter(list(fdata_trans['sampleID']))
# Finding items with more than one occurrence
repeat_items = [item for item, count in item_counts.items() if count > 1]
print("Repeated items:", repeat_items)



mappingdf = pd.read_excel('result/taxonIDmappingdict_fGenus.xlsx')
df = pd.read_csv('./Datasets/processing_4_microbiome/microbiome_combine_790row_480col.csv')

def Directlyreplace_keyid(mapping, df):
    mapping['Family_attr'] = mapping['attribute'].apply(lambda x: x.split(',')[:-1]).apply(lambda x: ','.join(x))
    mapping['taxonID'] = mapping['taxonID'].apply(lambda x: x.strip('[]'))
    newmapping = mapping.groupby(['Family_attr'], as_index=False)[['newtaxonID', 'taxonID']].agg(lambda x: ','.join(x)).reset_index()

    sequence = [f"{i:03}" for i in range(1, newmapping.shape[0]+1)]
    newmapping['FamilytaxonID'] = list(sequence)[:newmapping.shape[0]]
    newmapping['FamilytaxonID'] = 'FAMILY' + newmapping['FamilytaxonID'].astype(str)
    newmappingdf = newmapping[['FamilytaxonID', 'newtaxonID', 'taxonID', 'Family_attr']]
    newmappingdf = newmappingdf.reset_index(drop=True)
    print("lenght:", len(list(mapping['newtaxonID'])))

    result = {}
    for col in list(mapping['newtaxonID']):
        df[col] = df[col].astype(float)
        for j in range(len(newmappingdf)):
            if col in newmappingdf.loc[j, 'newtaxonID']:
                if newmappingdf.loc[j, 'FamilytaxonID'] in result:
                    result[newmappingdf.loc[j, 'FamilytaxonID']] += np.array(df[col])
                else:
                    result[newmappingdf.loc[j, 'FamilytaxonID']] = np.array(df[col])

    finaldf = pd.DataFrame(result)
    finaldf['sampleID'] = df['sampleID']
    print("finaldf shape:", finaldf.shape)
    transcols = list(set(finaldf.columns) - set(['sampleID']))
    print("length cols:", len(transcols))
    finaltrans = finaldf.fillna(0)
    finaltrans = finaltrans[transcols].div(finaltrans[transcols].sum(axis=1), axis=0)
    finaldf_trans = pd.concat([finaldf['sampleID'], finaltrans], axis=1)
    finaldata = finaldf_trans.groupby(['sampleID'], as_index=False)[transcols].mean()
    print("finaldata shape:", finaldata.shape)
    finaldata = finaldata.reset_index(drop=True)
    finaldata = pd.concat([finaldata['sampleID'], np.arcsin(np.sqrt(finaldata[transcols]))], axis=1)
    finaldata[transcols] = finaldata[transcols].apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=0), axis=1)
    return finaldf, finaldata, newmappingdf

fdata, fdata_trans, mappingdf = Directlyreplace_keyid(mappingdf, df)

print(fdata.shape)
print(fdata.columns)
print(fdata_trans.shape)

mappingdf.to_excel('result/taxonIDmappingdict_fFamily.xlsx', index=False)
fdata.to_csv('Datasets/processing_4_micrFamily/micro_combineFamily_790row_191col.csv', index=False)
fdata_trans.to_csv('Datasets/processing_4_micrFamily/microFamily_trans.csv', index=False)



