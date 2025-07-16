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

path=r'C:\Users\yuqingw1\Workfolder'


df448_449 = pd.read_csv(path + '\Data\DFI 16S rRNA MiSeq Sequencing Result (MMF.16S.448 + MMF.16S.449)\MMF.16S.448_MMF.16S.449_JiajunLuo.taxTable_rdp.csv')
print(df448_449.shape)
df404_415 = pd.read_csv(path + '\Data\DFI 16S rRNA MiSeq Sequencing Result (MMF.16S.404 + MMF.16S.408-12 + MMF.16S.415-16)\MMF.16S.404_MMF.16S.408_MMF.16S.409_MMF.16S.410_MMF.16S.411_MMF.16S.412_MMF.16S.415_MMF.16S.416_JiajunLuo.taxTable_rdp.csv')
print(df404_415.shape)
df426_428 = pd.read_csv(path + '\Data\DFI 16S rRNA MiSeq Sequencing Result (MMF.16S.426-428)\MMF.16S.426_MMF.16S.427_MMF.16S.428_JiajunLuo.taxTable_rdp.csv')
print(df426_428.shape)
dfcachet = pd.read_csv(path + '\Data\Microbiome_CACHET\MMF.16S.109_MMF.16S.110_MMF.16S.111_SairaTasmin.taxTable_rdp.csv')
print(dfcachet.shape)



### ensure the large dataset has no repeat key
print(len(set(df448_449['taxonID'].apply(lambda x: x.split(';')[0]))))
print(df448_449.shape)
print(len(set(df404_415['taxonID'].apply(lambda x: x.split(';')[0]))))
print(df404_415.shape)
print(len(set(df426_428['taxonID'].apply(lambda x: x.split(';')[0]))))
print(df426_428.shape)
print(len(set(dfcachet['taxonID'].apply(lambda x: x.split(';')[0]))))
print(dfcachet.shape)

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
        listofdataframe[i][['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']] = listofdataframe[i][['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']].apply(lambda x: x.str.replace(' ', '_'))
        listofdataframe[i]['attribute'] = listofdataframe[i][['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']].astype(
            str).apply(lambda x: ','.join(x), axis=1)
        idlist.extend(list(listofdataframe[i]['taxonID'].apply(lambda x: x.split(';')[0]+'_%s' %listofdfname[i])))
        rstlist.extend(list(listofdataframe[i]['attribute']))
    # print("rstlist:", rstlist[0])
    # print("idlist:", idlist[0])
    multikeylist = []
    for i in set(rstlist):
        matchedkey = ''
        for key, value in zip(idlist, rstlist):
            if value == i:
                matchedkey += str(key)+' '
        multikeylist.append([i, matchedkey])
    return multikeylist

result= match_taxonID([df448_449, df404_415, df426_428, dfcachet], ['df448_449', 'df404_415','df426_428','dfcachet'])


rsdf = pd.DataFrame(result, columns=['attribute', 'taxonID'])
# rsdf.to_csv(path + '/result/multikeymapping.csv', index=False)

# name = pd.read_excel(path + '/ProcessingData/processing_4_microbiome/taxonIDmappingdict.xlsx')
# name = name[['attribute', 'newtaxonID']]
# rsdf = rsdf.merge(name, on=['attribute'])
rsdf['newtaxonID'] = [f"UPDASV{i:03}" for i in range(1, 486)]
rsdf.shape

#### replace according to the same key
df448_449_m = pd.read_csv(path + '\Data\DFI 16S rRNA MiSeq Sequencing Result (MMF.16S.448 + MMF.16S.449)\MMF.16S.448_MMF.16S.449_JiajunLuo.otuTable_rdp.csv')
print(df448_449_m.shape)
df404_415_m = pd.read_csv(path + '\Data\DFI 16S rRNA MiSeq Sequencing Result (MMF.16S.404 + MMF.16S.408-12 + MMF.16S.415-16)\MMF.16S.404_MMF.16S.408_MMF.16S.409_MMF.16S.410_MMF.16S.411_MMF.16S.412_MMF.16S.415_MMF.16S.416_JiajunLuo.otuTable_rdp.csv')
print(df404_415_m.shape)
df426_428_m = pd.read_csv(path + '\Data\DFI 16S rRNA MiSeq Sequencing Result (MMF.16S.426-428)\MMF.16S.426_MMF.16S.427_MMF.16S.428_JiajunLuo.otuTable_rdp.csv')
print(df426_428_m.shape)
dfcachet_m = pd.read_csv(path + '\Data\Microbiome_CACHET\MMF.16S.109_MMF.16S.110_MMF.16S.111_SairaTasmin.otuTable_rdp.csv')
print(dfcachet_m.shape)



### check sampleid if nan
print(sum(pd.isnull(df448_449_m['sampleID'])))
print(sum(pd.isnull(df404_415_m['sampleID'])))
print(sum(pd.isnull(df426_428_m['sampleID'])))
print(sum(pd.isnull(dfcachet_m['sampleID'])))

print(df448_449_m['sampleID'])
print(df404_415_m['sampleID'])
print(df426_428_m['sampleID'])
print(dfcachet_m['sampleID'])

df404_415_m['sampleID'] = df404_415_m['sampleID'].astype(str)
df426_428_m['sampleID'] = df426_428_m['sampleID'].apply(lambda x: x.split('_')[-1]).astype(str)
dfcachet_m['sampleID'] = dfcachet_m['sampleID'].apply(lambda x: x.split('_')[-1]).astype(str)
df448_449_m['sampleID'] = df448_449_m['sampleID'].apply(lambda x: x.split('_')[-1]).astype(str)


# df404_415_m[df404_415_m['sampleID']=='102853']['ASV46880;seqs=25;samples=1']
# sum(df426_428_m['sampleID']=='102853')
# sum(df426_428_m['sampleID']=='102853')

def replace_keyid(listofdataframe, listofdfname, multimapping):

    # multimapping = multimapping.sort_values(by=['attribute', 'taxonID']).reset_index(drop=True)
    print("multimapping shape:", multimapping.shape[0])
    # unique_numbers = np.random.choice(range(10000, 100000), size=multimapping.shape[0]+nonmapping.shape[0], replace=False)
    # sequence = [f"{i:03}" for i in range(1, multimapping.shape[0]+1)]
    # multimapping['newtaxonID'] = list(sequence)
    # multimapping['newtaxonID'] = 'UPDASV' + multimapping['newtaxonID'].astype(str)

    mappingdf = multimapping
    finaldf=pd.DataFrame()
    for i in range(len(listofdataframe)):
        result = {}
        all_columns = list(set(listofdataframe[i].columns) - set(['sampleID']))
        print("all_columns", all_columns)
        # print("all_columns", all_columns[:5])
        for col in all_columns:
            # print("col", col)
            tmpcol = col.split(';')[0] + '_' + listofdfname[i]
            # print("tmpcol", tmpcol)
            for j in range(len(mappingdf)):
                # print("j", j)
                # print("check1:",  mappingdf.loc[j, 'taxonID'].split(' ')[:-1])
                if tmpcol in mappingdf.loc[j, 'taxonID'].split(' ')[:-1]:
                    # print("result:", result)
                    if mappingdf.loc[j, 'newtaxonID'] in result:
                        # print("array:", np.array(listofdataframe[i][col]))
                        result[mappingdf.loc[j, 'newtaxonID']] += np.array(listofdataframe[i][col])
                    else:
                        result[mappingdf.loc[j, 'newtaxonID']] = np.array(listofdataframe[i][col])
        eachdata = pd.DataFrame(result)
        eachdata['sampleID'] = listofdataframe[i]['sampleID']
        finaldf = pd.concat([finaldf, eachdata], axis=0, join='outer', ignore_index=True)
    print("finaldf shape:", finaldf.shape)
    finaldf = finaldf.drop_duplicates().reset_index(drop=True)
    print("finaldf first drop duplicates shape:", finaldf.shape)
    transcols = list(set(finaldf.columns) - set(['sampleID']))
    print("length cols:", len(transcols))
    finaltrans = finaldf.copy()
    finaltrans = finaltrans.fillna(0)
    finaltrans = finaltrans[transcols].div(finaltrans[transcols].sum(axis=1), axis=0)
    finaldf_trans = pd.concat([finaldf['sampleID'], finaltrans], axis=1)
    # print("check2finaltrans:", finaldf_trans.head(10))
    # from collections import Counter
    # item_counts = Counter(list(finaldf_trans['sampleID']))
    # repeat_items = [item for item, count in item_counts.items() if count > 1]
    # print("Repeated items:", repeat_items)
    finaldata = finaldf_trans.groupby(['sampleID'], as_index=False)[transcols].mean()
    # print("finaldata shape:", finaldata.head(10))
    finaldata = finaldata.reset_index(drop=True)
    finaldata[transcols] = np.arcsin(np.sqrt(finaldata[transcols]))
    # print("check2finaltrans:", finaldata[transcols].describe())
    finaldata[transcols] = finaldata[transcols].apply(lambda x: (x - np.mean(x)) / np.std(x, ddof=0), axis=1)
    # print("check3finaltrans:", finaldata.head(10))
    mappingdf['taxonID'] = mappingdf['taxonID'].apply(lambda x: x.split(' ')[:-1])
    return finaldf, finaldata, mappingdf

fdata, fdata_trans, mappingdf = replace_keyid([df448_449_m, df404_415_m, df426_428_m, dfcachet_m], ['df448_449', 'df404_415','df426_428','dfcachet'], rsdf )

# fdata, fdata_trans, mappingdf = replace_keyid([df404_415_m[['sampleID','ASV11981;seqs=543;samples=37','ASV85372;seqs=2;samples=1','ASV47538;seqs=24;samples=2']]], ['df404_415'], rsdf.sort_values(['newtaxonID']).iloc[:11] )

print(fdata.shape)
print(fdata.columns)
print(fdata_trans.shape)
fdata_trans.head()

fdata.to_csv(path + '/ProcessingData/processing_4_microbiome_new/microbiome_combine_879row_485col.csv', index=False)
fdata.groupby(['sampleID'], as_index=False).mean().to_csv(path + '/ProcessingData/processing_4_microbiome_new/microbiome_combine_806row_485col.csv', index=False)
fdata_trans.to_csv(path + '/ProcessingData/processing_4_microbiome_new/microbiome_combine_trans_806_485.csv', index=False)
mappingdf.to_excel(path + '/ProcessingData/processing_4_microbiome_new/taxonIDmappingdict.xlsx', index=False)



# ### check the repeat ID
# from collections import Counter
# item_counts = Counter(list(fdata['sampleID']))
# # Finding items with more than one occurrence
# repeat_items = [item for item, count in item_counts.items() if count > 1]
# print("Repeated items:", repeat_items)
#
# fdata[np.in1d(fdata['sampleID'], repeat_items)].to_csv('repeated_items.csv', index=False)
#
# ####
# df_repeat = fdata_trans[(fdata_trans['sampleID']=='104301') | (fdata_trans['sampleID']=='104277') | (fdata_trans['sampleID']=='104390') | (fdata_trans['sampleID']=='106701')]
# # df_repeat.to_csv('result/incorrectmicrobiome.csv', index=False)
#
# mappingdf = pd.read_csv('result/taxonIDmappingdict.csv')
# # mappingdf[mappingdf['newtaxonID']=='NEWASV146']['taxonID'].to_list()
# # mappingdf[mappingdf['newtaxonID']=='NEWASV147']
#
# mappingdf.to_excel('result/taxonIDmappingdict.xlsx', index=False)
#
#
#
# import pandas as pd
# import numpy as np
# trans = pd.read_csv('Datasets/microbiome_trans_compass_pm25.csv')
# selectcols=['NEWASV040', 'NEWASV041', 'NEWASV047', 'NEWASV051', 'NEWASV061',
#  'NEWASV062', 'NEWASV069', 'NEWASV074', 'NEWASV102', 'NEWASV113', 'NEWASV128',
#  'NEWASV131', 'NEWASV132', 'NEWASV179', 'NEWASV198', 'NEWASV211', 'NEWASV217',
#  'NEWASV218', 'NEWASV219', 'NEWASV228', 'NEWASV240', 'NEWASV241', 'NEWASV243',
#  'NEWASV244', 'NEWASV255', 'NEWASV264', 'NEWASV274', 'NEWASV275', 'NEWASV279',
#  'NEWASV280', 'NEWASV286', 'NEWASV289', 'NEWASV292', 'NEWASV312', 'NEWASV316', 'NEWASV321']
# trans[selectcols] = trans[selectcols].apply(lambda x: -x)
# for i in selectcols:
#     trans.rename(columns={i: i+'_reverse'}, inplace=True)
# micro = [f"NEWASV{i:03}" for i in range(1, 481)]
# micro = [item for item in micro if item not in selectcols]
# trans.drop(micro, axis=1, inplace=True)
# trans.shape
# trans.columns
# trans.to_csv('Datasets/microbiome_trans_compass_pm25_reverse.csv')