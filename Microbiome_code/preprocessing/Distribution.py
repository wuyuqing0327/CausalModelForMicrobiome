#### comprehensive data
import pandas as pd
import numpy as np
path=r'C:\Users\yuqingw1\Workfolder\ProcessingData'

micro_org=pd.read_csv(path+'\microbiome_allvariable_806_original.csv')

micro_org = micro_org[ ~pd.isnull(micro_org['Resistin'] )]
# micro_org = micro_org[ ~pd.isnull(micro_org['Insulin'] )]
# micro_org['Insulin'].describe()

print(micro_org.shape)

print(sum(pd.isnull(micro_org['sampleID'])))

print(sorted(micro_org.columns)[-100:])

'CKMB', 'Cpeptide', 'Ghrelin', 'HBP', 'IL6', 'IL6_gp', 'Insulin', 'Leptin',
'age', 'bmi', 'date', 'date2', 'diabetes2', 'diastolic', 'education', 'gender', 'household_income', 'race', 'sampleID', 'smoking_status', 'systolic'


micro_org['age'].describe()

micro_org['gender'] = np.where(pd.isnull(micro_org['gender']), 'missing', micro_org['gender'])
micro_org['gender'].value_counts()

micro_org['smoking_status'] = np.where(pd.isnull(micro_org['smoking_status']), 'missing', micro_org['smoking_status'])
micro_org['smoking_status'].value_counts()


cutPoint = [18.5, 25]
micro_org['bmi_binned2'] = pd.cut(micro_org[~pd.isnull(micro_org['bmi'])]['bmi'],bins=[-np.inf] + cutPoint + [np.inf]).astype(str)
micro_org['bmi_binned2'] = np.where(pd.isnull(micro_org['bmi_binned2']), 'missing', micro_org['bmi_binned2'])
print(micro_org['bmi_binned2'].value_counts())


micro_org['diabetes2'] = np.where(pd.isnull(micro_org['diabetes2']), 'missing', micro_org['diabetes2'])
micro_org['diabetes2'].value_counts()


micro_org['race'] = np.where(pd.isnull(micro_org['race']), 'missing', micro_org['race'])
micro_org['race'].value_counts()

micro_org['Resistin'].describe()
micro_org['Ghrelin'].describe()
micro_org['Leptin'].describe()

shannon = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\processing_4_microbiome_new\microbiome_combine_806row_485col_shannon.csv')
shannon.shape
shannon.columns
shannon = shannon.merge(micro_org['sampleID'])
print(shannon.shape)
shannon['Shannon_Diversity'].describe()



micro_cols = [f"UPDASV{i:03}" for i in range(1, 486)]
data = micro_org
print(data.shape)

# data.head()
taxonIDname = pd.read_excel(r'C:\Users\yuqingw1\Workfolder\ProcessingData\processing_4_microbiome_new\taxonIDmappingdict.xlsx')

taxonIDname[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus']] = taxonIDname['attribute'].str.split(',', expand=True)
#
taxonIDname['category'] = np.where(taxonIDname['Family']=='Lachnospiraceae', 'Lachnospiraceae (Family)',
    np.where(taxonIDname['Family']=='Oschillospiraceae', 'Oschillospiraceae (Family)',
    np.where(taxonIDname['Class'] == 'Other Clostridia', 'Other Clostridia (Class)',
    np.where(taxonIDname['Phylum'] == 'Bacteroidetes', 'Bacteroidetes (Phylum)',
    np.where(taxonIDname['Phylum'] == 'Actinobacteria', 'Actinobacteria (Phylum)',
    np.where(taxonIDname['Phylum'] == 'Proteobacteria', 'Proteobacteria (Phylum)',
    np.where(taxonIDname['Family'] == 'Erysipelotrichaceae', 'Erysipelotrichaceae (Family)',
    np.where(taxonIDname['Genus'] == 'Parvimonas', 'Parvimonas  (Genus)',
    np.where(taxonIDname['Genus'] == 'Tissierella', 'Tissierella (Genus)',
    np.where(taxonIDname['Genus'] == 'Pectinatus', 'Pectinatus (Genus)',
    'Other Bacteria (Kingdom)'))))))))))


print(taxonIDname['category'].value_counts())

taxonIDname[np.in1d(taxonIDname['newtaxonID'], ['UPDASV480', 'UPDASV053', 'UPDASV051', 'UPDASV186', 'UPDASV236'])][['Class','Genus', 'category']]

data[micro_cols] = data[micro_cols].fillna(0)
resultdf = data[micro_cols].sum(axis=1)
rate_df = data[micro_cols].div(resultdf, axis=0)
rate_df['sampleID'] = data['sampleID']
rate_df['sampleID'] = rate_df['sampleID'].astype(int).astype(str)

statdf = rate_df[list(set(rate_df.columns) - set(['sampleID']))].mean()
statdf = statdf.reset_index()
statdf.columns=['sampleID', 'mean']

statdf = statdf.merge(taxonIDname[['newtaxonID', 'Genus']], left_on=['sampleID'], right_on=['newtaxonID'])
statdf
statdf.sort_values(by=['mean'], ascending=False)

rate_df[['UPDASV480', 'UPDASV053', 'UPDASV051', 'UPDASV186', 'UPDASV236']].mean()

rate_df[['Streptococcus', 'Prevotella', 'Haemophilus', 'Veillonella', 'Fusobacterium']]
