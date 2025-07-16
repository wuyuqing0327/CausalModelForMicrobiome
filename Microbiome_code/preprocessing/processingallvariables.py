#### comprehensive data
import pandas as pd
import numpy as np
path=r'C:\Users\yuqingw1\Workfolder'
path2 = path + '\ProcessingData\processing_4_microbiome_new'

micro_org=pd.read_csv(path2+'\microbiome_combine_806row_485col.csv')
micro_trans=pd.read_csv(path2+'\microbiome_combine_trans_806_485.csv')
# shannon = pd.read_csv(path2+'\microbiome_combine_806row_485col_shannon.csv')
# bray_curtis = pd.read_csv(path2+'\microbiome_combine_806row_bray_curtis.csv')
# bray_curtis.rename(columns={'Unnamed: 0':'sampleID'}, inplace=True)

### concatenate all variables
micro_org = micro_org.groupby(['sampleID'], as_index=False).mean()
basicdetails = pd.read_csv(path+'/Data/compass_deidentified.csv')
print("basecdetail:", basicdetails.shape)
basicdetails = basicdetails[~pd.isnull(basicdetails['specimen_id'])].reset_index(drop=True)
print("basecdetail:", basicdetails.shape)
pm25 = pd.read_csv(path+ '/Data/compass_pm25_Mvalue.csv')
print("pm25:", pm25.shape)
pm25 = pm25[~pd.isnull(pm25['specimen_id'])].reset_index(drop=True)
print("pm25:", pm25.shape)
compass_cardio =pd.read_csv(path +'/Data/compass_cardiometabolic_marker.csv')
print("compass_cardio:", compass_cardio.shape)
compass_cardio = compass_cardio[~pd.isnull(compass_cardio['specimen_id'])].reset_index(drop=True)
print("compass_cardio:", compass_cardio.shape)
insulin = pd.read_csv(path + '/Data/cardiometabolic_insulin.csv')
print("insulin:", insulin.shape)
insulin = insulin[~pd.isnull(insulin['specimen_id'])].reset_index(drop=True)
print("insulin:", insulin.shape)


micro_org['sampleID'] = micro_org['sampleID'].astype(str)
micro_trans['sampleID'] = micro_trans['sampleID'].astype(str)
# shannon['sampleID'] = shannon['sampleID'].astype(str)
# bray_curtis['sampleID'] = bray_curtis['sampleID'].astype(str)

pm25['specimen_id'] = pm25['specimen_id'].astype(int).astype(str)
compass_cardio['specimen_id'] = compass_cardio['specimen_id'].astype(int).astype(str)
insulin['specimen_id'] = insulin['specimen_id'].astype(int).astype(str)
basicdetails['specimen_id'] = basicdetails['specimen_id'].astype(int).astype(str)

### processing detail
basicdetails.loc[basicdetails['specimen_id']=='106503', 'age'] = (53.36 + 56.77)/2
print(basicdetails.shape)

def adjustbloodpressure(data, sampleID:str):
    """
    Regualations for adjusting systolic and diastolic
    Regarding systolic and diastolic
    1) LOOK INTO systolic_1, systolic_2 & diastolic_1, diastolic_2
    do turn over its value if it is needed
    2) if all the value in the same value then delete
    3)

    """
    ### exchange systolic and diastolic
    equalcon = (data['systolic']==data['diastolic']) & (data['systolic']==data['diastolic_2']) & (data['systolic']==data['diastolic_1']) & (data['systolic']== data['systolic_1']) & (data['systolic']==data['systolic_2'])
    print(sum(equalcon))
    mask1 = data['systolic'] <= data['diastolic']
    print( " incorrrect SBD AND DBD 1:", sum(~equalcon & mask1))

    mask2 = ( data['systolic_1'] <= data['diastolic_1']) | ( data['systolic_2'] <= data['diastolic_2'] )
    print( " incorrrect SBD AND DBD 2:", sum(~mask1 & mask2 ))

    mask3 = ( (data['systolic'] - data['diastolic']) <=15 )
    print( " incorrrect SBD AND DBD 2:", sum( ~mask1 & ~mask2 & mask3 ))

    data.loc[data[sampleID]=='100283', 'diastolic'] = 99.0
    data.loc[data[sampleID]=='100995', 'diastolic'] = 78.0
    data.loc[data[sampleID]=='101252', 'diastolic'] = 63.0
    data.loc[data[sampleID]=='102389', 'diastolic'] = 76.0
    data.loc[data[sampleID]=='102746', 'diastolic'] = 90.0
    data.loc[data[sampleID]=='104672', 'diastolic'] = 78.0
    data.loc[data[sampleID]=='107830', 'diastolic'] = 108.0
    data.loc[data[sampleID]=='108747', 'diastolic'] = 91.0

    data.loc[data[sampleID]=='101748', 'systolic'] = 162.0
    data.loc[data[sampleID]=='101779', 'systolic'] = 103.0
    data.loc[data[sampleID]=='101549', 'systolic'] = 165.0
    data.loc[data[sampleID]=='103150', 'systolic'] = 120.0


    data.loc[data[sampleID]=='101633', 'systolic'] = (153.0 + 162.0)/2
    data.loc[data[sampleID]=='101633', 'diastolic'] = (100.0 + 90.0 + 94.0 + 100)/4

    data.loc[data[sampleID]=='100738', 'systolic'] = (194.0 + 156.0 + 196.0 + 155.0)/4
    data.loc[data[sampleID]=='100738', 'diastolic'] = (146.0 + 105.0 + 147.0 + 100.0)/4


    data.loc[data[sampleID]=='102500', 'systolic'] = (156.0 + 157.0)/2
    data.loc[data[sampleID]=='102500', 'diastolic'] = (103.0 + 93.0)/2


    data.loc[data[sampleID]=='106984', 'systolic'] = (139.0 + 136.0)/2
    data.loc[data[sampleID]=='106984', 'diastolic'] = (78.0 + 80.0)/2

    data.loc[data[sampleID]=='107730', 'systolic'] = (156.0 + 149.0 + 152.0)/3
    data.loc[data[sampleID]=='107730', 'diastolic'] = 98.0

    data.loc[data[sampleID]=='107956', 'systolic'] = (140.0 + 134.0)/2
    data.loc[data[sampleID]=='107956', 'diastolic'] = (86.0 + 85.0)/2

    data.loc[data[sampleID]=='108578', 'systolic'] = (113.0 + 108.0)/2
    data.loc[data[sampleID]=='108578', 'diastolic'] = (66.0 + 69.0)/2

    data.loc[data[sampleID]=='108601', 'systolic'] = 107.0
    data.loc[data[sampleID]=='108601', 'diastolic'] = (81.0 + 80.0)/2

    data.loc[data[sampleID]=='108894', 'systolic'] = (142.0 + 129.0)/2
    data.loc[data[sampleID]=='108894', 'diastolic'] = (76.0 + 75.0)/2

    data.loc[data[sampleID]=='103535', 'systolic'] = 100.0
    data.loc[data[sampleID]=='103535', 'diastolic'] = 80.0

    data.loc[data[sampleID]=='100346', 'systolic'] = (182.0 + 177.0)/2
    data.loc[data[sampleID]=='100346', 'diastolic'] = 123.0

    data.loc[data[sampleID]=='100713', 'systolic'] = (176.0 + 183.0)/2
    data.loc[data[sampleID]=='100713', 'diastolic'] = 110.0

    data.loc[data[sampleID]=='101169', 'systolic'] = 165.0
    data.loc[data[sampleID]=='101169', 'diastolic'] = 109.0


    data.loc[data[sampleID]=='101127', 'systolic'] = (141.0 + 143.0)/2
    data.loc[data[sampleID]=='101127', 'diastolic'] = (118.0 + 112.0)/2


    data.loc[data[sampleID]=='101495', 'systolic'] = (133.0 + 123.0)/2
    data.loc[data[sampleID]=='101495', 'diastolic'] = (114.0 + 81.0)/2

    data.loc[data[sampleID]=='101785', 'systolic'] = 134.0
    data.loc[data[sampleID]=='101785', 'diastolic'] = 98.0

    data.loc[data[sampleID]=='102259', 'systolic'] = 140.0
    data.loc[data[sampleID]=='102259', 'diastolic'] = 90.0

    data.loc[data[sampleID]=='102523', 'systolic'] = (116.0 + 117.0)/2
    data.loc[data[sampleID]=='102523', 'diastolic'] = 79.0

    data.loc[data[sampleID]=='102704', 'systolic'] = (157.0 + 160.0)/2
    data.loc[data[sampleID]=='102704', 'diastolic'] = 85.0

    data.loc[data[sampleID]=='102738', 'systolic'] = 138.0
    data.loc[data[sampleID]=='102738', 'diastolic'] = 91.0

    data.loc[data[sampleID]=='102606', 'systolic'] = (138.0 + 145)/2
    data.loc[data[sampleID]=='102606', 'diastolic'] = 105.0

    data.loc[data[sampleID]=='103262', 'systolic'] = (138.0 + 134)/2
    data.loc[data[sampleID]=='103262', 'diastolic'] = 84.0

    data.loc[data[sampleID]=='103112', 'systolic'] = (151.0 + 129)/2
    data.loc[data[sampleID]=='103112', 'diastolic'] = 82.0

    data.loc[data[sampleID]=='103345', 'systolic'] = (132.0 + 131)/2
    data.loc[data[sampleID]=='103345', 'diastolic'] = (80.0 + 99)/2

    data.loc[data[sampleID]=='103554', 'systolic'] = 178
    data.loc[data[sampleID]=='103554', 'diastolic'] = 106

    data.loc[data[sampleID]=='102484', 'systolic'] = (105.0 + 114)/2
    data.loc[data[sampleID]=='102484', 'diastolic'] = (69.0 + 73)/2

    data.loc[data[sampleID]=='106969', 'systolic'] = (125.0 + 134 + 122 + 126)/4
    data.loc[data[sampleID]=='106969', 'diastolic'] = (79.0 + 83 + 78 + 80)/4

    data.loc[data[sampleID]=='108143', 'systolic'] = (132.0 + 133)/2
    data.loc[data[sampleID]=='108143', 'diastolic'] = 96.0

    data.loc[data[sampleID]=='108145', 'systolic'] = (127.0 + 135)/2
    data.loc[data[sampleID]=='108145', 'diastolic'] = 75.0

    data.loc[data[sampleID]=='106906', 'systolic'] = (126.0 + 122 )/2
    data.loc[data[sampleID]=='106906', 'diastolic'] = (78 + 80)/2

    data.loc[data[sampleID]=='102865', 'systolic'] = 159
    data.loc[data[sampleID]=='102865', 'diastolic'] = (125 + 122)/2

    data.loc[data[sampleID]=='103233', 'systolic'] = 103
    data.loc[data[sampleID]=='103233', 'diastolic'] = 47

    cols = [sampleID, 'age', 'bmi', 'gender', 'systolic', 'diastolic', 'smoking_status', 'diabetes2']
    data = data.drop_duplicates(cols)
    print("final data:", data.shape)
    data = data.reset_index(drop=True)
    return data

basicdetails2 = adjustbloodpressure(basicdetails, 'specimen_id')


covariates = ['age', 'bmi', 'gender','systolic', 'diastolic','smoking_status', 'diabetes2', 'race', 'household_income', 'education']
basicdetails2 = basicdetails2.drop_duplicates(['specimen_id'] +covariates, keep='first')
print(basicdetails2.shape)
basicdetails3 = basicdetails2.drop_duplicates(['specimen_id'], keep='first')
print(basicdetails3.shape)



#### merge
data_org = micro_org.merge(pm25, left_on=['sampleID'], right_on=['specimen_id'], how='left').drop(['specimen_id'],axis=1)\
          .merge(basicdetails3[['specimen_id']+covariates], left_on=['sampleID'], right_on=['specimen_id'], how='left').drop(['specimen_id'],axis=1)\
         .merge(compass_cardio, left_on=['sampleID'], right_on=['specimen_id'], how='left').drop(['specimen_id'],axis=1)\
         .merge(insulin, left_on=['sampleID'], right_on=['specimen_id'], how='left').drop(['specimen_id'], axis=1)
print(data_org.shape)
print(data_org.columns)
data_org['HBP'] = np.where(pd.isnull(data_org['systolic']) | pd.isnull(data_org['diastolic']), np.nan,
                           np.where( (data_org['systolic'] >= 140) | (data_org['diastolic'] >= 90), 1, 0))


data_trans = micro_trans.merge(pm25, left_on=['sampleID'], right_on=['specimen_id'], how='left').drop(['specimen_id'],axis=1)\
          .merge(basicdetails3[['specimen_id']+covariates], left_on=['sampleID'], right_on=['specimen_id'], how='left').drop(['specimen_id'],axis=1)\
         .merge(compass_cardio, left_on=['sampleID'], right_on=['specimen_id'], how='left').drop(['specimen_id'],axis=1)\
         .merge(insulin, left_on=['sampleID'], right_on=['specimen_id'], how='left').drop(['specimen_id'], axis=1)
print(data_trans.shape)
print(data_trans.columns)
print(sum(pd.isnull(data_trans['UPDASV139'])))
print(sum(pd.isnull(data_trans['Insulin'])))

data_trans['HBP'] = np.where(pd.isnull(data_trans['systolic']) | pd.isnull(data_trans['diastolic']), np.nan,
                           np.where( (data_trans['systolic'] >= 140) | (data_trans['diastolic'] >= 90), 1, 0))



data_org.to_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\microbiome_allvariable_806_original.csv', index=False)
data_trans.to_csv(r'C:\Users\yuqingw1\Workfolder\ProcessingData\microbiome_allvariable_806_transformation.csv', index=False)

#
#
# ### check
# olddata = pd.read_csv(path + '\Data\oldmicrobiome.csv')
# print(olddata.shape)
# print(len(set(olddata['sampleID'])))
# print(sum(pd.isnull(olddata['sampleID'])))
#
# newdata = pd.read_csv(path + '\ProcessingData\processing_4_microbiome\microbiome_combine_879row_493col.csv')
# print(newdata.shape)
# print(len(set(newdata['sampleID'])))
# print(sum(pd.isnull(newdata['sampleID'])))
#
#
# ### check the repeat ID
# from collections import Counter
# item_counts = Counter(list(newdata['sampleID']))
# # Finding items with more than one occurrence
# repeat_items = [item for item, count in item_counts.items() if count > 1]
# print("Repeated items:", repeat_items)
#
# repeat_in_old = olddata[np.in1d(olddata['sampleID'], repeat_items)]
# print(repeat_in_old.shape)
# print(len(repeat_items))