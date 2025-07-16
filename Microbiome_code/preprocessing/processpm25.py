import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import rasterio
from rasterio.plot import show
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from statsmodels.graphics.correlation import plot_corr
import os

folder_path = '/Users/yuqingwu/UChicago_MSDS/Part-time/BSD/Datasets/PM25/pm25_2/'
csv_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.csv')]
csv_files.sort()  # Optional: Sort the files if needed

merged_df = pd.DataFrame()
for file in csv_files:
    # Read the current CSV file into a DataFrame
    current_df = pd.read_csv(file)
    current_df = current_df[~pd.isnull(current_df['specimen_id'])]
    current_df = current_df.reset_index(drop=True)
    current_df['specimen_id'] = current_df['specimen_id'].astype(int).astype(str)
    current_df['lon'] = current_df['lon'].round(6)
    current_df['lat'] = current_df['lat'].round(6)
    current_df = current_df.groupby(['specimen_id', 'lon', 'lat'], as_index=False).mean()
    print("shape:", current_df.shape)
    # If merged_df is empty, copy current_df to it; else, merge current_df into merged_df
    if merged_df.empty:
        merged_df = current_df
    else:
        merged_df = merged_df.merge(current_df, on=['specimen_id', 'lat', 'lon'])

print("merged_df shape:", merged_df.shape)
merged_df = merged_df.drop_duplicates(['specimen_id'])
merged_df = merged_df.reset_index(drop=True)
print("merged_df shape:", merged_df.shape)


#
# ######## check data ########
# compass_all = pd.read_csv('/Users/yuqingwu/UChicago_MSDS/Part-time/BSD/Datasets/compass_all_address.csv')
# compass_all = compass_all[~pd.isnull(compass_all['specimen_id'])]
# compass_all = compass_all.reset_index(drop=True)
# compass_all['specimen_id'] = compass_all['specimen_id'].astype(int).astype(str)
# print(compass_all[np.in1d(compass_all['specimen_id'], ['106525', '106881'])])
# print(compass_all.shape)
# compass_all = compass_all.drop_duplicates(['specimen_id'])
# compass_all = compass_all.reset_index(drop=True)
# print("compass_all shape:", compass_all.shape)
#
# ### check the repeat ID
# from collections import Counter
# item_counts = Counter(list(compass_all['specimen_id']))
# # Finding items with more than one occurrence
# repeat_items = [item for item, count in item_counts.items() if count > 1]
# print("basicdata Repeated items:", repeat_items)
# repeatid = ['100271', '101281', '101437', '101261', '101578', '102509', '102597', '102676', '103377', '104366', '104406', '103986', '104159', '104726', '105083', '105688', '106060', '106539', '106394', '106588', '106525', '106615', '106766', '106881', '107037', '106503', '107556', '107863', '108093', '107936', '108196', '108332', '108359', '108344', '108375']
# print('104277' in repeatid)
# print('104390' in repeatid)
# print('106701' in repeatid)
# print('104301' in repeatid)


compass_deidentified = pd.read_csv('/Users/yuqingwu/UChicago_MSDS/Part-time/BSD/Datasets/compass_deidentified.csv')
compass_deidentified = compass_deidentified[~pd.isnull(compass_deidentified['specimen_id'])]
compass_deidentified['specimen_id'] = compass_deidentified['specimen_id'].astype(int).astype(str)

# ### check the repeat ID
# from collections import Counter
# item_counts = Counter(list(compass_deidentified['specimen_id']))
# # Finding items with more than one occurrence
# repeat_items = [item for item, count in item_counts.items() if count > 1]
# print("basicdata Repeated items:", repeat_items)
# repeatid2 = ['100271', '101281', '101437', '101261', '101578', '102597', '102676', '103377', '104366', '104406', '103986', '104159', '104726', '105083', '105688', '106060', '106394', '106588', '106525', '106615', '106766', '106881', '107037', '106503', '107863', '108093', '107936', '108332', '108359', '108344']
# print(set(repeatid) - set(repeatid2)) # ['102509', '106539', '107556', '108196', '108375']

# df_repeat = compass_deidentified[(compass_deidentified['specimen_id']=='104726') | (compass_deidentified['specimen_id']=='106503') | (compass_deidentified['specimen_id']=='106615') | (compass_deidentified['specimen_id']=='108359')]


### replace value
compass_deidentified = compass_deidentified.drop_duplicates(['specimen_id'], keep='first')
compass_deidentified.loc[compass_deidentified['specimen_id']=='106503', 'age'] = (53.36 + 56.77)/2
print(compass_deidentified.shape)
select_cols = ['specimen_id', 'date', 'age', 'bmi', 'gender', 'smoking_status', 'diabetes2',
               'systolic', 'systolic_1', 'systolic_2', 'diastolic', 'diastolic_1', 'diastolic_2']


compass_ap_seqn = merged_df.merge(compass_deidentified[select_cols], how='left')
print("compass_ap_seqn.shape:", compass_ap_seqn.shape)
print("check null whether exists:", sum(pd.isnull(compass_ap_seqn['age'])), sum(pd.isnull(compass_ap_seqn['bmi'])), sum(pd.isnull(compass_ap_seqn['gender'])))
compass_ap_seqn['date2'] = pd.to_datetime(compass_ap_seqn['date']).dt.strftime('%Y%m')
compass_ap_seqn = compass_ap_seqn[~pd.isnull(compass_ap_seqn['date'])]
compass_ap_seqn = compass_ap_seqn.reset_index(drop=True)
print("compass_ap_seqn.shape2:", compass_ap_seqn.shape)
print("check null whether exists2:", sum(pd.isnull(compass_ap_seqn['age'])), sum(pd.isnull(compass_ap_seqn['bmi'])), sum(pd.isnull(compass_ap_seqn['gender'])))


from datetime import datetime
from dateutil.relativedelta import relativedelta

datesunique = list(set(compass_ap_seqn['date2']))
basicdata = pd.DataFrame()
for i in datesunique:
    currentdata = compass_ap_seqn[compass_ap_seqn['date2'] == i][select_cols+['date2']]
    print("length:", currentdata.shape)
    enrollmonth_str = str(list(currentdata['date2'])[0])
    date_format = "%Y%m"
    enroll_month = datetime.strptime(enrollmonth_str, date_format)
    max_col = 'X_' + str(enrollmonth_str)
    all_columns = compass_ap_seqn.columns.tolist()
    for j in range(1, 37):
        months_to_subtract = j
        min_date = enroll_month - relativedelta(months=months_to_subtract)
        min_col = 'X_' + str(min_date.strftime('%Y%m'))
        filtered_columns = [col for col in all_columns if min_col < col <= max_col]
        # Select the filtered columns from the DataFrame
        currentdata['M%s' %j] = np.array(compass_ap_seqn[compass_ap_seqn['date2'] == i][filtered_columns].mean(axis=1))
    if basicdata.empty:
        basicdata = currentdata
    else:
        basicdata = pd.concat([basicdata, currentdata], axis=0)

print(basicdata.shape)
print("check null whether exists:", sum(pd.isnull(basicdata['age'])), sum(pd.isnull(basicdata['bmi'])), sum(pd.isnull(basicdata['gender'])))

### check the repeat ID
from collections import Counter
item_counts = Counter(list(basicdata['specimen_id']))
repeat_items = [item for item, count in item_counts.items() if count > 1]
print("basicdata Repeated items:", repeat_items)

basicdata.to_csv('Datasets/compass_pm25value_basic_detail.csv', index=False)

