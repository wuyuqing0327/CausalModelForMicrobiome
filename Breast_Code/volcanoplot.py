import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
from scipy.stats import t


df = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\Data\Breast_Cancer\breastcancer_pls_data.csv')
vip = pd.read_csv(r'C:\Users\yuqingw1\Workfolder\result\Breast_Cancer\vip_matrix.csv')
vip.columns= ['MetabolomicsID', 'White_alone_not_Hispanic_or_Latino_rate_binary_vip',
       'high_school_and_lessthen_rate_binary_vip',
       'poverty_below_100_rate_binary_vip']

df['StudyID'] = df['StudyID'].astype(str)
df = df.set_index(['StudyID'])
data = df[list(vip['MetabolomicsID'])].T
data = data.reset_index()
data.rename(columns={'index':'MetabolomicsID'}, inplace=True)
data = data.merge(vip, on=['MetabolomicsID'])
data.shape




def plot_volcano(labeldata, data, label):

    # Step 1: Calculate mean values for each binary group
    labeldata['binary_label'] = np.where(labeldata[label] <= labeldata[label].quantile(1 / 3), 1, 0)
    log_fold_changes = []
    for metabolite in data['MetabolomicsID']:
        mean_group_0 = labeldata[labeldata['binary_label'] == 0][metabolite].mean()
        mean_group_1 = labeldata[labeldata['binary_label'] == 1][metabolite].mean()
        # Compute log2 fold change (handle division by zero)
        if mean_group_0 > 0 and mean_group_1 > 0:
            log_fc = np.log2(mean_group_1 / mean_group_0)
        else:
            log_fc = 0  # Avoid invalid log operations
        log_fold_changes.append(log_fc)

    # Step 2: Add log2 fold change and -log10(p-value) to the data
    data['p_value'] = np.where(data['%s_binary_vip' % label] >= 2, 0.03, 2)
    data['log2_fold_change'] = log_fold_changes
    data['neg_log_p'] = -np.log10(data['p_value'])

    # Step 3: Define thresholds for significance
    fold_change_threshold = 1.0  # Log2 fold change threshold
    p_threshold = 0.05  # p-value threshold

    # Determine significance
    data['status'] = 'Not Significant'
    data.loc[(data['log2_fold_change'] > fold_change_threshold) & (
                data['p_value'] < p_threshold), 'status'] = 'Upregulated'
    data.loc[(data['log2_fold_change'] < -fold_change_threshold) & (
                data['p_value'] < p_threshold), 'status'] = 'Downregulated'

    # Step 4: Plot the volcano plot
    plt.figure(figsize=(10, 7))

    # Scatter plot
    colors = {'Upregulated': 'magenta', 'Downregulated': 'blue', 'Not Significant': 'gray'}
    for status, color in colors.items():
        subset = data[data['status'] == status]
        plt.scatter(subset['log2_fold_change'], subset['neg_log_p'], label=status, color=color)

    # Add threshold lines
    plt.axvline(fold_change_threshold, color='black', linestyle='--')
    plt.axvline(-fold_change_threshold, color='black', linestyle='--')
    plt.axhline(-np.log10(p_threshold), color='black', linestyle='--')

    # Add labels and legend
    plt.title("Volcano Plot of Metabolomics Data")
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-log10(p-value)")
    plt.legend()
    plt.grid()
    plt.show()




def plot_volcano(labeldata, data, label, coef_col):

    labeldata['binary_label'] = np.where(labeldata[label] <= labeldata[label].quantile(1/3), 1, 0)
    binary_1 = labeldata[labeldata['binary_label']==1].index
    binary_0 = labeldata[labeldata['binary_label']==0].index
    data['binary_1_value'] = data[binary_1].mean(axis=1)
    data['binary_0_value'] = data[binary_0].mean(axis=1)

    # non-significant (grey)
    data['p-value'] = np.where(data['%s_binary_vip' %label]>=2, 0.03, 2)
    not_significant = data['p-value'] > 0.05
    grey_points = data[not_significant]

    # Binary_bin < 1 and significant (green)
    significant_pos = (data['binary_1_value']/data['binary_0_value'] < 1) & (data['p-value'] <= 0.05)
    green_points = data[significant_pos]

    # Binary_bin >= 1 and significant (red)
    significant_neg = (data['binary_1_value']/data['binary_0_value'] >=1) & (data['p-value'] <= 0.05)
    red_points = data[significant_neg]

    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.scatter(grey_points['binary_1_value']/grey_points['binary_0_value'], grey_points['p-value'],
                alpha=0.75, c='grey', label='Non-significant')
    plt.scatter(np.log2(green_points['binary_1_value']/green_points['binary_0_value']), green_points['p-value'],
                alpha=0.75, c='green', label='Significant (Response Rate below than 1)')
    plt.scatter(red_points['binary_1_value']/red_points['binary_0_value']), np.log2(red_points['p-value']),
                alpha=0.75, c='red', label='Significant (Response Rate greater than 1)')

    basecoef = min(data[label])
    basepvalue = data['p-value']

    # zero coefficient
    plt.axvline(basecoef, color="black", linestyle="dashed", linewidth=1.2, label="Zero Coefficient")

    #  Significance line (p-value = 0.05)
    plt.axhline(basepvalue, color="blue", linestyle="--", linewidth=1.2, label='p-value threshold (0.05)')

    # Labels and title
    plt.xlabel(label)
    plt.ylabel("P-value")
    plt.title('%s_Volcano_Plot' %label)

    # Customize plot
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.xticks(rotation=45, ha="right")

    # Optional legend (can be commented out)
    plt.legend(loc="upper right")
    # Display the plot
    plt.tight_layout()
    plt.show()


ylist = ['White_alone_not_Hispanic_or_Latino_rate',
       'high_school_and_lessthen_rate',
       'poverty_below_100_rate']

# Create volcano plots for each biomarker
for label in ylist:
    plot_volcano(data,  label, list(set(data.columns) - set(['MetabolomicsID'])))

