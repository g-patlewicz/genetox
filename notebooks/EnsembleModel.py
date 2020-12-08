# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:41:38 2020

@author: ppradeep

#######################################################################################################
This script analyses the genetox data used in the TSCA effort to:
1. Calculate the predicitive ability of QSAR tools and OECD alerts against Williams at el. calls
2. Develop and validate an ensemble model for genotoxicity
#######################################################################################################

"""

#%%
###########################################################################
## Set working directory 
###########################################################################
import os
clear = lambda: os.system('cls')
clear()
path = 'Y:/Projects/GeneTox/'
path = 'C:/Users/ppradeep/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/GeneTox/'
os.chdir(path)

#%%
###########################################################################
## Import libraries
###########################################################################
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix
import MLfunctions as ml
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations
from sklearn.model_selection import KFold


#%%
#######################################################################################################
# 1. Read the original data file
#######################################################################################################
genetox_data = pd.read_excel('data/genetox_160120.xlsx')

#%%
#######################################################################################################
# 2. Read DeMarini Calls
#######################################################################################################
#-------------------------------------------------------------------------------#
## 1. Read DeMarini Calls
#-------------------------------------------------------------------------------#
df_demarini = pd.read_csv('output/DeMariniCalls.csv', index_col = 'CASRN')
# Convert the calls to binary 0/1
demarini = df_demarini.copy()
demarini.replace('non gentox', 0, inplace = True)
demarini.replace('not clastogen' , 0, inplace = True)
demarini.replace('gentox', 1, inplace = True)
demarini.replace('clastogen', 1, inplace = True)
demarini.replace('inconclusive', np.nan, inplace = True)

#%%
#######################################################################################################
# 3. Create master dictionary of smiles, casrn and dtxsid 
#######################################################################################################
smiles = pd.read_csv(path+'data/smiles.csv')
smiles.index = smiles['CASRN']
smiles['DTXSID'] = df_demarini['DTXSID']


#%%
#######################################################################################################
# 4. Read Predictions and Alerts
####################################################################################################### 
metrics_tools=pd.DataFrame()
#-------------------------------------------------------------------------------#
# 4.1 EPA TEST Predictions (T1)
#-------------------------------------------------------------------------------#
df_test = pd.read_csv('data/TEST_Batch_Mutagenicity_Consensus.txt', sep = '|')[['ID','Pred_Result']]
df_test.rename(columns={'ID':'CASRN','Pred_Result':'TEST'}, inplace = True)
df_test.index = df_test['CASRN']

df_test.replace('Mutagenicity Positive', 1, inplace = True)
df_test.replace('Mutagenicity Negative', 0, inplace = True)

df_test['DTXSID'] = smiles['DTXSID']
df_test.drop_duplicates(inplace=True)

# Calculate Prediction Metrics - compare against DeMarini calls
compare = pd.concat([demarini['DeMarini_call'], df_test['TEST']], axis = 1, sort = True).dropna()
cnf_matrix = confusion_matrix(compare['DeMarini_call'], compare['TEST'])
total, acc, sens, spec, ba, kappa = ml.calcMetrics(cnf_matrix)
metrics_tools.loc['EPA TEST|Williamsetal', 'Truth'] = 'Williamns et al. 2019 Calls'
metrics_tools.loc['EPA TEST|Williamsetal', 'Total'] = total 
metrics_tools.loc['EPA TEST|Williamsetal', 'Accuracy'] = acc 
metrics_tools.loc['EPA TEST|Williamsetal', 'Sensitivity'] = sens 
metrics_tools.loc['EPA TEST|Williamsetal', 'Specificity'] = spec 
metrics_tools.loc['EPA TEST|Williamsetal', 'Balanced Accuracy'] = ba 
metrics_tools.loc['EPA TEST|Williamsetal', 'Kappa'] = kappa
# Rename Column ##
df_test.rename(columns={'TEST':'T1'}, inplace = True)

#%%
#-------------------------------------------------------------------------------#
# 4.2 LAZAR Predictions (T2)
#-------------------------------------------------------------------------------#
lazar1 = pd.read_csv('data/2020-01-14_lazar_batch_prediction_smiles1.csv')
lazar1.index = lazar1.CASRN

lazar2 = pd.read_csv('data/2020-01-14_lazar_batch_prediction_smiles2.csv')
lazar2.index = lazar2.CASRN

# Combine and convert into binary predictions 
df_lazar = pd.concat([lazar1, lazar2])
df_lazar.replace('mutagenic', 1, inplace = True)
df_lazar.replace('non-mutagenic', 0, inplace = True)

df_lazar = df_lazar[['CASRN', 'Prediction']]
df_lazar['DTXSID'] = smiles['DTXSID']
df_lazar.rename(columns={'Prediction': 'Lazar'}, inplace = True)

df_lazar.drop_duplicates(inplace=True)

## Calculate Prediction Metrics - compare with DeMarini calls
compare = pd.concat([demarini['DeMarini_call'], df_lazar['Lazar']], axis = 1, sort = True).dropna()
cnf_matrix = confusion_matrix(compare['DeMarini_call'], compare['Lazar'])
total, acc, sens, spec, ba, kappa = ml.calcMetrics(cnf_matrix)
metrics_tools.loc['Lazar|Williamsetal', 'Truth'] = 'Williamns et al. 2019 Calls'
metrics_tools.loc['Lazar|Williamsetal', 'Total'] = total 
metrics_tools.loc['Lazar|Williamsetal', 'Accuracy'] = acc 
metrics_tools.loc['Lazar|Williamsetal', 'Sensitivity'] = sens 
metrics_tools.loc['Lazar|Williamsetal', 'Specificity'] = spec 
metrics_tools.loc['Lazar|Williamsetal', 'Balanced Accuracy'] = ba 
metrics_tools.loc['Lazar|Williamsetal', 'Kappa'] = kappa
# Rename Column ##
df_lazar.rename(columns={'Lazar':'T2'}, inplace = True)

#%%
#-------------------------------------------------------------------------------#
# 4.3 OECD Toolbox Predictions (Save the oecd.csv as tab delimited text file) 
#-------------------------------------------------------------------------------#
df_oecd = pd.read_csv('data/oecd.txt', sep = '\t', encoding = "ISO-8859-1")
#df_oecd.rename(columns={'CAS Number': 'CASRN'}, inplace = True)

for idx in df_oecd.index:
    cas = df_oecd.loc[idx,'CASRN']
    df_oecd.loc[idx,'DTXSID'] = smiles[smiles['CASRN']==cas]['DTXSID'].values

df_oecd.index = df_oecd['CASRN']
df_oecd = df_oecd[['DNA alerts for AMES by OASIS',\
            'DNA alerts for CA and MNT by OASIS',\
            'Protein binding alerts for Chromosomal aberration by OASIS',\
            'in vitro mutagenicity (Ames test) alerts by ISS',\
            'in vivo mutagenicity (Micronucleus) alerts by ISS']]

df_oecd['DNA alerts for AMES by OASIS'] = df_oecd['DNA alerts for AMES by OASIS'].apply(lambda x: 0 if x == 'No alert found' else 1)
df_oecd['DNA alerts for CA and MNT by OASIS'] = df_oecd['DNA alerts for CA and MNT by OASIS'].apply(lambda x: 0 if x == 'No alert found' else 1)
df_oecd['Protein binding alerts for Chromosomal aberration by OASIS'] = df_oecd['Protein binding alerts for Chromosomal aberration by OASIS'].apply(lambda x: 0 if x == 'No alert found' else 1)
df_oecd['in vitro mutagenicity (Ames test) alerts by ISS'] = df_oecd['in vitro mutagenicity (Ames test) alerts by ISS'].apply(lambda x: 0 if x == 'No alert found' else 1)
df_oecd['in vivo mutagenicity (Micronucleus) alerts by ISS'] = df_oecd['in vivo mutagenicity (Micronucleus) alerts by ISS'].apply(lambda x: 0 if x == 'No alert found' else 1)

# Rename alerts A1, A2 etc ##
df_oecd.rename(columns={'DNA alerts for AMES by OASIS':'A1',\
                     'DNA alerts for CA and MNT by OASIS':'A2',\
                     'Protein binding alerts for Chromosomal aberration by OASIS':'A3',\
                     'in vitro mutagenicity (Ames test) alerts by ISS': 'A4',\
                     'in vivo mutagenicity (Micronucleus) alerts by ISS': 'A5',}, inplace = True)

# Drop duplciated CASRNs
df_oecd = df_oecd.reset_index().drop_duplicates(subset='CASRN', keep='last').set_index('CASRN')

# Calculate Prediction Metrics - compare with DeMarini calls
compare = pd.concat([demarini['DeMarini_call'], df_oecd], axis = 1, sort = True).dropna()
cnf_matrix = confusion_matrix(compare['DeMarini_call'], compare['A1'])
total, acc, sens, spec, ba, kappa = ml.calcMetrics(cnf_matrix)
metrics_tools.loc['OECD A1|Williamsetal', 'Truth'] = 'Williamns et al. 2019 Calls'
metrics_tools.loc['OECD A1|Williamsetal', 'Total'] = total 
metrics_tools.loc['OECD A1|Williamsetal', 'Accuracy'] = acc 
metrics_tools.loc['OECD A1|Williamsetal', 'Sensitivity'] = sens 
metrics_tools.loc['OECD A1|Williamsetal', 'Specificity'] = spec 
metrics_tools.loc['OECD A1|Williamsetal', 'Balanced Accuracy'] = ba 
metrics_tools.loc['OECD A1|Williamsetal', 'Kappa'] = kappa

cnf_matrix = confusion_matrix(compare['DeMarini_call'], compare['A2'])
total, acc, sens, spec, ba, kappa = ml.calcMetrics(cnf_matrix)
metrics_tools.loc['OECD A2|Williamsetal', 'Truth'] = 'Williamns et al. 2019 Calls'
metrics_tools.loc['OECD A2|Williamsetal', 'Total'] = total 
metrics_tools.loc['OECD A2|Williamsetal', 'Accuracy'] = acc 
metrics_tools.loc['OECD A2|Williamsetal', 'Sensitivity'] = sens 
metrics_tools.loc['OECD A2|Williamsetal', 'Specificity'] = spec 
metrics_tools.loc['OECD A2|Williamsetal', 'Balanced Accuracy'] = ba 
metrics_tools.loc['OECD A2|Williamsetal', 'Kappa'] = kappa

cnf_matrix = confusion_matrix(compare['DeMarini_call'], compare['A3'])
total, acc, sens, spec, ba, kappa = ml.calcMetrics(cnf_matrix)
metrics_tools.loc['OECD A3|Williamsetal', 'Truth'] = 'Williamns et al. 2019 Calls'
metrics_tools.loc['OECD A3|Williamsetal', 'Total'] = total 
metrics_tools.loc['OECD A3|Williamsetal', 'Accuracy'] = acc 
metrics_tools.loc['OECD A3|Williamsetal', 'Sensitivity'] = sens 
metrics_tools.loc['OECD A3|Williamsetal', 'Specificity'] = spec 
metrics_tools.loc['OECD A3|Williamsetal', 'Balanced Accuracy'] = ba 
metrics_tools.loc['OECD A3|Williamsetal', 'Kappa'] = kappa

cnf_matrix = confusion_matrix(compare['DeMarini_call'], compare['A4'])
total, acc, sens, spec, ba, kappa = ml.calcMetrics(cnf_matrix)
metrics_tools.loc['OECD A4|Williamsetal', 'Truth'] = 'Williamns et al. 2019 Calls'
metrics_tools.loc['OECD A4|Williamsetal', 'Total'] = total 
metrics_tools.loc['OECD A4|Williamsetal', 'Accuracy'] = acc 
metrics_tools.loc['OECD A4|Williamsetal', 'Sensitivity'] = sens 
metrics_tools.loc['OECD A4|Williamsetal', 'Specificity'] = spec 
metrics_tools.loc['OECD A4|Williamsetal', 'Balanced Accuracy'] = ba 
metrics_tools.loc['OECD A4|Williamsetal', 'Kappa'] = kappa

cnf_matrix = confusion_matrix(compare['DeMarini_call'], compare['A5'])
total, acc, sens, spec, ba, kappa = ml.calcMetrics(cnf_matrix)
metrics_tools.loc['OECD A5|Williamsetal', 'Truth'] = 'Williamns et al. 2019 Calls'
metrics_tools.loc['OECD A5|Williamsetal', 'Total'] = total 
metrics_tools.loc['OECD A5|Williamsetal', 'Accuracy'] = acc 
metrics_tools.loc['OECD A5|Williamsetal', 'Sensitivity'] = sens 
metrics_tools.loc['OECD A5|Williamsetal', 'Specificity'] = spec 
metrics_tools.loc['OECD A5|Williamsetal', 'Balanced Accuracy'] = ba 
metrics_tools.loc['OECD A5|Williamsetal', 'Kappa'] = kappa

metrics_tools.to_csv('output/metrics_tools.csv', index_label = 'Tool|Truth Combination')

#%%
#######################################################################################################
# 5. Simple Bayesian combinatorial model for 2 tools and 5 alerts 
####################################################################################################### 

#-------------------------------------------------------------------------------#
## 5.1. Combine all the tools into a single dataframe and check the correlation
#-------------------------------------------------------------------------------#
x = pd.concat([df_test['T1'], df_lazar['T2'], df_oecd], axis = 1, sort=True)

# Evaluate and plot correlation
corr_pearson = x.corr(method='pearson')
plt.figure(figsize=(12, 8), dpi = 200)
mask = np.zeros_like(corr_pearson)
mask[np.triu_indices_from(mask)] = True
with sns.axes_style("white"):
    ax = sns.heatmap(corr_pearson, mask=mask, vmax=1, square=True, cmap="Blues", 
                     annot=True, cbar_kws={'label': 'Pearson Correlation'}, annot_kws={"size": 20})
ax.figure.axes[-1].yaxis.label.set_size(24)    
plt.xticks(fontsize = 24)
plt.yticks(fontsize = 24)
plt.savefig(path+'/output/ToolsCorr.png', bbox_inches='tight')    
#sns.cubehelix_palette(8)

#------------------------------------------------------------------------------------------#
# 5.2 Remove columns with >80% correlation to be used as descriptors in the ensemble model
#------------------------------------------------------------------------------------------#
x = ml.correlation(x, 0.80)

#------------------------------------------------------------------------------------------#
# 5.3. Calculate combination for all chemicals (Combination 1 and 2)
#------------------------------------------------------------------------------------------#
ntools = len(x.columns)
x_combns = pd.DataFrame(index=x.index)
for n in range(2,ntools+1):
    cc = np.array(list(combinations(x.columns, n)))
    for c in cc:
        # combination 1
        x_combns['comb1_%d-%s'%(n,''.join(c))] = x[c].sum(axis=1)
        # combination 2
        s = 0
        for p, tool in enumerate(c):
            s=s+2**p*x[tool]
        x_combns['comb2_%d-%s'%(n,''.join(c))] = s

#------------------------------------------------------------------------------------------#
# 5.4 Save the list of combinations
#------------------------------------------------------------------------------------------#

combns = x_combns.columns.tolist()
#x_combns.dropna(inplace=True) #remove all null data

#%%
#------------------------------------------------------------------------------------------#
# 5.5 5-fold CV Predictions
#------------------------------------------------------------------------------------------#
# Divide into test and training datasets
kf = KFold(n_splits=5) # Define the split - into 5 folds 
kf.get_n_splits(x_combns) # returns the number of splitting iterations in the cross-validator

# Make the predictions
predictions_x_combns = pd.DataFrame()
for train_index, test_index in kf.split(x_combns):
    train_set, test_set = x_combns.iloc[train_index], x_combns.iloc[test_index]
    #Calculate posteriors from training set
    post_prob = {}
    for combn in combns:
        for comb in train_set[combn].dropna().unique():
            try:
                n_comb = len(train_set[train_set[combn]==comb])
                n_1_comb = len(train_set[(train_set[combn]==comb) & (train_set['demarini']==1)])
                post_prob[comb] = round(n_1_comb/n_comb,2)
            except:
                post_prob[comb] = np.nan
                
        # Calculate new predictions for test set
        for idx in test_set.index:
            comb = test_set.loc[idx, combn]
            for prob_cutoff in [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55]:
                try:
                    if post_prob[comb] >= prob_cutoff:
                        predictions_x_combns.loc[idx, '%s_PredCutOff%0.2f' %(combn,prob_cutoff)] = 1
                    else:
                        predictions_x_combns.loc[idx, '%s_PredCutOff%0.2f' %(combn,prob_cutoff)] = 0
                except:
                    if comb==0:
                        predictions_x_combns.loc[idx, '%s_PredCutOff%0.2f' %(combn,prob_cutoff)] = 0
                    elif comb==sum([2**x for x in range(1,ntools+1)]):
                        predictions_x_combns.loc[idx, '%s_PredCutOff%0.2f' %(combn,prob_cutoff)] = 1   
                    else:
                        predictions_x_combns.loc[idx, '%s_PredCutOff%0.2f' %(combn,prob_cutoff)] = np.nan #'Cannot predict'    


#---------------------------------------------------------------------------------------------#
# 5.6. Calculate Prediction Metrics 
#---------------------------------------------------------------------------------------------#
x_combns['demarini'] = demarini['DeMarini_call']
metrics_x_combns = pd.DataFrame(columns = ['EnsembleCombination', 'Truth', 'Prediction_CutOff', 'nTools', 'Total', 'Accuracy', 'Sensitivity', 'Specificity', 'Balanced Accuracy', 'Kappa']) 
for combn in combns:
    # vary over cut-offs used in prediction
    for prob_cutoff in [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55]:
        xy = pd.concat([predictions_x_combns['%s_PredCutOff%0.2f' %(combn,prob_cutoff)], x_combns['demarini']], axis = 1).dropna()
        cnf_matrix = confusion_matrix(xy['%s_PredCutOff%0.2f' %(combn,prob_cutoff)],xy['demarini'])
        n = int(combn.split('_')[1][0])
        total, acc, sens, spec, ba, kappa = ml.calcMetrics(cnf_matrix)
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'EnsembleCombination'] = combn
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'Truth'] = 'Williamns et al. 2019 Calls'
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'Prediction_CutOff'] = prob_cutoff
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'nTools'] = 'n=%d' %n
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'Total'] = total
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'Accuracy'] = acc
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'Sensitivity'] = sens
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'Specificity'] = spec
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'Balanced Accuracy'] = ba
        metrics_x_combns.loc['%s|%s' %('%s_PredCutOff%0.2f' %(combn,prob_cutoff), 'DeMarini'), 'Kappa'] = kappa
metrics_x_combns.to_csv('output/metrics_combns.csv', index_label = 'Combination')


#######################################################################################################
# 6. Select the best cut-off/metrics for each combination and plot
#######################################################################################################
metrics_toplot = pd.DataFrame(columns=metrics_x_combns.columns)
for combn in combns:
    metrics_toplotsubset1 = metrics_x_combns[metrics_x_combns['EnsembleCombination']==combn]
    metrics_toplotsubset2 = metrics_toplotsubset1[metrics_toplotsubset1['Balanced Accuracy']==metrics_toplotsubset1['Balanced Accuracy'].max()]
    metrics_toplotsubset3 = metrics_toplotsubset2[metrics_toplotsubset2['Prediction_CutOff']==metrics_toplotsubset2['Prediction_CutOff'].max()]
    metrics_toplot=metrics_toplot.append(metrics_toplotsubset3)


#---------------------------------------------------------------------------------------------#
# 6.1 Plot the Prediction Metrics
#---------------------------------------------------------------------------------------------#
# Set the labels for x-axis - tool combinations
xaxis_label_dup = [combn.split('-')[1] for combn in combns]
xaxis_label = xaxis_label_dup[1::2]

# Combination 1
metrics_x_combns1 = metrics_toplot[metrics_toplot['EnsembleCombination'].str.contains('comb1')]
plt.figure(figsize=(20, 8), dpi = 200)
metrics_x_combns1.nTools=metrics_x_combns1.nTools.astype('category')

ax = sns.scatterplot(x=range(len(metrics_x_combns1)), y = 'Balanced Accuracy', hue='nTools', data=metrics_x_combns1, marker = 's', s=100, label = 'Balanced Accuracy', legend = False)
ax = sns.scatterplot(x=range(len(metrics_x_combns1)), y = 'Sensitivity', hue='nTools', data=metrics_x_combns1, marker = '+', s=100, label = 'Sensitivity', legend = False)
ax = sns.scatterplot(x=range(len(metrics_x_combns1)), y = 'Specificity', hue='nTools', data=metrics_x_combns1, marker = 'x', s=100, label = 'Specificity')

metrics_x_combns1['Balanced Accuracy'] = metrics_x_combns1['Balanced Accuracy'].astype('float')
maxBA = metrics_x_combns1.loc[metrics_x_combns1['Balanced Accuracy'].idxmax()]['Balanced Accuracy']
pcutoff = metrics_x_combns1.loc[metrics_x_combns1['Balanced Accuracy'].idxmax()]['Prediction_CutOff']
maxAccPos = metrics_x_combns1.index.get_loc(metrics_x_combns1['Balanced Accuracy'].idxmax())
plt.plot(maxAccPos,maxBA,marker='o', markersize=20, color = 'k', fillstyle='none')
plt.axvline(maxAccPos, ls = '--', c = 'k', alpha = 0.25)
txt = '%0.1f (%0.2f)' %(maxBA, pcutoff)
plt.annotate(txt, [maxAccPos+0.25,maxBA+3.25], fontsize = 16)

plt.xticks(range(len(metrics_x_combns1)),xaxis_label, fontsize = 16, rotation = 90)
plt.yticks(fontsize = 22)
plt.xlabel('Tool Combination', size = 30)
plt.ylabel('Percentage', size = 30, labelpad = -5)
plt.legend(fontsize=14, loc='lower right')
plt.ylim([0,100])
plt.xlim([-1,len(metrics_x_combns1)])
plt.grid(which='major', axis='y', alpha=0.5)
plt.savefig('output/MetricsCombn1.png', bbox_inches='tight')

# Combination 2
metrics_x_combns2 = metrics_toplot[metrics_toplot['EnsembleCombination'].str.contains('comb2')]
plt.figure(figsize=(20, 8), dpi = 200)
metrics_x_combns2.nTools=metrics_x_combns2.nTools.astype('category')

ax = sns.scatterplot(x=range(len(metrics_x_combns2)), y = 'Balanced Accuracy', hue='nTools', data=metrics_x_combns2, marker = 's', s=100, label = 'Balanced Accuracy', legend = False)
ax = sns.scatterplot(x=range(len(metrics_x_combns2)), y = 'Sensitivity', hue='nTools', data=metrics_x_combns2, marker = '+', s=100, label = 'Sensitivity', legend = False)
ax = sns.scatterplot(x=range(len(metrics_x_combns2)), y = 'Specificity', hue='nTools', data=metrics_x_combns2, marker = 'x', s=100, label = 'Specificity')

metrics_x_combns2['Balanced Accuracy'] = metrics_x_combns2['Balanced Accuracy'].astype('float')
maxBA = metrics_x_combns2.loc[metrics_x_combns2['Balanced Accuracy'].idxmax()]['Balanced Accuracy']
pcutoff = metrics_x_combns2.loc[metrics_x_combns2['Balanced Accuracy'].idxmax()]['Prediction_CutOff']
maxAccPos = metrics_x_combns2.index.get_loc(metrics_x_combns2['Balanced Accuracy'].idxmax())
plt.plot(maxAccPos,maxBA,marker='o', markersize=20, color = 'k', fillstyle='none')
plt.axvline(maxAccPos, ls = '--', c = 'k', alpha = 0.25)
txt = '%0.1f (%0.2f)' %(maxBA, pcutoff)
plt.annotate(txt, [maxAccPos-5,maxBA+3.25], fontsize = 16)

plt.xticks(range(len(metrics_x_combns2)),xaxis_label, fontsize = 16, rotation = 90)
plt.yticks(fontsize = 22)
plt.xlabel('Tool Combination', size = 30)
plt.ylabel('Percentage', size = 30, labelpad = -5)
plt.legend(fontsize=14, loc='lower right')
plt.ylim([0,100])
plt.xlim([-1,len(metrics_x_combns2)])
plt.grid(which='major', axis='y', alpha=0.5)
plt.savefig('output/MetricsCombn2.png', bbox_inches='tight')

#%%
#######################################################################################################
# 7. Build the final model on complete data
########################################################################################################
final_comb = 'comb2_5-T1T2A1A3A4'
final_comb_post_prob = {}

for comb in x_combns[final_comb].dropna().unique():
    try:
        n_comb = len(x_combns[x_combns[final_comb]==comb])
        n_1_comb = len(x_combns[(x_combns[final_comb]==comb) & (x_combns['demarini']==1)])
        final_comb_post_prob[comb] = round(n_1_comb/n_comb,2)
    except:
        final_comb_post_prob[comb] = np.nan


# Convert dictionary to dataframe and save
final_model = pd.DataFrame.from_dict(final_comb_post_prob, orient='index', columns=['Probability'])
final_model.to_csv('output/FinalModel.csv', index_label='Combination')