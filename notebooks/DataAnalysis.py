# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:26:24 2019
@author: ppradeep

#######################################################################################################
This script analyses the genetox data used in the TSCA effort to:
1. Derive a genetox classification for each chemical based on Williams at el. 
2. Evaluate the distribution of Ames and Clastogen assays for each chemical
#######################################################################################################

"""
#%%
#######################################################################################################
## Set working directory 
#######################################################################################################
import os
clear = lambda: os.system('cls')
clear()

## For manually inputting the working directory - for now using a fixed directory
#workingdir = input("Enter project working directory:")
#path = workingdir
#os.chdir(path)
path = 'C:/Users/ppradeep/OneDrive - Environmental Protection Agency (EPA)/Profile/Desktop/GeneTox/'
os.chdir(path)
#%%
#########################################################################################################  
## Import libraries
#########################################################################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#%%
#########################################################################################################
# 1. Read input files
#########################################################################################################  
genetox_data = pd.read_excel('data/genetox_160120.xlsx')

#%%
#########################################################################################################    
# 2. Derive DeMarini Calls - assign calls based on DeMarini logic
######################################################################################################### 
df_demarini = pd.DataFrame(index = genetox_data['casrn'].unique())
#%%
for i, group in genetox_data.groupby('casrn'):
    cas = group['casrn'].unique()
    dtxsid = group['dsstox_substance_id'].unique()
    # check which CASRNs are associated with multiple DTXSIDs
    if len(dtxsid)>1:
        print(cas)
    subset = group[['simple_aggregate', 'assay_outcome']]
    tuples = [tuple(x) for x in subset.values]
    if ('Ames', 1) in tuples:
        df_demarini.loc[cas, 'DeMarini_call']='gentox'
    elif ('clastogen', 1) in tuples:
        df_demarini.loc[cas, 'DeMarini_call']='clastogen'  
    elif ('Ames', 2) in tuples:
        df_demarini.loc[cas, 'DeMarini_call']='inconclusive'    
    elif ('clastogen', 2) in tuples:
        df_demarini.loc[cas, 'DeMarini_call']='inconclusive' 
    elif ('Ames', 0) in tuples:
        df_demarini.loc[cas, 'DeMarini_call']='non gentox'
    elif ('clastogen', 0) in tuples:
        df_demarini.loc[cas, 'DeMarini_call']='not clastogen'   
# below snippet for predictions of these Assays and not actual assay data from logic used in TSCA work
#    elif ('pAmes', 1) in tuples:
#        df_demarini.loc[cas, 'DeMarini_call']='pred gentox'
#    elif ('pclastogen', 1) in tuples:
#        df_demarini.loc[cas, 'DeMarini_call']='pred clastogen'
#    elif ('pAmes', 0) in tuples:
#        df_demarini.loc[cas, 'DeMarini_call']='pred non gentox'
#    elif ('pclastogen', 0) in tuples:
#        df_demarini.loc[cas, 'DeMarini_call']='pred non clastogen'
    else:
        None
    df_demarini.loc[cas, 'DTXSID']=dtxsid[0]

df_demarini.to_csv('output/DeMariniCalls.csv', index_label = 'CASRN')

#%%
#########################################################################################################    
# 3. Plot the distribution of DeMarini assignments (Figure not in manuscript)
######################################################################################################### 
# create a copy of the dataframe for further manipulations
demarini = df_demarini.copy()

plt.figure(figsize=(8, 6), dpi = 200)
y = demarini.DeMarini_call.value_counts()
ax = y.plot.barh()
for i, v in enumerate(y):
    ax.text(v-190, i + .3, str(v), color='purple', fontsize = 16, fontweight = 'bold')
plt.xticks(fontsize = 16, family = 'serif')
ticklabels = ['Non-mutagenic', 'Mutagenic', 'Clastogen', 'Inconclusive', 'Non-clastogen']
plt.yticks(range(0,5), ticklabels, fontsize = 16, family = 'serif')
plt.xlabel('Number of Chemicals', size = 24, labelpad = 10, family = 'serif')
plt.ylabel('Classification', size = 24, labelpad = -10, family = 'serif')
plt.savefig('output/ClassificationCount.png', bbox_inches='tight')

#%%
#########################################################################################################    
# 4. Generate data for Table 1 and 2 in manuscript 
######################################################################################################### 

## Convert DeMarini classifications to binary Genotoxicity (0/1) call
demarini.replace('non gentox', 0, inplace = True)
demarini.replace('not clastogen' , 0, inplace = True)
demarini.replace('gentox', 1, inplace = True)
demarini.replace('clastogen', 1, inplace = True)
demarini.replace('inconclusive', np.nan, inplace = True)

## Create a new dataframe for summary assay statistics per chemical
genetox_data_summary = pd.DataFrame(index = demarini.index, columns = ['nAmes', 'nAmesPos', 'nAmesNeg', 'nAmesInc', '%AmesPos', \
                                                                'nClas', 'nClasPos','nClasNeg', 'nClasInc', '%ClasPos', 'DeMarini'])
genetox_data_summary['DeMarini'] = demarini['DeMarini_call']
# Update each row in the dataframe
for casrn in genetox_data_summary.index:
    # Ames
    genetox_data_summary.loc[casrn, 'nAmes'] = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'Ames')])
    genetox_data_summary.loc[casrn, 'nAmesPos'] = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'Ames')\
                            & (genetox_data['assay_outcome'] == 1)])
    genetox_data_summary.loc[casrn, 'nAmesNeg'] = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'Ames')\
                            & (genetox_data['assay_outcome'] == 0)])
    genetox_data_summary.loc[casrn, 'nAmesInc'] = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'Ames')\
                            & (genetox_data['assay_outcome'] == 2)])
    if genetox_data_summary.loc[casrn, 'nAmes'] == 0:
        genetox_data_summary.loc[casrn, '%AmesPos'] = np.nan
    else:
        genetox_data_summary.loc[casrn, '%AmesPos'] = 100*genetox_data_summary.loc[casrn, 'nAmesPos']/genetox_data_summary.loc[casrn, 'nAmes']
    # Update inconclusives
    n1 = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'Ames')\
                            & (genetox_data['assay_outcome'] == 2)])
    if genetox_data_summary.loc[casrn, '%AmesPos'] == 0 and n1 > 0:
        genetox_data_summary.loc[casrn, '%AmesPos'] = -9 #inconclusive is indicated by -9
   # Clastogens
    genetox_data_summary.loc[casrn, 'nClas'] = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'clastogen')])
    genetox_data_summary.loc[casrn, 'nClasPos'] = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'clastogen')\
                            & (genetox_data['assay_outcome'] == 1)])
    genetox_data_summary.loc[casrn, 'nClasNeg'] = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'clastogen')\
                            & (genetox_data['assay_outcome'] == 0)])    
    genetox_data_summary.loc[casrn, 'nClasInc'] = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'clastogen')\
                            & (genetox_data['assay_outcome'] == 2)])        
    if genetox_data_summary.loc[casrn, 'nClas'] == 0:
        genetox_data_summary.loc[casrn, '%ClasPos'] = np.nan
    else:
        genetox_data_summary.loc[casrn, '%ClasPos'] = 100*genetox_data_summary.loc[casrn, 'nClasPos']/genetox_data_summary.loc[casrn, 'nClas']
    # Update inconclusives 
    n2 = len(genetox_data[(genetox_data['casrn'] == casrn) \
                            & (genetox_data['simple_aggregate'] == 'clastogen')\
                            & (genetox_data['assay_outcome'] == 2)])    
    if genetox_data_summary.loc[casrn, '%ClasPos'] == 0 and n2 > 0:
        genetox_data_summary.loc[casrn, '%ClasPos'] = -9 #inconclusive is indicated by -9

#%%
#########################################################################################################    
# 5. Get counts for Table 1. in manuscript (Analysis of chemicals tested in the Ames assay and  
# their classification as genotoxic, non-genotoxic or inconclusive using the DeMarini scheme)
######################################################################################################### 
#Ames
n10 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] != -9) & (genetox_data_summary['%AmesPos'] < 50)])
n20 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] != -9) & (genetox_data_summary['%AmesPos'] >= 50)]) 
n30 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] == -9)]) 
print("Count PercAmesPos <50 = %d" %n10)
print("Count PercAmesPos >=50 = %d" %n20)
print("Count Ames Inconclusives = %d" %n30)

#Ames and DeMarini negative
n11 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] != -9) & (genetox_data_summary['%AmesPos'] < 50) & (genetox_data_summary['DeMarini'] ==0)])
n21 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] != -9) & (genetox_data_summary['%AmesPos'] >= 50) & (genetox_data_summary['DeMarini'] ==0)])
n31 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] == -9) & (genetox_data_summary['DeMarini'] ==0)])
print("Count DeMarini negative and PercAmesPos <50 = %d" %n11)
print("Count DeMarini negative and PercAmesPos >=50 = %d" %n21)
print("Count DeMarini negative and Ames Inconclusives = %d" %n31)

#Ames and DeMarini positive
n12 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] != -9) & (genetox_data_summary['%AmesPos'] < 50) & (genetox_data_summary['DeMarini'] ==1)])
n22 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] != -9) & (genetox_data_summary['%AmesPos'] >= 50) & (genetox_data_summary['DeMarini'] ==1)])
n32 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] == -9) & (genetox_data_summary['DeMarini'] ==1)])
print("Count DeMarini positive and PercAmesPos <50 = %d" %n12)
print("Count DeMarini positive and PercAmesPos >=50 = %d" %n22)
print("Count DeMarini positive and Ames Inconclusives = %d" %n32)

#Ames and DeMarini Inconclusive
n13 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] != -9) & (genetox_data_summary['%AmesPos'] < 50) & (genetox_data_summary['DeMarini'].isnull())])
n23 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] != -9) & (genetox_data_summary['%AmesPos'] >= 50) & (genetox_data_summary['DeMarini'].isnull())])
n33 = len(genetox_data_summary[(genetox_data_summary['%AmesPos'] == -9) & (genetox_data_summary['DeMarini'].isnull())])
print("Count DeMarini Inconclusive and PercAmesPos <50 = %d" %n13)
print("Count DeMarini Inconclusive and PercAmesPos >=50 = %d" %n23)
print("Count DeMarini Inconclusive and Ames Inconclusives = %d" %n33)

#%%
#########################################################################################################    
# 6. Get counts for Table 2 in manuscript (ACharacterization of number of chemicals tested 
# in less than K (1, 2, 3, 4, 5 and 10) or less number of Ames or Clastogen assays)
######################################################################################################### 
for k in range(11):
    a = genetox_data_summary[genetox_data_summary['nAmes'] <= k]['nAmes'].count()
    c = genetox_data_summary[genetox_data_summary['nClas'] <= k]['nClas'].count()
    print("Number of chemicals tested in <= %d Ames assays = %d" %(k,a))
    print("Number of chemicals tested in <= %d Clastogen assays = %d" %(k,c))
