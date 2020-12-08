# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 16:26:24 2019
@author: ppradeep

#######################################################################################################
This script analyses the genetox data used in the TSCA effort to evaluate the variability in 
genetox classification for each chemical based on Williams at el. using a sampling analysis
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
from sklearn import metrics 
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

#%%
#######################################################################################################
# 1. Read the original data file
#######################################################################################################
genetox_data = pd.read_excel('data/genetox_160120.xlsx')

#%%
#######################################################################################################
# 1. Bootstrapping: Sample 1 ames and 1 clastogen activity for each chemical from the experimental 
# data (underlying distribution). 
# - For each sample derive the Genetox classification as per the Williams et al. scheme
# - Calculate the overall variance/accuracy of agreement as a measure of reproducibility
# - Calculate the unalikability index as a measure of overall agreeability of experimental 
#   values/chemical
#######################################################################################################
n_bootstraps = 500
AmesSubsetforBootstrap = pd.DataFrame(index=[], \
                                      columns=['CASRN','AmesOutcomeArraynoInc', '%AmesPos', '%AmesNeg','AmesUnalikeability',\
                                               'AmesBootstrappedArray', '%bAmesPos', '%bAmesNeg','AmesBUnalikeability',\
                                               'ClasOutcomeArraynoInc', '%ClasPos', '%ClasNeg','ClasUnalikeability',\
                                               'ClasBootstrappedArray', '%bClasPos', '%bClasNeg','ClasBUnalikeability',\
                                               'DeMariniBootstrapped'])
idx = 0
for cas, group in genetox_data.groupby('casrn'):
    #cas = group['casrn'].unique()
    AmesOutcomeArray = np.array(group[group['simple_aggregate']=='Ames']['assay_outcome'])
    ClasOutcomeArray = np.array(group[group['simple_aggregate']=='clastogen']['assay_outcome'])
    # Remove all the inconclusive outcomes
    AmesOutcomeArraynoInc = AmesOutcomeArray[~np.isin(AmesOutcomeArray,(2))]
    ClasOutcomeArraynoInc = ClasOutcomeArray[~np.isin(ClasOutcomeArray,(2))]
    if len(AmesOutcomeArraynoInc) > 2:
        # create a bootstraped value array
        AmesBootstrapvalues = np.array(np.random.choice(AmesOutcomeArraynoInc, replace=True, size=n_bootstraps))
        try:
            ClasBootstrapvalues = np.array(np.random.choice(ClasOutcomeArraynoInc, replace=True, size=n_bootstraps))
        except:
            ClasBootstrapvalues = np.array([])
        # Calculate the % of 0,1 and 2 (non-mutagenic, mutagenic and inconclusives) in original data
        AmesPercPos = 100*(AmesOutcomeArraynoInc==1).sum()/len(AmesOutcomeArraynoInc)
        AmesPercNeg = 100*(AmesOutcomeArraynoInc==0).sum()/len(AmesOutcomeArraynoInc)
        Amesunalikeability = 2*AmesPercPos*AmesPercNeg/10000 #keep it as a number between 0 and 1
        #PercInc = 100*(AmesOutcomeArraynoInc==2).sum()/len(AmesOutcomeArraynoInc)
        # Calculate the % of 0,1 and 2 (clastogenicity) in original data
        ClasPercPos = 100*(ClasOutcomeArraynoInc==1).sum()/len(ClasOutcomeArraynoInc)
        ClasPercNeg = 100*(ClasOutcomeArraynoInc==0).sum()/len(ClasOutcomeArraynoInc)
        Clasunalikeability = 2*ClasPercPos*ClasPercNeg/10000 #keep it as a number between 0 and 1
        
        # Calculate the % of 0,1 and 2 (non-mutagenic, mutagenic and inconclusives) in bootsstrapped data
        AmesBPercPos = 100*(AmesBootstrapvalues==1).sum()/n_bootstraps
        AmesBPercNeg = 100*(AmesBootstrapvalues==0).sum()/n_bootstraps
        AmesBunalikeability = 2*AmesBPercPos*AmesBPercNeg/10000 #keep it as a number between 0 and 1
        #bPercInc = 100*(bootstrapvalues==2).sum()/n_bootstraps
        # Calculate the % of 0,1 and 2 (clastogenicity) in bootsstrapped data
        ClasBPercPos = 100*(ClasBootstrapvalues==1).sum()/n_bootstraps
        ClasBPercNeg = 100*(ClasBootstrapvalues==0).sum()/n_bootstraps
        ClasBunalikeability = 2*ClasBPercPos*ClasBPercNeg/10000 #keep it as a number between 0 and 1
               
        # Derive DeMarini Call for each bootstrap
        bootstrapdemarini = []
        for i, ames in enumerate(AmesBootstrapvalues):
            if ames == 1:
                bootstrapdemarini.append(1)
            else:
                try:
                    clas = ClasBootstrapvalues[i]
                    if clas == 1:
                        bootstrapdemarini.append(1)
                    else:
                        bootstrapdemarini.append(0)   
                except:
                    bootstrapdemarini.append(0)
            
        # Update the DataFrame
        AmesSubsetforBootstrap.loc[idx,'CASRN'] = cas
        #Ames
        AmesSubsetforBootstrap.loc[idx, 'AmesOutcomeArraynoInc'] = AmesOutcomeArraynoInc
        AmesSubsetforBootstrap.loc[idx, '%AmesPos']= AmesPercPos
        AmesSubsetforBootstrap.loc[idx, '%AmesNeg'] = AmesPercNeg
        AmesSubsetforBootstrap.loc[idx, 'AmesUnalikeability'] = Amesunalikeability
        AmesSubsetforBootstrap.loc[idx, 'AmesBootstrappedArray'] = AmesBootstrapvalues
        AmesSubsetforBootstrap.loc[idx, '%bAmesPos']= AmesBPercPos
        AmesSubsetforBootstrap.loc[idx, '%bAmesNeg'] = AmesBPercNeg    
        AmesSubsetforBootstrap.loc[idx, 'AmesBUnalikeability'] = float(AmesBunalikeability)
        # Clas
        AmesSubsetforBootstrap.loc[idx, 'ClasOutcomeArraynoInc'] = ClasOutcomeArraynoInc
        AmesSubsetforBootstrap.loc[idx, '%ClasPos']= ClasPercPos
        AmesSubsetforBootstrap.loc[idx, '%ClasNeg'] = ClasPercNeg
        AmesSubsetforBootstrap.loc[idx, 'ClasUnalikeability'] = Clasunalikeability
        AmesSubsetforBootstrap.loc[idx, 'ClasBootstrappedArray'] = ClasBootstrapvalues
        AmesSubsetforBootstrap.loc[idx, '%bClasPos']= ClasBPercPos
        AmesSubsetforBootstrap.loc[idx, '%bClasNeg'] = ClasBPercNeg    
        AmesSubsetforBootstrap.loc[idx, 'ClasBUnalikeability'] = float(ClasBunalikeability)
        # DeMarini calls
        AmesSubsetforBootstrap.loc[idx, 'DeMariniBootstrapped'] = bootstrapdemarini
        # Update the index
        idx=idx+1

AmesSubsetforBootstrap.index = AmesSubsetforBootstrap['CASRN']
AmesSubsetforBootstrap.to_csv('output/AmesBootstrapAnalysis.csv', index_label = 'CASRN')     

#%%
#######################################################################################################
# 2. Read DeMarini Calls
#######################################################################################################
df_demarini = pd.read_csv('output/DeMariniCalls.csv', index_col = 'CASRN')['DeMarini_call']
# Convert DeMarini classifications to binary Genotoxicity (0/1) call
demarini = df_demarini.copy()
demarini.replace('non gentox', 0, inplace = True)
demarini.replace('not clastogen' , 0, inplace = True)
demarini.replace('gentox', 1, inplace = True)
demarini.replace('clastogen', 1, inplace = True)
demarini.replace('inconclusive', np.nan, inplace = True)

#######################################################################################################
# 3. Calculate equivalence|accuracy of bootstrap DeMarini against original DeMarini calls
#######################################################################################################
evaluationdataset = pd.concat([AmesSubsetforBootstrap['DeMariniBootstrapped'], demarini], axis = 1).dropna()
AmesSubsetforBootstrap['DeMarini'] = demarini#['DeMarini_call']
bootstrapacc = []
for idx in evaluationdataset.index:
    x = np.empty(n_bootstraps)
    x.fill(evaluationdataset.loc[idx, 'DeMarini_call'])
    y = evaluationdataset.loc[idx, 'DeMariniBootstrapped']
    acc = 100*np.round(metrics.accuracy_score(x, y),2)
    evaluationdataset.loc[idx, 'DeMariniBootstrappedAccuracy'] = acc

#%%
#######################################################################################################
# 4. Plots
#######################################################################################################

# a. Plot Unalikeability
#-----------------------
plt.figure(figsize=(8, 6), dpi = 200)
Amestotaldata = len(AmesSubsetforBootstrap['AmesBUnalikeability'].dropna())
AmesSubsetforBootstrap['AmesBUnalikeability'].astype('float64').hist(bins=20,  weights=np.ones(Amestotaldata)/(Amestotaldata), color='g', alpha=0.5, label='Ames')

Clastotaldata = len(AmesSubsetforBootstrap['ClasBUnalikeability'].dropna())
AmesSubsetforBootstrap['ClasBUnalikeability'].astype('float64').hist(bins=20,  weights=np.ones(Clastotaldata)/(Clastotaldata), color='r', alpha=0.5, label='Clastogen')

plt.legend(prop={'size':18,'family':'serif'})
plt.xlabel('Unalikeability', fontsize=18, family = 'serif')
plt.ylabel('Fraction of Data', fontsize=18, family = 'serif')
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.savefig('output/unalikeability.png', bbox_inches='tight')

# b. Plot reproducibility
#------------------------
plt.figure(figsize=(8, 6), dpi = 200)
plt.grid(True)
totaldata = len(evaluationdataset['DeMariniBootstrappedAccuracy'])
plt.hist(evaluationdataset['DeMariniBootstrappedAccuracy'], bins=20,  weights=np.ones(totaldata)/(totaldata))
#plt.legend(prop={'size':18,'family':'serif'})
plt.xlabel('Reproducibility', fontsize=18, family = 'serif')
plt.ylabel('Frequency', fontsize=18, family = 'serif')
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.gca().xaxis.set_major_formatter(PercentFormatter(100))
plt.savefig('output/demariniReproducibility.png', bbox_inches='tight')
#%%