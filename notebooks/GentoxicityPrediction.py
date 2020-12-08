# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 00:17:40 2020

@author: ppradeep

#######################################################################################################
This script analyses implements an ensemble model to predict genotoxicity for a test chemical
Inputs: Prediction from EPA TEST, Lazar, and OECD alerts
Output: Genotoxicity classification
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

###########################################################################
## Import libraries
###########################################################################
import pandas as pd

#%%
###########################################################################
## Define functions
###########################################################################


#%%
#######################################################################################################
# 1. Read the original data file
#######################################################################################################
genetox_data = pd.read_excel('data/genetox_160120.xlsx')

#%%
######################################################################################################
# 2. Load posterior probability distribution from available data and set the cut-off value
######################################################################################################
# combination = 'comb2_5-T1T2A1A3A4'
final_model =  pd.read_csv('output/FinalModel.csv', index_col='Combination')
cut_off = 0.25

#%%
######################################################################################################
# 3. Input data for new chemical and predict 
######################################################################################################
print("Enter the tool predictions for the test chemical")
t1 = int(input("EPA TEST Prediction:"))
t2 = int(input("EPA TEST Prediction:"))
a1 = int(input("OECD Alert 1 Prediction:"))
a3 = int(input("OECD Alert 3 Prediction:"))
a4 = int(input("OECD Alert 4 Prediction:"))

# calculate the combination number 
comb = 2**0*t1+2**1*t2+2**2*a1+2**3*a3+2**4*a4

# compare the posterior probability with the cut-off and make a prediction
try:
   post_prob = final_model.loc[comb][0]
   if post_prob>cut_off:
       print("The chemical is genotoxic")
   else:
        print("The chemical is non-genotoxic")
except:
    print("The chemical is out of domain for prediction")
    