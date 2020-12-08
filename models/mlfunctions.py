# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 22:20:53 2020

@author: ppradeep
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 12:19:40 2018

@author: ppradeep

----------------------------------------------------------------------------------------
This script contains all the user-defined functions related to machine learning for use 
in this project.
----------------------------------------------------------------------------------------

"""
#%%
## Import libraries
import pandas as pd
import numpy as np

# Classifiers
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.svm import SVR

# Machine learning relevant 
import sklearn.metrics as sm
from sklearn.feature_selection import SelectPercentile, f_classif
from sklearn.feature_selection import RFE
from sklearn.feature_selection import VarianceThreshold 
from sklearn import preprocessing
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import GridSearchCV

#%%
## User-defined functions
def selectFeatures_VarThresh(X, threshold):
    sel = VarianceThreshold(threshold=(threshold * (1 - threshold))) 
    X_sel = sel.fit_transform(X)
    # Convert it into a dataframe 
    x_tr = pd.DataFrame(X_sel, index = X.index) 
    x_tr.columns = X.columns[sel.get_support(indices = True)]
    return x_tr

## Remove culumns with >80% correlation
def correlation(dataset, threshold):
    col_corr = set() # Set of all the names of deleted columns
    corr_matrix = dataset.corr()
    for i in range(len(corr_matrix.columns)):
        for j in range(i):
            if corr_matrix.iloc[i, j] >= threshold:
                colname = corr_matrix.columns[i] # getting the name of column
                col_corr.add(colname)
                if colname in dataset.columns:
                    del dataset[colname] # deleting the column from the dataset
    return dataset

# Normalize descriptors: Transform variables to mean=0, variance=1
def normalizeDescriptors(X):
    scaler = preprocessing.StandardScaler().fit(X)
    transformed = scaler.transform(X)
    x_norm = pd.DataFrame(transformed, index = X.index) 
    x_norm.columns = X.columns
    return(scaler, x_norm)

   
# Univariate feature selection with F-test for feature scoring
def selectFeatures_perc(X, Y, percentile):
    model = SelectPercentile(f_classif, percentile)
    model = model.fit(X, Y) #convert datatype for use in the fit function
    scores = -np.log10(model.pvalues_)
    scores /= scores.max()
    X_tr = model.transform(X)
    ## Convert it into a dataframe 
    X_tr = pd.DataFrame(X_tr, index = X.index) 
    X_tr.columns = X.columns[model.get_support(indices=True)] 
    return X_tr

def selectFeatures_RFE(X, Y, n_features_to_select):
    model = LinearRegression()
    rfe = RFE(model, n_features_to_select )
    rfe = rfe.fit(X, Y) #convert datatype for use in the fit function
    X_tr = rfe.transform(X)
    ## Convert it into a dataframe 
    X_tr = pd.DataFrame(X_tr, index = X.index) 
    X_tr.columns = X.columns[rfe.get_support(indices=True)] 
    return X_tr

def train_knn(n_fold, X, Y):
    parameters = {'weights':['uniform', 'distance'], 'n_neighbors':[5, 6, 7, 8, 9, 10, 15], 'algorithm': ['auto', 'ball_tree', 'kd_tree', 'brute']}
    clf = KNeighborsRegressor()
    grid_search = GridSearchCV(clf, cv = n_fold, param_grid = parameters)
    grid_search.fit(X, Y)
    knn_params = grid_search.best_params_
    knn_score = grid_search.best_score_
    print(knn_params)
    return(grid_search, knn_score, knn_params)

def train_svm(n_fold, X, Y):
    parameters = {'kernel':['linear', 'rbf'], 'C':[0.1, 1, 10], 'gamma':[0.01, 0.1, 1], 'epsilon': [0.1, 1]}
    clf = SVR()
    grid_search = GridSearchCV(clf, cv = 5, param_grid = parameters)
    grid_search.fit(X, Y)
    svm_params = grid_search.best_params_
    svm_score = grid_search.best_score_
    print(svm_params)
    return(grid_search, svm_score, svm_params)
    
def train_rf(n_fold, X, Y):
    parameters = {'n_estimators': [500, 750, 1000, 1500, 2000], 'max_features': ['auto', 'sqrt'], 'random_state':[5]} #'sqrt'
    clf = RandomForestRegressor()
    grid_search = GridSearchCV(clf, cv = n_fold, param_grid = parameters)
    grid_search.fit(X, Y)
    rf_params = grid_search.best_params_
    rf_score = grid_search.best_score_
    print(rf_params)
    return(grid_search, rf_score, rf_params)
   
def predict_y(clf, X, Y, n_fold):
    y = cross_val_predict(clf, X = X, y = Y, cv = n_fold)
    return y

def predict_test_y(clf, X, Y, X_test):
    clf = clf.fit(X,Y)
    y = clf.predict(X_test)
    return y

## Calculate metrics
# Classification
def calcMetrics(cnf_matrix):
    tn, fp, fn, tp = cnf_matrix.ravel()
    total = float(tp + tn + fp + fn)
    acc = round(100*float(tp + tn)/float(total),2)
    sens = round(100*float(tp)/float(tp + fp),2)
    spec = round(100*float(tn)/float(tn + fn),2)
    ba = round((sens+spec)/2,2)
    p_o = float(tp + tn)/total
    p_e = ((tp + fn)/total)*((tp + fp)/total) + ((fp + tn)/total)*((fn + tn)/total)
    kappa = round(((p_o - p_e)/(1 - p_e)), 2)
    return total, acc, sens, spec, ba, kappa
# Regression
def reg_metrics(true, pred):
    rmse = np.round(np.sqrt(sm.mean_squared_error(true, pred)),2)
    r2 = np.round(sm.r2_score(true, pred),2)
    return(rmse, r2)
  
#%%
## Clustering related functions
from sklearn.cluster import KMeans 
from scipy.spatial.distance import cdist, pdist
import matplotlib.pyplot as plt 

def n_clusters(X, nRange, stepsize):
    k_range = range(stepsize, nRange, stepsize) # Number of clusters
    k_means_var = [KMeans(init='k-means++', n_clusters = k).fit(X) for k in k_range]
    centroids = [v.cluster_centers_ for v in k_means_var]
    # determine cluster distances    
    k_euclid = [cdist(X, cent, 'euclidean') for cent in centroids]
    dist = [np.min(ke, axis=1) for ke in k_euclid]
    wcss = [sum(d**2) for d in dist]
    tss = sum(pdist(X)**2)/X.shape[0]
    bss = tss - wcss
    #percentage variance explained
    perc_var = bss/tss*100
    ## Number of clusters chosen kIdx
    for i, val in enumerate(perc_var):
        if val > 50:
            kIdx = i
            print(kIdx, val)
            break
    return(k_range, bss, tss, kIdx, val)
    
def fitKmeans(X, k):
    kmeans = KMeans(init='k-means++', n_clusters = k) 
    kmeans.fit(X) 
    return kmeans

def clusterboxplot(df_kmeans, k, y_var, count):
    df = df_kmeans
    upper_quartile = df.groupby(by='Cluster').quantile(q=0.75)
    ax = df_kmeans.boxplot(column = y_var, by = 'Cluster', showfliers = False, figsize=[20,8])
    ## Annotate each boxplot with number of chemicals in each cluster
    for x in upper_quartile.index:
        ax.annotate(count[x], (x+1,upper_quartile.loc[x,y_var]+0.05), fontsize = 14, color = 'r')
  
    plt.title("") 
    plt.suptitle("")  
    xtick_label = range(1,k+1)
    plt.xticks(range(1,k+1),xtick_label, fontsize = 22, rotation = 'vertical')
    plt.yticks(fontsize = 22)
    plt.ylabel('LogD (pH7.4)',  size = 30, labelpad = 10)
    plt.xlabel('Cluster Number',  size = 30, labelpad = 20)
    plt.grid(False)
    return plt


    