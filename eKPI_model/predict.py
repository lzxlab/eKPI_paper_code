import pandas as pd
from sklearn.naive_bayes import GaussianNB,MultinomialNB
from sklearn.datasets import load_digits 
from sklearn.model_selection import train_test_split
from sklearn import metrics
import pandas as pd 
from sklearn.decomposition import PCA 
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.svm import SVC
from sklearn.feature_selection import VarianceThreshold
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
import os
import sys
import numpy as np
from sklearn.utils import shuffle
import joblib
from sklearn.metrics import roc_auc_score
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score,classification_report,confusion_matrix
from sklearn.model_selection import KFold, cross_val_score,cross_validate
from imblearn.over_sampling import SMOTE
from sklearn.svm import SVC
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier
import xgboost as xgb
from sklearn.metrics import make_scorer
from xgboost import plot_importance
import dask.dataframe as dd
import warnings

def del_rows(data):
    threshold = int(2/3 * data.shape[1])
    non_zero_counts = data.isna().sum(axis=1)
    data_filtered = data[non_zero_counts <= threshold]
    return data_filtered

kinase_family = sys.argv[1]
input_file = sys.argv[2]
output_path = sys.argv[3]

# read input file
dask_df = dd.read_csv(input_file)

# load model
XGB = joblib.load('./Model/train_model_XGB_'+kinase_family+'.m')

# impute missing values
dask_df = dask_df.fillna(0.0)

# convert KPS column to string
dask_df['KPS'] = dask_df['KPS'].astype(str)
dask_df = dask_df.set_index('KPS')

# predict
predicted = XGB.predict(dask_df.compute())
y_scores = XGB.predict_proba(dask_df.compute())[:, 1]

# save result
result_df = pd.DataFrame({
    'KPS': dask_df.index.compute(),
    'prob': y_scores,
    'predicted_label': predicted
})

result_df.to_csv(os.path.join(output_path, 'result.csv'), index=False)

import time 

with open(output_path+"/runInfo.txt","a") as f:
    f.write("Dear user,\n\nThanks for your interest in eKPI!\n\nYour task has compeleted.\n\n")
    f.write(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
