import pandas as pd 
from pandas import MultiIndex, Int16Dtype
import numpy as np
import sklearn
import xgboost as xgb
from PRTCI_main import PRTCI

#path to features + real Tms (.csv file) for training
X_data = 'X_clean.csv'
y_data = 'y.csv'

#number of features for feature selection
number_features=24
#path to model
final_model="kemp_model_24.pickle"

print('preparing data (60 seconds)...')
X,y=PRTCI.prepare_data(X_data,y_data)
print('selecting features...')
selected_features,X_selected=PRTCI.select_features(X,y,number_features)
print('training model...')
PRTCI.train(X_selected,y,final_model)

