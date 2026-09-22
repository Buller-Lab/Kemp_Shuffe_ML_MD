import pandas as pd 
from pandas import MultiIndex, Int16Dtype
import numpy as np
import sklearn
import xgboost as xgb
from featureselector_main import featureselector


X_data = 'X_clean.csv'
y_data = 'y.csv'
final_model="prelim_model_kemp.json"
eval_=pd.DataFrame([['r2','r2_stdv','rmse','rmse_stdv']])
eval_.to_csv('feature_selection_monitor_kemp.csv',sep=',', header=False, index=False)
for number_features in range(1,100):
    X,y=featureselector.prepare_data(X_data,y_data)
    selected_features,X_selected=featureselector.select_features(X,y,number_features)
    r2_rmse = featureselector.train(X_selected,y,final_model)
    eval_df = pd.DataFrame([r2_rmse])
    eval_df.to_csv('feature_selection_monitor_kemp.csv', mode='a',sep=',', header=False, index=False)
