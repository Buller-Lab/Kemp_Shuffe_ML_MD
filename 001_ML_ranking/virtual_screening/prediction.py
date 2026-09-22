import pandas as pd 
from pandas import MultiIndex, Int16Dtype
import numpy as np
import sklearn
import xgboost as xgb
from prediction_main import PRED

#path to trained model
final_model='../training_shuffle_dataset/kemp_model_24.pickle'

#path to features 
X = 'feature_generation/features_combinatorial_clean.csv'

#path to selected features
selected_features = '../training_shuffle_dataset/Top_24_selected_features.tsv'

#path to resulting prediction file
prediction = 'y_preds_combinatorial.tsv'


X=PRED.prepare_data(X)
df_sel_features=pd.read_csv(selected_features,sep='\t')
selected_features=df_sel_features.values.flatten().tolist()
del selected_features[0]
selected_features = [str(x) for x in selected_features]
X_selected = X.loc[:, selected_features]
PRED.apply(final_model,X_selected,prediction)