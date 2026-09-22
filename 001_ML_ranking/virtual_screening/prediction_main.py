import warnings
warnings.filterwarnings("ignore")
import pandas as pd 
import numpy as np
import sklearn
import xgboost as xgb
from sklearn.impute import KNNImputer
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_classif
from sklearn import model_selection, metrics
import joblib
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_validate
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
import random
from scipy.stats import uniform
from sklearn.model_selection import RandomizedSearchCV
import pickle

class PRED:
    def __init__(self):
        self.not_numeric=[]      
    def prepare_data(self,X):
        X = pd.read_csv(X, sep=';')
        return(X)
    def apply(self,final_model,X_full,prediction):        
        model = pickle.load(open(final_model, "rb"))
        y_pred = pd.DataFrame(model.predict(X_full))
        y_pred = y_pred.rename(columns={0: 'pred'})
        y_pred.to_csv(prediction, sep="\t")
        print('Successful Prediction saved as {}'.format(prediction))
PRED=PRED()        

