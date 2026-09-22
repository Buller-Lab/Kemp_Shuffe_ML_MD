import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import KFold, cross_val_score

class featureselector:
    def __init__(self):
        self.not_numeric=[]      
    def prepare_data(self,X,y):
        X = pd.read_csv(X, sep=';')
        X = X.iloc[: , 1:]
        y = pd.read_csv(y, sep=';')
        return(X,y)
    def select_features(self,X_full,y,number_features):
        clf = xgb.XGBRegressor()
        model = clf.fit(X_full, y)
        importance = model.feature_importances_
        df_importance = pd.DataFrame([importance])
        df_importance = pd.DataFrame(data=df_importance.values,columns=X_full.columns).sort_values(0, axis=1, ascending=False)
        Top_X  = df_importance.iloc[: , :number_features]
        selected_features = list(Top_X.columns)
        df_sel_features = pd.DataFrame([selected_features])
        df_sel_features.to_csv('Top_{}_selected_features.tsv'.format(number_features), sep="\t")
        X_selected = X_full.loc[:, selected_features]
        print("Selected Top {} Features:".format(number_features))
        return(selected_features,X_selected)
    def train(self,X,y,final_model):
        clf = xgb.XGBRegressor()
        model=clf.fit(X, y)
        cv = KFold(n_splits=10, shuffle=True, random_state=0)
        rmse = cross_val_score(model, X, y, scoring='neg_root_mean_squared_error', cv=cv)
        print("RMSD: {} +/- {}".format(rmse.mean(), rmse.std()))         
        r2 = cross_val_score(model, X, y, scoring='r2', cv=cv)
        print("R2: {} +/- {}".format(r2.mean(), r2.std()))
        return([r2.mean(),r2.std(),rmse.mean(), rmse.std()])
featureselector=featureselector()        
