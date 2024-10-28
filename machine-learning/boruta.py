import csv
import pandas as pd
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.model_selection import train_test_split, cross_val_score,StratifiedKFold
#from dna_featuregenerator import create_feature_vectors
from balanced_class import balanced_subsample
from sklearn.metrics import f1_score,confusion_matrix
from sklearn.feature_extraction.text import CountVectorizer
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from BorutaShap import BorutaShap, load_data
a=[]
readfile= "ACC-hex_dataset.csv"

#create data structure for SVM
data = pd.read_csv(readfile)
accuracy=[]
X = data.drop('states', axis=1)
y = data['states']
xs,ys=balanced_subsample(X,y)
og_oob_scores = []
fil_oob_scores = []
rf = RandomForestClassifier(bootstrap=False, max_depth=20, max_features=2,
                      min_samples_split=4, n_estimators=500, random_state=42)
fs = BorutaShap(model=rf, importance_measure='shap', classification=True)
fs.fit(xs, ys, n_trials=100, sample=False, train_or_test = 'train', normalize=True, verbose=True)
fs.results_to_csv()
df1= pd.DataFrame(pd.read_csv('feature_importance.csv'))

for i in range(100):
    data = pd.read_csv(readfile)
    accuracy=[]
    X = data.drop('states', axis=1)
    y = data['states']
    xs,ys=balanced_subsample(X,y)
    og_oob_scores = []
    fil_oob_scores = []
    rf = RandomForestClassifier(bootstrap=False, max_depth=20, max_features=2,
                      min_samples_split=4, n_estimators=500, random_state=42)
    fs = BorutaShap(model=rf, importance_measure='shap', classification=True)
    fs.fit(xs, ys, n_trials=100, sample=False, train_or_test = 'train', normalize=True, verbose=True)
    fs.results_to_csv()
    df2= pd.DataFrame(pd.read_csv('feature_importance.csv'))
    df1= pd.concat([df1,df2])

df1.to_excel('final_mulliken_boruta.xlsx')

