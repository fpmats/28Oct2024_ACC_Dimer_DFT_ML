import csv
import pandas as pd
import numpy as np
from sklearn.svm import LinearSVC
from sklearn.model_selection import train_test_split, cross_val_score,StratifiedKFold
from balanced_class import balanced_subsample
from sklearn.metrics import f1_score,confusion_matrix
from sklearn.feature_extraction.text import CountVectorizer
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from BorutaShap import BorutaShap, load_data
from statistics import mean, stdev


data = pd.read_csv("occupied_1000.csv")
accuracy=[]
f1=[]
X = data.drop('states', axis=1)
y = data['states']

for i in range(10):
    xs,ys=balanced_subsample(X,y)
    xdata, X_test, ydata, y_test = train_test_split(X, y, test_size=0.1 ,random_state=15+i)
    clf = RandomForestClassifier(bootstrap=False, max_depth=20, max_features=2,
                      min_samples_split=4, n_estimators=500, random_state=42)
    clf.fit(xs,ys)
    scores = cross_val_score(clf, xs, ys, cv=10)
    b=(sum(scores)/len(scores))
    accuracy.append(b)
    acc_mean=mean(accuracy)
    clf.fit(xdata,ydata)
    predictions = clf.predict(X_test)
    cm = confusion_matrix(y_test, predictions, labels=clf.classes_)
    #print('Score: ', clf.score(X_test, y_test))
    f1.append(f1_score(y_test,predictions,average="macro"))
    print(b)
acc_stdev=stdev(accuracy)
f1_stdev=stdev(f1)
print(acc_mean,acc_stdev,mean(f1),f1_stdev)

#accuracy=[]
#std=[]
#f1=[]
#for i in range(100):
#    data = pd.read_csv("occupied.csv")
#    accuracy=[]
#    X = data.drop('states', axis=1)
#    y = data['states']
#    og_oob_scores = []
#    fil_oob_scores = []
#    xs,ys=balanced_subsample(X,y)
#    X_train, X_test, y_train, y_test = train_test_split(xs, ys, test_size=0.1,random_state=i)
#    rf = RandomForestClassifier(n_estimators=100)
#    rf.fit(xs,ys)
#    scores= cross_val_score(rf, xs, ys, cv=10)
#    skf = StratifiedKFold(n_splits=10)
#    skf.get_n_splits(xs, ys)
#    StratifiedKFold(n_splits=10, random_state=None, shuffle=False)#

#    for train_index, test_index in skf.split(xs, ys):
#         X_train, X_test = xs.iloc[train_index], xs.iloc[test_index]
#         y_train, y_test = ys.iloc[train_index], ys.iloc[test_index]
#         rf.fit(X_train, y_train)
#         y_pred=rf.predict(X_test)
#         a=confusion_matrix(rf.predict(X_test),y_test)
#         y_pred=rf.predict(X_test)
#         a=confusion_matrix(rf.predict(X_test),y_test)
#         #print(a[0,0])
#         #print(a[1,0])
#         #print(a[0,1])
#         #print(a[1,1])
#         f1.append(f1_score(rf.predict(X_test),y_test))
#         print(sum(f1)/len(f1))
#    #print(scores[0])
#    #print(scores[1])
#    #print(scores[2])
#    #print(scores[3])
#    #print(scores[4])
#    #print(scores[5])
#    #print(scores[6])
#    #print(scores[7])
#    #print(scores[8])
#    #print(scores[9])
#    accuracy.append(sum(scores)/len(scores))
#    std.append(statistics.stdev(scores))
#    for i in range(len(accuracy)):
#        print(accuracy[i],std[i],f1[i])
