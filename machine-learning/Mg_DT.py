import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn.tree import DecisionTreeClassifier
import dtreeviz

# Load your custom dataset from Excel
df = pd.read_excel('final_physical.xlsx')

# Assume 'twins' is the target column and the rest are features
X = df.drop('states', axis=1)
y = df['states']

# Initialize lists to store accuracy and F1 scores
accuracy = []
f1 = []

# Initialize and fit the decision tree classifier
clf = DecisionTreeClassifier(max_depth=5, min_samples_leaf=8,random_state=42)

# Calculate cross-validated scores for accuracy
scores = cross_val_score(clf, X, y, cv=10)
print("Cross-validated accuracy scores:", scores)
b = (sum(scores) / len(scores))
print("Average accuracy:", b)
accuracy.append(b)

# Calculate cross-validated scores for F1 score
f1_scores = cross_val_score(clf, X, y, cv=10, scoring='f1_macro')
print("Cross-validated F1 scores:", f1_scores)
c = (sum(f1_scores) / len(f1_scores))
print("Average F1 score:", c)
f1.append(c)

# Fit the classifier on entire dataset to visualize the tree
clf.fit(X, y)

# Feature names array
feats = X.columns.to_numpy()

# Visualize the decision tree with dtreeviz
viz = dtreeviz.model(clf,X_train=X, y_train=y,
               target_name="states",
               feature_names=feats,
               class_names=['greater', 'fewer'])
v = viz.view()
v.show()
v.save("decision_tree_visualization_5_3.svg")
