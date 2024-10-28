from pprint import pprint
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
import pandas as pd
import numpy as np

def evaluate(model, test_features, test_labels):
    predictions = model.predict(test_features)
    errors = abs(predictions - test_labels)
    mape = 100 * np.mean(errors / test_labels)
    accuracy = 100 - mape
    print('Model Performance')
    print('Average Error: {:0.4f} degrees.'.format(np.mean(errors)))
    print('Accuracy = {:0.2f}%.'.format(accuracy))

    return accuracy


df = pd.read_csv("mulliken_2.csv")
X = df.drop('states', axis=1)
y = df['states']
rf = RandomForestRegressor(random_state = 42)

base_model = RandomForestRegressor(n_estimators = 100, random_state = 42)
base_model.fit(X, y)
base_accuracy = evaluate(base_model,X, y)

print('Parameters currently in use:\n')
pprint(rf.get_params())

param_grid = {
    'bootstrap': [True,False],
    'max_depth': [15,20,60,80,None],
    'max_features': [2, 3, 4, 6,"sqrt","log2"],
    'min_samples_leaf': [2, 3, 4, 6,1],
    'min_samples_split': [2, 3, 4, 6],
    'n_estimators': [100, 200, 300, 500]}

grid_search = GridSearchCV(estimator = rf, param_grid = param_grid, cv = 10, n_jobs = -1, verbose = 2)
# Fit the grid search to the data
grid_search.fit(X, y)
grid_search.best_params_

best_grid = grid_search.best_estimator_
grid_accuracy = evaluate(best_grid, X, y)
print(best_grid)
print('Improvement of {:0.2f}%.'.format( 100 * (grid_accuracy - base_accuracy) / base_accuracy))
