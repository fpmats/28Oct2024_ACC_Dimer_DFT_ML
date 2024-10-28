# Purpose of this script is to combine all the outputs from
# spherical_by_index.py and the later supplied SASAs csv file
# into one giant features csv file.

# Importing relevant packages
import pandas as pd
import os
import numpy as np

# Getting a list of files in current working directory.
ls = os.listdir('.')

# Extracting the files of interest into a single list.
A_data_files = [x for x in ls if x.startswith('A-')]
B_data_files = [x for x in ls if x.startswith('B-')]
data_files = A_data_files + B_data_files

# Initializing empty list of soon to be pandas DataFrames
df_list = []

# Looping through the features.
for feature in data_files:
    # Extracting the data.
    data = np.loadtxt(feature,skiprows=2)
    # Making a DataFrame from the feature (data is combined into a single column)
    df = pd.read_table(feature,header=1)
    # Extracting the columns from the DataFrame above.
    columns = df.columns[0].split(',')

    # Creation of the proper DataFrame.
    df = pd.DataFrame(data, columns=columns)
    # Removing the Snapshot Number column.
    df_no_idx = df.drop(columns=columns[0])
    # Adding the atoms of interest as a prefix for identification.
    df_ready = df_no_idx.add_prefix(feature)
    # Appending to list of DataFrames.
    df_list.append(df_ready)

# Combining all the DataFrames in the list of DataFrames.
master_df = pd.concat(df_list, ignore_index=False,axis=1)

# We need to add the data from SASAs.csv.
sasas = pd.read_csv('SASAs.csv',)
sasas = sasas.drop(sasas.columns[0], axis=1)
master_df = master_df.join(sasas)

# Writing the massive DataFrame to a csv.
master_df.to_csv('all_features.csv',index=False)
print('Done!')

