# Purpose of this script is to take the giant list of features
# and splitting them into types of interactions that are
# examined.

# Importing pandas
import pandas as pd

# Creating pandas DataFrame from the features csv
DF = pd.read_csv('/path/to/feature/csv')

# Creating lists containing the column names for each category of feature
charges_and_dipoles = [
    *[x for x in DF.columns if 'CA' in x and 'NZ' in x],
    *[x for x in DF.columns if 'HZ' in x and 'NZ' in x],
    *[x for x in DF.columns if 'CE' in x and 'NZ' in x],
    *[x for x in DF.columns if 'CA' in x and 'CD' in x],
    *[x for x in DF.columns if 'CG' in x and 'CD' in x],
    *[x for x in DF.columns if 'OE' in x and 'CD' in x],
    *[x for x in DF.columns if '1-N' in x and 'CD' in x],
    *[x for x in DF.columns if '1-N' in x and 'CA' in x],
    *[x for x in DF.columns if 'OE' in x and 'HT' in x],
    *[x for x in DF.columns if 'NZ' in x and 'CD' in x],
    *[x for x in DF.columns if 'HZ' in x and 'OE' in x],
    *[x for x in DF.columns if 'CA' in x and 'LYS-29' in x and 'NZ' not in x],
]

backbone = [
    *[x for x in DF.columns if 'HN' in x and 'O' in x],
    *[x for x in DF.columns if 'Ramachandran' in x],
    ]

hydrophobic_core = [
    *[x for x in DF.columns if 'CA' in x and 'ILE' in x],
    *[x for x in DF.columns if 'CA' in x and 'PHE' in x],
    *[x for x in DF.columns if 'PHE' in x and 'CA' not in x and 'HN' not in x and 'Ramachandran' not in x and 'SASA' not in x],
]

sasa = [
    *[x for x in DF.columns if 'SASA' in x],
]

# Creating smaller csv files for the machine-learning
names = ['charges_and_dipoles','backbone','hydrophobic_core','sasa']
my_list = [charges_and_dipoles,backbone,hydrophobic_core,sasa]
for i in range(len(my_list)):
    DF.to_csv(f'{names[i]}.csv',columns=my_list[i], index=False)

