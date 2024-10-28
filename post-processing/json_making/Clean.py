# The purpose of this script is to take the output from
# Read.sh and clean it up such that each line is an orbital

# Importing
import os
from glob import iglob

# Inputs: rootdir_glob = directory pointing at *_Gaussian.txt files
#         new_tail = suffix for the output of this script

rootdir_glob = f'/projectnb/fpmats/CCAM/conductivity_workflow/ACC_Dimer/3_post-processing/json_making/*_Gaussian.txt'
new_tail = '_Gaussian_cleaned.txt'

# Creates a list of paths to process.
file_list = [f for f in iglob(rootdir_glob, recursive=True) if os.path.isfile(f)]

# Looping over every path in the file_list.
for file in file_list:
    # Reading in the file.
    with open(file,'r') as text:
        # Removing any virtual orbitals and any "Atomic contribution" lines.
        orbital = [line.replace('\n','') for line in text if not "vir" in line and not "Atomic" in line]
    # This is to get the indices where we have a new orbital starting
    beginnings = [x for x in range(len(orbital)) if 'occ' in orbital[x]]
    # For the function of the next loop, we append None to this list.
    # Without the None, we miss the last occupied orbital.
    beginnings.append(None)
    # looping to create the one-liners for each orbital
    for i in range(len(beginnings)-1):
        orbital[beginnings[i]]=orbital[beginnings[i]:beginnings[i+1]]
    # Getting rid of strings not created by the above (we created a list of string(s) for each orbital).
    orbital = [x for x in orbital if type(x)==list]
    # Looping over this new list of lists for each orbital.
    for i in range(len(orbital)):
        # Joining the strings in the list for each orbital, while removing unnecessary white space.
        orbital[i] = ''.join([orbital[i][1:] for orbital[i] in orbital[i]])
    # Writing to file
    with open(file.replace('_Gaussian.txt',new_tail),'w') as output:
        for i in orbital:
            output.write(i+'\n')

