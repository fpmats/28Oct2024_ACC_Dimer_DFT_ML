# Purpose of this script is to take specified hydrogen bonds and
# then sort them so that we keep track of the bonds by length rather
# than by atom index

# Importing pandas and numpy
import numpy as np
import pandas as pd

# File list of features for the specific hydrogen bonds of interest
file_list='A-GLU-1-HT1_A-GLU-1-OE1 A-GLU-1-HT1_A-GLU-1-OE2 A-GLU-1-HT2_A-GLU-1-OE1 A-GLU-1-HT2_A-GLU-1-OE2 A-GLU-1-HT3_A-GLU-1-OE1 A-GLU-1-HT3_A-GLU-1-OE2 B-GLU-1-HT1_B-GLU-1-OE1 B-GLU-1-HT1_B-GLU-1-OE2 B-GLU-1-HT2_B-GLU-1-OE1 B-GLU-1-HT2_B-GLU-1-OE2 B-GLU-1-HT3_B-GLU-1-OE1 B-GLU-1-HT3_B-GLU-1-OE2 B-LYS-14-HZ1_B-GLU-15-OE1 B-LYS-14-HZ1_B-GLU-15-OE2 B-LYS-14-HZ2_B-GLU-15-OE1 B-LYS-14-HZ2_B-GLU-15-OE2 B-LYS-14-HZ3_B-GLU-15-OE1 B-LYS-14-HZ3_B-GLU-15-OE2'

file_list=file_list.split(' ')

# Electrostatics in A1, B1, and between B14 and B15
A1 , B1, B1415 = np.split(np.array(file_list),3)

# A1
data = []
for i in A1:
    data.append(np.loadtxt(i,usecols=1,skiprows=2))

for n in range(1,7):
    with open(f'A-GLU-1-HTX_A-GLU-1-OEX_{n}','w') as f:
        f.write('# Output from sorter.py format is below\n# Snapshot Number, Distance (Angstroms)\n')
        f.write(f'#A-GLU-1-HTX_A-GLU-1-OEX_{n}\n')
        for snapshot,distance in enumerate(np.sort(data)[n-1]):
            f.write(f'{snapshot:04d}    {distance:7.3f}\n')
# B1
data = []
for i in B1:
    data.append(np.loadtxt(i,usecols=1,skiprows=2))

for n in range(1,7):
    with open(f'B-GLU-1-HTX_B-GLU-1-OEX_{n}','w') as f:
        f.write('# Output from sorter.py format is below\n# Snapshot Number, Distance (Angstroms)\n')
        f.write(f'#B-GLU-1-HTX_B-GLU-1-OEX_{n}\n')
        for snapshot,distance in enumerate(np.sort(data)[n-1]):
            f.write(f'{snapshot:04d}    {distance:7.3f}\n')

# B14 with B15
data = []
for i in B1415:
    data.append(np.loadtxt(i,usecols=1,skiprows=2))

for n in range(1,7):
    with open(f'B-LYS-14-HZX_B-GLU-15-OEX_{n}','w') as f:
        f.write('# Output from sorter.py format is below\n# Snapshot Number, Distance (Angstroms)\n')
        f.write(f'#B-LYS-14-HZX_B-GLU-15-OEX_{n}\n')
        for snapshot,distance in enumerate(np.sort(data)[n-1]):
            f.write(f'{snapshot:04d}    {distance:7.3f}\n')

