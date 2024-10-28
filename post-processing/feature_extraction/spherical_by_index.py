# The purpose of this script is to take two atom indices
# and then generate the vector connecting the two in
# spherical coordinates

# Importing relevant packages
import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import distances
import time
import sys

start_time = time.time()

# Residue Numbers
indices = [int(sys.argv[-2]), int(sys.argv[-1])]

# Function to get a variable's name as a string
def namestr(obj):
    return [name for name in globals() if globals()[name] is obj][0]

# Function to convert Cartesian coordinates to spherical
def asSpherical(xyz):
    #takes list xyz (single coord)
    x       = xyz[0]
    y       = xyz[1]
    z       = xyz[2]
    r       =  np.sqrt(x*x + y*y + z*z)
    theta   =  np.arccos(z/r)*180/ np.pi #to degrees
    phi     =  np.sign(y)*np.arctan2(y,x)*180/ np.pi
    return np.array([r,theta,phi])

# Function that prints the atoms selected
def locate(index1, index2):
    # Input are two atom indices that will be located
    temp_universe =  mda.Universe('/path/to/PDB/file')
    temp_atoms = temp_universe.select_atoms(f'index {index1} {index2}')
    # Prints atom information by "column"
    print(f'Atoms selected:\n',
      f'Chains {temp_atoms.chainIDs}\n',
      f'Residues {temp_atoms.resnames}\n',
      f'Atoms {temp_atoms.names}\n',
      f'Indices {temp_atoms.indices}\n',)

# Function to get a feature name
def feat_name(index1, index2):
    # Input are two atom indices that will be used for the name
    temp_universe =  mda.Universe('/path/to/PDB/file')
    temp_atoms1 = temp_universe.select_atoms(f'index {index1}')
    temp_atoms2 = temp_universe.select_atoms(f'index {index2}')
    name = f'{temp_atoms1[0].chainID}-{temp_atoms1[0].resname}-{temp_atoms1[0].resid}-{temp_atoms1[0].name}_{temp_atoms2[0].chainID}-{temp_atoms2[0].resname}-{temp_atoms2[0].resid}-{temp_atoms2[0].name}'
    return name

# Function that returns the spherical coordinates
def dists(n, index1, index2):
    # Input n: snapshot number (expects integer from 0 to 1100 inclusive)
    # Input index1 index2: atom indeces of interest (expects integer greater or equal to 0)

    # Reading in the snapshot.
    if n>1000:
        universe = mda.Universe(f'/path/to/PDB/file/for/a/snapshot/number/above/1000')
    else:
        universe =  mda.Universe(f'/path/to/PDB/file/for/a/snapshot/number/under/1000')

    # Selecting atoms from given indices
    i1 = universe.select_atoms(f'index {index1}')
    i2 = universe.select_atoms(f'index {index2}')

    # Acquiring displacement vector from index1 to index2
    vector = i2[0].position-i1[0].position

    #Converting to spherical coordinates using function above
    return asSpherical(vector)


locate(indices[0],indices[1])

filename = feat_name(indices[0], indices[1])

with open(filename,'w') as f:
    f.write('# Output from spherical_by_index.py format is below \n# Snapshot Number, Distance (Angstroms), Phi (degrees), Theta (degrees)\n')
    f.write(f'#{feat_name(indices[0],indices[1])}\n')
with open(filename,'a') as f:
    for n in range(0,1101):
        temp = dists(n,indices[0],indices[1])
        f.write(f'{n:04d}    {temp[0]:7.3f}    {temp[1]:7.3f}   {temp[2]:7.3f}\n')
with open(filename,'a') as f:
    f.write(f'# time elapsed: {time.time() - start_time}')
