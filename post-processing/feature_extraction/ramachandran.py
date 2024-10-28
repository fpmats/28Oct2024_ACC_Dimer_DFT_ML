# Purpose of this script is to obtain all the Ramachandran
# angles for each residue in the ACC-Dimer.

# Importing relevant python packages
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import time
import sys

start_time = time.time()

# Chain ID
c = sys.argv[-2]
# Residue ID
r = int(sys.argv[-1])

# Function to acquire a variable's name as a string
def namestr(obj):
    return [name for name in globals() if globals()[name] is obj][0]

# Function to create a feature name
def feat_name(chainID, resid):
    # Input is location info for a residue
    temp_universe =  mda.Universe('/path/to/some/PDB/file')
    temp_atoms = temp_universe.select_atoms(f'chainID {chainID} and resid {resid}')
    name = f'{temp_atoms[0].chainID}-{temp_atoms[0].resname}-{temp_atoms[0].resid}-Ramachandran-Angles'
    return name

# Function that will store the Ramachandran angles for a peptide
def ramarama(n, resid, chainID):
    # Input n: snapshot number (expects integer from 0 to 1100 inclusive)
    # Input resid: Residue ID as per the PDB file (expects integer from 1 to 29 inclusive)
    # Input chainID: Chain ID as per the PDB file (expects either 'A' or 'B')
    
    # Reading in the snapshot.
    if n>1000:
        universe = mda.Universe(f'/path/to/PDB/file/for/snapshots/numbered/over/1000/')
    else:
        universe =  mda.Universe(f'/path/to/PDB/file/for/snapshots/numbered/under/1000/')

    # Selecting residue from given inputs
    res = universe.select_atoms(f'chainID {chainID} and resid {resid}')

    # Acquiring Ramachandran angles
    x = Ramachandran(res).run()
    ramachandran_angles = x.results.angles

    # Returning time
    return ramachandran_angles[0][0]

filename=feat_name(c,r)

with open(filename,'w') as f:
    f.write('# Output from ramachandran.py format is below \n# Snapshot Number, Phi (degrees), Psi (degrees)\n')
    f.write(f'#{filename}\n')
with open(filename,'a') as f:
    for n in range(0,1101):
        temp = ramarama(n=n,resid=r,chainID=c)
        f.write(f'{n:04d}    {temp[0]:7.3f}    {temp[1]:7.3f}\n')
with open(filename,'a') as f:
    f.write(f'# time elapsed: {time.time() - start_time}')

