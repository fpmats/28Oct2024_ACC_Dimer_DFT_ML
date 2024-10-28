# Purpose of this script is similar to spherical_by_index.py
# however the centroid of the phenyl groups belonging to
# phenylalanines are not pre-defined, so we create them here.

# Importing relevant packages
import MDAnalysis as mda
import numpy as np
import time
import sys

start_time = time.time()

# Phenylalanine 1 as "{ChainID} {ResidueNumber}"
c1 = sys.argv[-4]
r1 = int(sys.argv[-3])
# Phenylalanine 2 as "{ChainID} {ResidueNumber}"
c2 = sys.argv[-2]
r2 = int(sys.argv[-1])

# Function to get the string for a variable name
def namestr(obj):
    return [name for name in globals() if globals()[name] is obj][0]

# Function to get the centroid of a phenyl group
def get_centroid(AtomGroup):
    positions = np.array(AtomGroup.positions, dtype='float64')
    centroid = np.mean(positions, axis=0)
    return centroid

# Function to define the vector that is normal
# to the plane containing a phenyl group
def get_norm2plane_vec(AtomGroup):
    positions = np.array(AtomGroup.positions, dtype='float64')
    v0 = positions[0] - positions[3]
    v1 = positions[1] - positions[5]
    v2 = positions[2] - positions[4]
    c0 = np.cross(v1,v0)/np.linalg.norm(np.cross(v1,v0))
    c1 = np.cross(-v0,v2)/np.linalg.norm(np.cross(-v0,v2))
    c2 = np.cross(v1,v2)/np.linalg.norm(np.cross(v1,v2))
    norm2plane_vec = np.mean([c0,c1,c2],axis=0)/np.linalg.norm(np.mean([c0,c1,c2],axis=0)) # Extra normalization was necessary for accuracy.
    return norm2plane_vec

# Function to get the dihedral angle between two phenyl groups
def get_dihedral_degrees(Ring1, Ring2):
    dihedral = np.arccos(np.dot(get_norm2plane_vec(Ring1),get_norm2plane_vec(Ring2)))*180/np.pi
    return dihedral

# Function to transform xyz coordinates to spherical coordinates
def asSpherical(xyz):
    #takes list xyz (single coord)
    x       = xyz[0]
    y       = xyz[1]
    z       = xyz[2]
    r       =  np.sqrt(x*x + y*y + z*z)
    theta   =  np.arccos(z/r)*180/ np.pi #to degrees
    phi     =  np.sign(y)*np.arctan2(y,x)*180/ np.pi
    return np.array([r,theta,phi])

# Function that returns the feature name
def feat_name(chainID_1, resid_1, chainID_2, resid_2):
    # Input is location info for two residues
    temp_universe =  mda.Universe('/projectnb/fpmats/common/all_1000_acc_dimer/snapshots/test_0.pdb')
    temp_atoms_1 = temp_universe.select_atoms(f'chainID {chainID_1} and resid {resid_1}')
    temp_atoms_2 = temp_universe.select_atoms(f'chainID {chainID_2} and resid {resid_2}')

    name = f'{temp_atoms_1[0].chainID}-{temp_atoms_1[0].resname}-{temp_atoms_1[0].resid}-{temp_atoms_2[0].chainID}-{temp_atoms_2[0].resname}-{temp_atoms_2[0].resid}'
    return name

# Function to return the vector between two phenyl group
# centroids and the dihedral created by them 
def centcent(n, chainID_1, resid_1, chainID_2, resid_2):
    if n>1000:
        universe = mda.Universe(f'/projectnb/fpmats/acc_dimer_share/All_hetero_solvation_RDF_val_100/{n-1000}_snapshot/test_{n-1000}.pdb')
    else:
        universe =  mda.Universe(f'/projectnb/fpmats/common/all_1000_acc_dimer/snapshots/test_{n}.pdb')

    ring1 = universe.select_atoms(f'chainID {chainID_1} and resid {resid_1} and type C and not (name CA or name CB or name C)')
    ring2 = universe.select_atoms(f'chainID {chainID_2} and resid {resid_2} and type C and not (name CA or name CB or name C)')
    return [*asSpherical(get_centroid(ring1)-get_centroid(ring2)),get_dihedral_degrees(ring1,ring2)]


filename = feat_name(c1, r1, c2, r2)
with open(filename, 'w') as f:
    f.write('# Output from FFcentroid.py and format is below'
            '\n# Snapshot Number, Distance (Angstroms), Phi (degrees), Theta (degrees), Dihedral (degrees)\n')
    f.write(f'#{filename}\n')
with open(filename,'a') as f:
    for n in range(0,1101):
        temp = centcent(n, c1, r1, c2, r2)
        f.write(f'{n:04d}    {temp[0]:7.3f}    {temp[1]:7.3f}    {temp[2]:7.3f}    {temp[3]:7.3f}\n')
with open(filename,'a') as f:
    f.write(f'# time elapsed: {time.time() - start_time}')
 
