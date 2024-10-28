# The purpose of this script is to take the output from Clean.py
# and compile it into a singular JSON file that can be used for
# other purposes.

# Importing things
import os
from glob import iglob
import MDAnalysis as mda
import numpy as np
import json

# Function to get the index for the HOMO
def get_homo_number(orbitals):
    states = [int(x[2]) for x in orbitals]
    HOMO = states[-1]
    # One-liner for HOMO: HOMO = [int(x[2]) for x in orbitals][-1]
    return HOMO

# Function to get the number of near-gap states
def get_number_near_gap(orbitals, OEs, cutoff):
    HOMO = [int(x[2]) for x in orbitals][-1]
    # Gets the orbital number corresponding to HOMO
    x = [y*27.2114 for y in OEs[get_homo_number(orbitals)-len(orbitals):get_homo_number(orbitals)] if np.abs(y*27.2114-OEs[get_homo_number(orbitals)-1]*27.2114)<=cutoff]
    # x is the list of energies for each occupied orbital in orbitals within cutoff eV (up until HOMO)
    my_list = len(x)
    return my_list

# Function to get the orbital character (e.g. s, p, d)
def get_character(orbital):
    character = [x[x.find('-')+1] for x in orbital[5:]]
    # get_character(orbital)[i] is the orbital character for entry i in an orbital.
    return character

def num_contributors(orbital):
    x = len(orbital[5:])
    return x

def get_contributions(orbital):
    contributions = [float(x[x.find('=')+1:]) for x in orbital[5:]]
    # get_contributions(orbital)[i] is the mulliken charge for a 
    return contributions
  
def get_indices(orbital):
    indices = [int(x.split('-')[0][1:]) for x in orbital[5:]] # Index = Gaussian indexing. Index-1 as below for MDAnalysis.
    return indices
  
def get_locations(orbital):
    indices = [int(x.split('-')[0][1:])-1 for x in orbital[5:]] # Matches atom with Gaussian, but needs to be -1 for MDAnalysis indexing.
    name_chain_res_resid = [ [system[i].name, system[i].chainID, system[i].resname, int(system[i].resid), system[i].type] for i in indices]
    # get_location(orbital)[0] = PDB Name for Atom
    # get_location(orbital)[1] = PDB Name for Peptide Chain
    # get_location(orbital)[2] = Three-letter name for Peptide Residue
    # get_location(orbital)[3] = Integer Index for peptide residue. Starts from 1 for the first entry in PDB file.
    # int()" was added because the previously returned value was numpy.int64 type and that could not be written to JSON.
    # get_location(orbital)[4] = Element. 
    return name_chain_res_resid

master = {}

# Below are the loops used to append to the master dictionary.
# Note there are two loops due to snapshots being split 1001 and 100.
for n in range(0,1001):
    occupied_20 = f'/path/to/current/working/directory/3_post-processing/json_making/{n}_Gaussian_cleaned.txt'
    occ_20 = [f for f in iglob(occupied_20, recursive=True) if os.path.isfile(f)]
    pdb_dir = f'/path/to/parent/dir/for/each/PDB/file/'
    pdb_file = [f for f in iglob(pdb_dir, recursive=True) if os.path.isfile(f)]
    OE_dir = f'/path/to/orbital_energies.out/file/'
    OE_file = [f for f in iglob(OE_dir, recursive=True) if os.path.isfile(f)]
    atoms = mda.Universe(pdb_file[0])
    system = atoms.select_atoms('not element Cl K')
    my_bigger_dictionary = {}
    with open(occ_20[0],'r') as Clean_Gauss:
        orbitals = [line.split() for line in Clean_Gauss]
    with open(OE_file[0],'r') as so_many_energies:
        all_of_them = [line.split() for line in so_many_energies]
    OEs = [float(OE) for row in all_of_them[1:-1] for OE in row]
    my_dictionary = {}
    for j in orbitals:
        my_smaller_dictionary = {}
        assert len(get_character(j))==num_contributors(j),"Previous definition for number of contributors does not equal new definition"
        for i in range(num_contributors(j)):
            my_smaller_dictionary[str(i)] = {
                'atom': get_indices(j)[i],
                'location': get_locations(j)[i],
                'character': get_character(j)[i],
                'mulliken': get_contributions(j)[i],
            }
            my_dictionary[f'HOMO-{get_homo_number(orbitals)-int(j[2])}'] = {
                'energy': OEs[int(j[2])-1]*27.2114,
                'composition': my_smaller_dictionary,
            }
    my_bigger_dictionary[n]=f'Snapshot {n}'
    my_bigger_dictionary['occ_near_gap']= get_number_near_gap(orbitals,OEs,0.2)
    my_bigger_dictionary['details']=my_dictionary
    master[n] = my_bigger_dictionary


for n in range(1,101):
    occupied_20 = f'/path/to/current/working/directory/3_post-processing/json_making/{n}_Gaussian_cleaned.txt'
    occ_20 = [f for f in iglob(occupied_20, recursive=True) if os.path.isfile(f)]
    pdb_dir = f'/path/to/parent/dir/for/each/PDB/file/'
    pdb_file = [f for f in iglob(pdb_dir, recursive=True) if os.path.isfile(f)]
    OE_dir = f'/path/to/orbital_energies.out/file/'
    OE_file = [f for f in iglob(OE_dir, recursive=True) if os.path.isfile(f)]
    atoms = mda.Universe(pdb_file[0])
    system = atoms.select_atoms('not element Cl K')
    my_bigger_dictionary = {}
    with open(occ_20[0],'r') as Clean_Gauss:
        orbitals = [line.split() for line in Clean_Gauss]
    with open(OE_file[0],'r') as so_many_energies:
        all_of_them = [line.split() for line in so_many_energies]
    OEs = [float(OE) for row in all_of_them[1:-1] for OE in row]
    my_dictionary = {}
    for j in orbitals:
        my_smaller_dictionary = {}
        assert len(get_character(j))==num_contributors(j),"Previous definition for number of contributors does not equal new definition"
        for i in range(num_contributors(j)):
            my_smaller_dictionary[str(i)] = {
                'atom': get_indices(j)[i],
                'location': get_locations(j)[i],
                'character': get_character(j)[i],
                'mulliken': get_contributions(j)[i],
            }
            my_dictionary[f'HOMO-{get_homo_number(orbitals)-int(j[2])}'] = { #THIS IS J FOR THE HOMO SHENANIGANS
                'energy': OEs[int(j[2])-1]*27.2114,
                'composition': my_smaller_dictionary,
            }
    my_bigger_dictionary[n+1000]=f'Snapshot {n+1000}'
    my_bigger_dictionary['occ_near_gap']= get_number_near_gap(orbitals,OEs,0.2)
    my_bigger_dictionary['details']=my_dictionary
    master[n+1000] = my_bigger_dictionary

with open('jason.json','w') as output:
    json.dump(master,output)

