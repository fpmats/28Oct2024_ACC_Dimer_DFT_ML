# Each script and its intended purpose.

## `buncher.py`
 Purpose of this script is to take the giant list of features
 and splitting them into types of interactions that are
 examined.

## `combiner.py`
 Purpose of this script is to combine all the outputs from
 `spherical_by_index.py` and the later supplied SASAs csv file
 into one giant features csv file.

## `FFcentroid.py`
 Purpose of this script is similar to `spherical_by_index.py`
 however the centroid of the phenyl groups belonging to
 phenylalanines are not pre-defined, so we create them here.

## `ramachandran.py`
 Purpose of this script is to obtain all the Ramachandran
 angles for each residue in the ACC-Dimer.


## `sorter.py`
 Purpose of this script is to take specified hydrogen bonds and
 then sort them so that we keep track of the bonds by length rather
 than by atom index

## `spherical_by_index.py`
 The purpose of this script is to take two atom indices
 and then generate the vector connecting the two in
 spherical coordinates

