# The purpose of this script is to read the Gaussian log files
# and then capture the near-gap orbitals reported.

#!/bin/bash
for (( j=0; j<=1000; j++ ))

do
    cat /path/to/gaussian/log/files/*/${j}_snapshot/*.1/*.1-g16.log | grep "Atomic contributions to Alpha molecular orbitals:" -A 10000 | grep "Alpha vir" -B 10000 -m 1 > ${j}_Gaussian.txt
done


n=1000
for (( j=1; j<=100; j++ ))

do
    cat /path/to/gaussian/log/files/*/${j}_snapshot/*.1/*.1-g16.log | grep "Atomic contributions to Alpha molecular orbitals:" -A 10000 | grep "Alpha vir" -B 10000 -m 1 > $((${n}+${j}))_Gaussian.txt
done

