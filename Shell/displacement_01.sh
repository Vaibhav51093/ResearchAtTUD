#!/bin/bash -l

echo 0
cd 0
cp POSCAR POSCAR_init
cp CONTCAR POSCAR
phonopy -d --dim="2 2 2"
#cd ..

for NUM in {1..5} #must be consistent with volumes.sh/directory names
do
        echo $NUM-
        cd $NUM-
        mv POSCAR POSCAR_init
        mv CONTCAR POSCAR
        phonopy -d --dim="2 2 2"
        cd ..
done

for NUM in {1..5} #must be consistent with volumes.sh/directory names
do
        echo $NUM+
        cd $NUM+
        mv POSCAR POSCAR_init
        mv CONTCAR POSCAR
        phonopy -d --dim="2 2 2"
        cd ..
done
