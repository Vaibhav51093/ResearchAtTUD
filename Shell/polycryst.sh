#!/bin/bash

# Already we have crystal structure unit cell 

# Remove old files 
rm *cfg *txt *xsf

# Write text file 

echo "box 250 250 0" >> poly.txt
echo "random 50" >> poly.txt

# Generate polycrsytal
atomsk --polycrystal nasicon_ordered_lisette_unit_cell_mini.cif poly.txt poly.cfg  -wrap
