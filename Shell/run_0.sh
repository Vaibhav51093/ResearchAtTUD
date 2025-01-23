#!/bin/bash

# Define the grid range and step size
min_val=0
max_val=99
step_size=1

# Make a copy of the input file with placeholders
cp control.inp input_file.in

# Loop through the grid
for i in $(seq $min_val $step_size $max_val); do
#for j in $(seq $min_val $step_size $max_val); do
  # Construct the file name for the current point in the grid
  file_name="structure_${i}_${i}.data"
  
  cp control.inp input_file.in

  # Replace the file name inside the input file
  sed -i "s/replace/$file_name/g" input_file.in
    
  # Run the simulation for the current file
  echo "Running simulation for $file_name"
    
  # Construct the log file name for the current point in the grid
  log_file_name="log_file_${i}_${i}.log"
    
  # Replace the command below with the command to run your simulation using LAMMPS
  lmp -i input_file.in -log $log_file_name #< $file_name > output_file.out
done
#done

