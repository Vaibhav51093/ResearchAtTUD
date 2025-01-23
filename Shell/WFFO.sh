#!/bin/bash

# Define the grid range and step size
min_val=0
max_val=100
step_size=1

# Loop through the grid
for i in $(seq $min_val $step_size $max_val); do
  for j in $(seq $min_val $step_size $max_val); do
    # Construct the file name for the current point in the grid
    file_name="structure_${i}_${j}.data"
    
    # Replace the file name inside the input file
    sed "s/IDp02_gb-structure-ID18\.data/$file_name/g" pedone_tabular_minimize_template.in > input_file.in
    
    # Run the simulation for the current file
    echo "Running simulation for $file_name"
    
    # Construct the log file name for the current point in the grid
    log_file_name="log_file_${i}_${j}.log"
    
    # Replace the command below with the command to run your simulation using LAMMPS
    lammps -in input_file.in -log $log_file_name #< $file_name > output_file.out
  done
done

