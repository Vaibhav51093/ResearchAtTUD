import sys
import glob
import os 

path = "/home/schulze/simulations/manual-structures/jochen_18/00-minimize/01-heat-523/02-diffusion/"
tot_list = glob.glob(f"{path}/structure.dump.*")
list = [file for file in tot_list if int(file.split(".dump.")[-1])%1000 == 0]  #Only every 10th dump-file loaded --> 500 of 5000 dump files 	

n_frames = len(list)
half_frame = n_frames//2
every_nth_frame = 1


for i in range(every_nth_frame, half_frame, every_nth_frame):
    #Make submit script:
    script = f"#!/bin/bash \n #PBS -q big \n #PBS -l nodes=1:ppn=1 \n #PBS -l walltime=00:30:00 \n #PBS -N id18-diff-smooth \n date \n # Change to working dir.\n cd $PBS_O_WORKDIR \n /home/schulze/software/ovito-pro-3.4.1-x86_64/bin/ovitos Quismas-Calculation-Diffusivity-smooth-auto.py {path} {i} {1000} \n date "
    #save string script to file
    with open(f"submit_smooth{i}.sh", "w") as text_file:
        text_file.write(script)
    os.system(f"qsub submit_smooth{i}.sh")
    print("OK")
    
