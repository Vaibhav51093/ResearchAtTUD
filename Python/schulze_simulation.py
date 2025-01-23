lmp_input = """\
# This script will run the calculations in one run in pyiron.
# official work of Vaibhav Arun Deshmukh, TU-Darmstadt

#------------------------------Cook and Quench method---------------------------#

#---------------------------- Atomic setup ------------------------------------#
units metal
dimension 3
boundary p p p
atom_style full
read_data structure.inp
include potential.inp
neigh_modify     delay 0
#------------------------------------------------------------------------------#


# cook and quencxh 

fix		1 all npt temp 1. 4200. 0.02 iso 0. 0. 2.        # 
run		2100000                                          #new aim 1K / ps  ##100000 now it is 500 ps
unfix 		1

fix		1 all nvt temp 4200. 4200. 0.02
run 		100000
unfix 1

restart	1000000 restart.*

fix		1 all npt temp 4200. 1. 0.02 iso 0. 0. 2.
run		2099500

write_data      npt-4200K_cooldown_1Kps.data
write_restart	npt-4200K_cooldown_1Kps.restart
"""