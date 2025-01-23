#!/bin/bash

nohup python msd_crystal_fra_2.py &
wait $!

nohup python msd_crystal_fra_3.py &
wait $!

nohup python msd_crystal_fra_4.py &
wait $!

nohup python msd_crystal_fra_5.py &
wait $!

#nohup python msd_amorph_fra_5.py &
#wait $!

