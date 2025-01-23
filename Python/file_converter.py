#!/nfshome/deshmukh/miniconda3/bin/python

from ase import io
import sys

#infile = sys.argv[1]
#outfile = sys.argv[2]
#conv_format = sys.argv[2]

atoms = io.read('%s'%sys.argv[1])
atoms.write('%s'%sys.argv[2], format = '%s'%sys.argv[3])
