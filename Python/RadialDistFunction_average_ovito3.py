##############
# General usage:
# ovitos path/to/file > RDF.log
#
# Usually:
# ovitos ~/scripts/Ovito_analysis/RadialDistFunction_average_ovito3.py > RDF.log
#
# Generated with   ovito-3.0.0-dev628-x86_64/bin/ovitos
#
##############


from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier
from ovito.modifiers import SelectTypeModifier
from ovito.modifiers import InvertSelectionModifier
from ovito.modifiers import DeleteSelectedModifier
import numpy

### Parameters ########################################################################################################################################
in_filename       = 'XDATCAR_S1_S2_renamed'
Prefix_filename   = 'RDF' # There will be scaled RDFs (all partial sum up to the total RDF) and unscaled (as if only those particles were in the box)
EquilibrationTime = 500   # Skip this number of STEPS from the XDATCAR. Not necessary fs -> remember NBLOCK from INCAR! 
cut               = 7.0   # cutoff
bins              = 700   # number of bins
AtomTypes         = ['Li','P','S2','S1','Br']  # same order as in XDATCAR!
######################################################################################################################################################## 


print("\nStart RDF analysis of:  {}\n".format(in_filename))
print("With atom types:")
for i in AtomTypes:
    print(i)

print("\n############\nParameters are:\n")
print("Skipped steps:      {}".format(EquilibrationTime))
print("Cutoff for RDF:     {}".format(cut))
print("Number of bins:     {}".format(bins))
print("############\n")


### Load a simulation trajectory consisting of several frames:
pipeline = import_file(in_filename)



### Extract the atom concentrations ( <Number of X atoms> / <Total number of atoms> )
concentrations    = []

sel_modifier = SelectTypeModifier(
        operate_on = "particles",
        property = "Particle Type")
pipeline.modifiers.append(sel_modifier)
pipeline.modifiers.append(InvertSelectionModifier())

sel_modifier.enabled = False
data = pipeline.compute()
Total_Atom_Number = data.particles.count

pipeline.modifiers.append(DeleteSelectedModifier())

sel_modifier.enabled = True
for i in AtomTypes:
    sel_modifier.types = {i}
    data = pipeline.compute()
    concentrations.append(data.particles.count/Total_Atom_Number)
sel_modifier.enabled = False



### Insert the modifier into the pipeline:
pipeline = import_file(in_filename)
modifier = CoordinationAnalysisModifier(cutoff = cut, number_of_bins = bins, partial=True)
pipeline.modifiers.append(modifier)



### Initialize array for accumulated RDF:
ArraySize     = int(   len(AtomTypes)*(len(AtomTypes)+1) / 2   + 1   )
partial_rdf   = numpy.zeros((modifier.number_of_bins, ArraySize))



### Iterate over all frames of the sequence and calculate RDF:
for frame in range(EquilibrationTime,pipeline.source.num_frames):
    if frame % 1000 == 0:
        print("Computing frame %d of %d" % (frame,pipeline.source.num_frames),flush=True)
    # Evaluate pipeline to let the modifier compute the RDF of the current frame:
    data = pipeline.compute(frame)
    # Accumulate RDF histograms:
    partial_rdf += data.tables['coordination-rdf'].as_table()



# Averaging the RDFs over number of frames:
frame_count = pipeline.source.num_frames-EquilibrationTime
partial_rdf /= frame_count



# Set up header for output file and the right factors that need to be multiplied with the rdf values.
header  = '#Distance  '
factors = numpy.array([1])  # Initialize with a 1 as factor for the distance
for i in range(len(AtomTypes)):
    for j in range(i,len(AtomTypes)):
        AtomCombination= AtomTypes[i] + '-' + AtomTypes[j] 
        header = header + AtomCombination.center(11)
        if i == j:
            factor = concentrations[i]**2
            factors = numpy.append(factors, factor)
        else:
            factor = 2*concentrations[i]*concentrations[j]
            factors = numpy.append(factors, factor)

header = header + '\n'

scaled_partial_rdf = partial_rdf * factors


# Export the average unscaled and scaled RDF to text files:
out_filename = Prefix_filename + '_unscaled.txt'
with open(out_filename, 'w') as f:
    f.write(header)
    numpy.savetxt(f, partial_rdf,fmt='%.4e')

out_filename = Prefix_filename + '_scaled.txt'
with open(out_filename, 'w') as f:
    f.write(header)
    numpy.savetxt(f, scaled_partial_rdf,fmt='%.4e')


