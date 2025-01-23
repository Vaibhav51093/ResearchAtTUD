# Boilerplate code generated by OVITO Pro 3.7.2
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *

# Data import:
pipeline = import_file('/nfshome/deshmukh/vaibhav/project/bulk_NaSICON_MD/dump/dump_amo_723.out', multiple_frames = True)

# Manual modifications of the imported data objects:
def modify_pipeline_input(frame: int, data: DataCollection):
    data.particles_.particle_types_.type_by_id_(1).color = (0.6666666865348816, 0.0, 1.0)
    data.particles_.particle_types_.type_by_id_(1).radius = 0.2
    data.particles_.particle_types_.type_by_id_(2).color = (1.0, 1.0, 0.0)
    data.particles_.particle_types_.type_by_id_(2).radius = 0.5
    data.particles_.particle_types_.type_by_id_(3).radius = 0.2
    data.particles_.particle_types_.type_by_id_(4).radius = 0.2
    data.particles_.particle_types_.type_by_id_(5).radius = 0.2
pipeline.modifiers.append(modify_pipeline_input)

# Wrap at periodic boundaries:
pipeline.modifiers.append(WrapPeriodicImagesModifier())

# Create bonds:
mod = CreateBondsModifier()
mod.mode = CreateBondsModifier.Mode.Pairwise
pipeline.modifiers.append(mod)
mod.set_pairwise_cutoff(3, 1, 2.5)

# Select type:
pipeline.modifiers.append(SelectTypeModifier(types = {2, 4, 5}))

# Delete selected:
pipeline.modifiers.append(DeleteSelectedModifier())

# Cluster analysis:
pipeline.modifiers.append(ClusterAnalysisModifier(
    neighbor_mode = ClusterAnalysisModifier.NeighborMode.Bonding, 
    sort_by_size = True))

# Expression selection:
pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'Cluster >= 10'))

# Delete selected:
pipeline.modifiers.append(DeleteSelectedModifier())

# Color coding:
pipeline.modifiers.append(ColorCodingModifier(
    property = 'Cluster', 
    start_value = 1.0, 
    end_value = 9.0, 
    gradient = ColorCodingModifier.Viridis()))

# Select type:
pipeline.modifiers.append(SelectTypeModifier(types = {3}))

# Coordination polyhedra:
pipeline.modifiers.append(CoordinationPolyhedraModifier())