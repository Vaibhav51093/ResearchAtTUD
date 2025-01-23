from ovito.io import *
from ovito.modifiers import *
from ovito.pipeline import *
from ovito.data import *
import numpy as np

def trajectory(path_line='dump_nvt_prod.out', outfile='last_struct.data', frame=100):

    file = path_line

    # Data import 
    pipeline = import_file(file)

    # Setting up the atom type 
    pipeline.add_to_scene()
    def setup_particle_types(frame, data):
        types = data.particles_.particle_types_
        types.type_by_id_(1).name = "Zr"
        types.type_by_id_(2).name = "Si"
        types.type_by_id_(3).name = "P"
        types.type_by_id_(4).name = "Na"
        types.type_by_id_(5).name = "O"
    pipeline.modifiers.append(setup_particle_types)

    def modify1(frame: int, data: DataCollection):

        charge_dict = {"Si": 2.4, "Zr": 2.4, "P": 3.0, "O":-1.2, "Na":0.6}     # Padone charges 
        mass_dict = {"Zr": 91.224,
                        "Si": 28.086,
                        "P": 30.974,
                        "O": 15.999,
                        "Na": 22.99}
        charge = data.particles_.create_property("Charge", data = np.ones(data.particles.count))
        mass = data.particles_.create_property("Mass", data = np.ones(data.particles.count))
        ptypes = data.particles_.particle_types_
        for i in range(data.particles.count):
            charge[i] = charge_dict[ptypes.type_by_id(ptypes[i]).name]
            mass[i] = mass_dict[ptypes.type_by_id(ptypes[i]).name]
        
    pipeline.modifiers.append(modify1)

    num = pipeline.source.num_frames

    r = range(num)
    print(r)

    #for i in r:
    #print(i)
    #data = pipeline.compute(i)
    export_file(pipeline, outfile,format='xyz', multiple_frames = True, columns =
        ["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"], every_nth_frame=frame)
    
# Call the function

trajectory(path_line='/nfshome/deshmukh/vaibhav/schulze_project/vaibhav_sun/crystall/623/trajectory/structure.dump.*', outfile='trajectory.xyz', frame=1)
