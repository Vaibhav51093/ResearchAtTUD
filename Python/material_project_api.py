import itertools as int 
from pyiron import Project
import numpy as np
import pandas as pd 
from mp_api import MPRester
import os
import requests
from pymatgen.io.cif import CifWriter


# Make sure that, you have the Materials project API key 
MP_API_KEY="2pmOKuKG8dXsfqzFS8WO9s1udVz9hQ1n"

# Data from API can be accessed using REST request 
response = requests.get("https://api.materialsproject.org/materials/mp-149/?fields=material_id%2Cstructure&all_fields=false", 
                                headers={"X-API-Key": MP_API_KEY})  
#print(response.text)
# We would like to interact to the MATERIALS PROJECT    
mpr = MPRester(MP_API_KEY)


class get_info():# How to get info using chemical elements from materials project 

    def get_mat_info_materials_project(element, hull_cut=0.00, cif_file=False):

        mat_id = []
        mat_compo = []
        mat_form = []
        mat_for = []
        mat_above_h = []
        name_cif = []
    
        for i in element:
            with MPRester(MP_API_KEY) as mpr:
                results = mpr.query(chemsys="%s"%i)
                rang = range(len(results))
                for i in rang:
                    doc = results[i]
                    hull = doc.energy_above_hull
                    if hull == 0 or hull <= hull_cut:
                        if cif_file == False:
                            material_id = doc.material_id
                            formation = doc.formation_energy_per_atom
                            compo = doc.composition_reduced
                            form = doc.composition.reduced_formula
                            mat_id.append(material_id)
                            mat_compo.append(compo)
                            mat_for.append(formation)
                            mat_above_h.append(hull)
                            mat_form.append(form)
                        else:
                            material_id = doc.material_id
                            structure = doc.structure
                            formation = doc.formation_energy_per_atom
                            compo = doc.composition_reduced
                            form = doc.composition.reduced_formula
                            mat_id.append(material_id)
                            mat_compo.append(compo)
                            mat_for.append(formation)
                            mat_above_h.append(hull)
                            cif = CifWriter(structure)
                            name = material_id+".cif"
                            name_cif.append(name)
                            cif.write_file(name)
                            mat_form.append(form)

                    
        return mat_id, mat_compo, mat_for, mat_above_h, mat_form, name_cif