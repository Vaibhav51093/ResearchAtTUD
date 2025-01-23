import os
import re



def species_and_nelect_from_POTCAR(pot_dir):
    with open(os.path.join(pot_dir, "POTCAR"), "r") as f:
        first_lines = [next(f) for l in range(24)] 
        nelect = float(first_lines[1].strip())
        ele_line = first_lines[3]
        if "VRHFIN" in ele_line:
            ele = re.split("=|:", ele_line)[1]
        return ele, nelect

def update_vasp_potential_list(LDA_dir="potpaw", PBE_dir="potpaw_PBE"): 
    LDA_potentials = [d for d in os.listdir(LDA_dir) if os.path.isdir(os.path.join(LDA_dir, d))]
    PBE_potentials = [d for d in os.listdir(PBE_dir) if os.path.isdir(os.path.join(PBE_dir, d))]  
    
    LDA_potentials.sort()
    PBE_potentials.sort()  
 
    with open("potentials_vasp.csv", "w") as f:
        f.write(",Filename,Model,Name,Species,n_elect\n")
        index=0
        for pot in PBE_potentials:
            element, nelect = species_and_nelect_from_POTCAR(os.path.join(PBE_dir, pot))
            f.write("{},potpaw_PBE/{}/POTCAR,gga-pbe,{}-gga-pbe,['{}'],{}\n".format(
                index,
                pot,
                pot,
                element,
                nelect,
                ))
            index += 1
        for pot in LDA_potentials:
            element, nelect = species_and_nelect_from_POTCAR(os.path.join(LDA_dir, pot))
            f.write("{},potpaw/{}/POTCAR,lda,{}-lda,['{}'],{}\n".format(
                index,
                pot,
                pot,
                element,
                nelect,
                ))
            index += 1
        

if __name__ == "__main__":
    update_vasp_potential_list()    


