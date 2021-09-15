from gemmi import cif
import pandas as pd
import numpy as np
import json

from pathlib import Path

DATA_PATH = str(Path(__file__).absolute().parent) + '/data'

from .atom_radius import atom_radius

with open(DATA_PATH+'/atomic_radii.json',mode='r') as handle:
    atomic_radii_lib = json.load(handle)

def secondlayer(filename):
    """
    Port of the 'secondlayer.f90' program in the published m-CBAC code
    This code determines atomic connectivity based on the m-CBAC algorithm
    and returns the two-layer connectivity relationships
    
    Inputs:
    f = filename of P 1 CIF
       
    Outputs:
    connectivity = a Pandas dataframe containing the connectivity layers
       column '0': the zero layer (just the constituent atoms)
       column '1': the first layer
       column '2': the second layer
    """
    cif_data = cif.read(filename).sole_block()
    # Unit Cell Parameters
    a = cif_data.find_values('_cell_length_a')
    b = cif_data.find_values('_cell_length_b')
    c = cif_data.find_values('_cell_length_c')
    box = np.array([float(x[0]) for x in [a,b,c]])
    alpha = cif_data.find_values('_cell_angle_alpha')
    beta = cif_data.find_values('_cell_angle_beta')
    gamma = cif_data.find_values('_cell_angle_gamma')
    angles = np.array([float(x[0]) for x in [alpha,beta,gamma]])
    volume = box[0] * box[1] * box[2] * np.sqrt(1. - np.cos(angles[0]*np.pi/180.)**2 
                                    - np.cos(angles[1]*np.pi/180.)**2
                                    - np.cos(angles[2]*np.pi/180.)**2
                                    + 2. * np.cos(angles[0]*np.pi/180.) 
                                    * np.cos(angles[1]*np.pi/180.) * np.cos(angles[2]*np.pi/180.)
                                )
    
    # Atom Symbols [e.g., Carbon=C, Oxygen=O, etc.]
    try:
        atom_symbols = [x for x in cif_data.find_loop('_atom_site_type_symbol')]
    except:
        atom_symbols = []
    # Atom Labels [i.e., the specific type of the atom]
    try:
        atom_labels = [x for x in cif_data.find_loop('_atom_site_label')]
    except:
        atom_labels = []
    
    # If necessary, convert atom_labels to atom_symbols
    if len(atom_symbols) == 0:
        atom_symbols = []
        for atom in atom_labels:
            atom_symbols.append(''.join(i for i in atom if not i.isdigit()))
        
    # Atom Positions
    fracx = np.array(cif_data.find_loop('_atom_site_fract_x'), dtype=float)
    fracy = np.array(cif_data.find_loop('_atom_site_fract_y'), dtype=float)
    fracz = np.array(cif_data.find_loop('_atom_site_fract_z'), dtype=float)
    Natoms = len(atom_symbols)
    
    #-----------------------------------------------
    # Massive rewrite possible here
    
    # Lattice Vectors
    a1 = box[0]
    a2 = 0.
    a3 = 0.
    b1 = box[1] * np.cos(angles[2]*np.pi/180.)
    b2 = box[1] * np.sin(angles[2]*np.pi/180.)
    b3 = 0.
    c1 = box[2] * np.cos(angles[1]*np.pi/180.)
    c2 = box[2] * (np.cos(angles[0]*np.pi/180.)-np.cos(angles[1]*np.pi/180.)*np.cos(angles[2]*np.pi/180.))/np.sin(angles[2]*np.pi/180.)
    c3 = volume / (box[0] * box[1] * np.sin(angles[2]*np.pi/180.))
      
    # Build an expanded 3x3x3 cell and a 1x1x1 primitive cell 
    Natoms_exp = Natoms * 27  # this is done to avoid PBC calculations
    atom_symbols_exp = []
    ratoms_exp = []
    ratoms = []
    for xrep in range(-1,2):
        for yrep in range(-1,2):
            for zrep in range(-1,2):
                for index in range(Natoms):
                    x = ((fracx[index] + float(xrep))*a1 + (fracy[index] + float(yrep))*b1
                        + (fracz[index] + float(zrep))*c1)
                    y = ((fracx[index] + float(xrep))*a2 + (fracy[index] + float(yrep))*b2
                        + (fracz[index] + float(zrep))*c2)
                    z = ((fracx[index] + float(xrep))*a3 + (fracy[index] + float(yrep))*b3
                        + (fracz[index] + float(zrep))*c3)
                    atom_symbols_exp.append(atom_symbols[index])
                    ratoms_exp.append([x,y,z])
                    if xrep == 0 and yrep == 0 and zrep == 0:
                        ratoms.append([x,y,z])
    #-----------------------------------------------


    # Calculate distances       
    rank_list = []
    num_list = []
    for index0 in range(Natoms):
        distance = [100.]*8
        rank = [None]*8
        num = [0]*8
        for index1 in range(Natoms_exp):
            dr = [ abs( ratoms[index0][idim] - ratoms_exp[index1][idim] ) for idim in range(3) ]
            dist = np.sqrt(sum([x**2 for x in dr]))
            Rcut = (atom_radius(atom_symbols[index0],atomic_radii_lib)
                    + atom_radius(atom_symbols_exp[index1],atomic_radii_lib))*1.25
           
            for rank_index in range(0,8,+1):
                if dist < distance[rank_index] and dist < Rcut:
                    for i in range(7,rank_index,-1):
                        distance[i] = distance[i-1]
                        rank[i] = rank[i-1]
                        num[i] = num[i-1]
                    distance[rank_index] = dist
                    rank[rank_index] = atom_symbols_exp[index1]
                    
                    flo = np.floor((index1+1)/Natoms)
                    # NOTE: +1 here compared to LCL because python starts at zero
                    if int(flo)*Natoms == (index1+1):
                        num[rank_index] = Natoms
                    else:
                        num[rank_index] = (index1+1) - int(flo)*Natoms
                    break
                            
#        print(rank,num)
        if index0%100 == 0:
            print('finished ', index0)
        rank_list.append(rank)
        num_list.append(num)
    
    # Rank the connectivity based on the m-CBAC algorithm
    list_2nd = []
    list_1st = []
    list_0th = []
    for index0 in range(Natoms):
        m = sum([1 for x in num_list[index0] if x != 0])

        # Data for "2nd.txt"
        line = [[ rank_list[num_list[index0][i]-1][j] for j in range(1,8)] for i in range(1,m)]
        new_line = []
        for part in line:
            new_line += [x for x in part if x is not None]
        new_line = ''.join(sorted(new_line))
        list_2nd.append(new_line)
            
        # Data for "1st.txt"
        line = [ rank_list[index0][i] for i in range(1,m)]
        list_1st.append(''.join(sorted(line)))
    
        # Data for "0th.txt"
        list_0th.append(rank_list[index0][0])

#         # Check for floating atoms
#         l = 0
#         for i in range(1,m):
#             n = sum([1 for x in num_list[num_list[index0][i]] if x != 0])
#             l = l + n -1
#             if m > l+1:
#                 print("float")
        
    print('finished determining connectivity')
    connectivity = pd.DataFrame(list(zip(list_0th,list_1st,list_2nd)), columns=['0', '1', '2'])
    
    return connectivity


