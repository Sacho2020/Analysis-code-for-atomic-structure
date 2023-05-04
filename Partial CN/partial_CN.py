
import numpy as np
from rich import print
from rich.progress import track
from ovito.io import import_file
from ovito.data import BondsEnumerator

# Load a dataset containing atoms and bonds.
pipeline = import_file("E:\\python_tool for structure analysis\\Calculate_partial_CN\\LAMMPS_data_file", atom_style='bond')

# Obtain pipeline results.
data = pipeline.compute()
bond_topology = data.particles.bonds.topology  # array with bond topology
ptypes = data.particles.particle_types
bonds_enum = BondsEnumerator(data.particles.bonds)

#Loop over atoms.
out = np.zeros((data.particles.count, 5), dtype = int)
for particle_index in track(range(data.particles.count)):
    out[particle_index][0] = particle_index
    out[particle_index][1] = ptypes[particle_index]
    atom_type = []
    for bond_index in bonds_enum.bonds_of_particle(particle_index):
        num = (np.where(bond_topology[bond_index] != particle_index))
        atom_type.append(ptypes[bond_topology[bond_index][num]])
        out[particle_index][2] = atom_type.count(1)
        out[particle_index][3] = atom_type.count(2)
        out[particle_index][4] = atom_type.count(3)

np.savetxt("E:\\python_tool for structure analysis\\Calculate_partial_CN\\partial_CN.txt", out, fmt = "%d" )

# Generate an empty list for each type of atom
atype = []
for i in range(3):
    atype.append([])
for i in range(np.array(out).shape[0]):
    if out[i][1] == 1:
        atype[0].append(out[i][2:])
    if out[i][1] == 2:
        atype[1].append(out[i][2:])
    if out[i][1] == 3:
        atype[2].append(out[i][2:])

# calculate the partial coordination number for each type of atom 
for i in range(len(atype)):
    print("原子类型{}的偏配位数分别为（1, 2, 3）： ".format(i+1))
    print(np.array(atype[i]).mean(axis = 0))




# #### code for OVITO’s Python programming interface
# from ovito.data import *
# from ovito.data import BondsEnumerator
# import numpy as np

# def modify(frame, data):
#     bond_topology = data.particles.bonds.topology # array with bond topology
#     ptypes = data.particles.particle_types
#     bonds_enum = BondsEnumerator(data.particles.bonds)
#     out = np.zeros((data.particles.count, 5), dtype = int)
#     for particle_index in range(data.particles.count):
#         out[particle_index][0] = particle_index
#         out[particle_index][1] = ptypes[particle_index]
#         atom_type = []
#         for bond_index in bonds_enum.bonds_of_particle(particle_index):
#             num = (np.where(bond_topology[bond_index] != particle_index))
#             atom_type.append(ptypes[bond_topology[bond_index][num]])
#             out[particle_index][2] = atom_type.count(1)
#             out[particle_index][3] = atom_type.count(2)
#             out[particle_index][4] = atom_type.count(3)
#     atype = []
#     for i in range(3):
#         atype.append([])
#     for i in range(np.array(out).shape[0]):
#         if out[i][1] == 1:
#             atype[0].append(out[i][2:])
#         if out[i][1] == 2:
#             atype[1].append(out[i][2:])
#         if out[i][1] == 3:
#             atype[2].append(out[i][2:])
#     for i in range(len(atype)):
#         print("原子类型{}的偏配位数分别为（1, 2, 3）： ".format(i+1))
#         print(np.array(atype[i]).mean(axis = 0))