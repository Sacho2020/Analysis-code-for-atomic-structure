# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import VoronoiAnalysisModifier, CommonNeighborAnalysisModifier, CreateBondsModifier

# Import NumPy module.
import numpy
import sys

# Load a simulation snapshot of a Cu-Zr metallic glass.
# pipeline = import_file("D:/INP_work/Analysis code/Spitial Correlation/Al78Cu.1069K.*")
pipeline = import_file("D:/INP_work/Analysis code/Spitial Correlation/al_pot.*.xyz")

# The LAMMPS dump file imported above contains only numeric atom type IDs but 
# no chemical element names or atom radius information. That's why we explicitly set the
# atomic radii of Cu & Zr atoms now (required for polydisperse Voronoi tessellation).
def assign_particle_radii(frame, data):
    atom_types = data.particles_.particle_types_
    atom_types.type_by_id_(1).radius = 1.43   # Al atomic radius assigned to atom type 1
    # atom_types.type_by_id_(2).radius = 1.28   # Cu atomic radius assigned to atom type 2
    atom_types.type_by_id_(2).radius = 1.6   # Zr atomic radius assigned to atom type 3
pipeline.modifiers.append(assign_particle_radii)

# # Create bonds based on the cutoff value
create_bonds_mod = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)

### Al-Cu-869K-cutoff
# create_bonds_mod.set_pairwise_cutoff("Al","Al",3.861)
# create_bonds_mod.set_pairwise_cutoff("Al","Cu",3.513)
# create_bonds_mod.set_pairwise_cutoff("Cu","Cu",3.455)

### Al-Zr-1821K-cutoff
# create_bonds_mod.set_pairwise_cutoff("Al","Al",3.492)
# create_bonds_mod.set_pairwise_cutoff("Al","Zr",3.917)
create_bonds_mod.set_pairwise_cutoff("Zr","Zr",4.342)

### Al-Cu-Zr-1634K-cutoff
# create_bonds_mod.set_pairwise_cutoff("Al","Al",3.809)
# create_bonds_mod.set_pairwise_cutoff("Al","Cu",3.632)
# create_bonds_mod.set_pairwise_cutoff("Al","Zr",4.104)
# create_bonds_mod.set_pairwise_cutoff("Cu","Cu",3.514)
# create_bonds_mod.set_pairwise_cutoff("Cu","Zr",3.986)
# create_bonds_mod.set_pairwise_cutoff("Zr","Zr",4.576)

pipeline.modifiers.append(create_bonds_mod)

# Set up the Voronoi analysis modifier.
voro = VoronoiAnalysisModifier(compute_indices = True, use_radii = True, edge_threshold = 0)
pipeline.modifiers.append(voro)

# Let OVITO compute the results.
for frame in range(pipeline.source.num_frames):
    data = pipeline.compute(frame)

    # Access computed Voronoi indices and the five-symmetry
    # This is an (N) x (M) array, where M is the maximum face order.
    voro_indices = data.particles['Voronoi Index'][:,2:6]
    f5 = []
    for i in range(len(voro_indices)):
        f5.append(voro_indices[i][2]/(voro_indices[i][0]+voro_indices[i][1]+voro_indices[i][2]+voro_indices[i][3]))
    
    # output the f5
    sys.stdout = open('D:/INP_work/Analysis code/Spitial Correlation/Atom_Five.'+str(frame+1)+'.dat','w')
    for particle_index in range(data.particles.count):
        if 8749 < particle_index:
        # # Print particle index (1-based).
            sys.stdout.write("%i\t%f " % (particle_index, f5[particle_index]))
            # End of particle line
            sys.stdout.write("\n")

    # output the average five-fold symmetry
    sys.stdout = open('D:/INP_work/Analysis code/Spitial Correlation/Average_Five_Symmetry.dat','a')
    # Print particle index (1-based).
    sys.stdout.write("%f" % (sum(f5)/(data.particles.count)))
    sys.stdout.write("\n")
    # output the bond lists
    Bond_torology = data.particles.bonds['Topology']
    sys.stdout = open('D:/INP_work/Analysis code/Spitial Correlation/Bond.'+str(frame+1)+'.dat','w')
    for i in range(len(Bond_torology)):
    # Print particle index (1-based).
        sys.stdout.write("%s" % (' '.join(str(j) for j in Bond_torology[i])))
    # End of particle line
        sys.stdout.write("\n")
