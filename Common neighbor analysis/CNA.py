# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import VoronoiAnalysisModifier, CommonNeighborAnalysisModifier, CreateBondsModifier
from ovito.data import BondsEnumerator

# Import standard Python and NumPy modules.
import numpy
import sys

####### Load the simulation dataset to be analyzed from AIMD simulation
# pipeline = import_file('D:/INP_work/RMC/Result/cfg-Al78Cu/869K/al_pot.3.xyz')
pipeline = import_file('D:/INP_work/RMC/Result/cfg-Al87Zr/1921K/al_pot.3.xyz')
# pipeline = import_file('D:/INP_work/RMC/Result/cfg-Al70Cu20Zr/1834K/al_pot.3.xyz')

# Create bonds based Voronoi method
# voro_bond = VoronoiAnalysisModifier(generate_bonds = True)
# pipeline.modifiers.append(voro_bond)

# # Create bonds based on the cutoff value
create_bonds_mod = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
#### Al-Cu
# create_bonds_mod.set_pairwise_cutoff("Al","Al",3.861)
# create_bonds_mod.set_pairwise_cutoff("Al","Cu",3.513)
# create_bonds_mod.set_pairwise_cutoff("Cu","Cu",3.455)
#### Al-Zr
create_bonds_mod.set_pairwise_cutoff("Al","Al",3.524)
create_bonds_mod.set_pairwise_cutoff("Al","Zr",3.806)
create_bonds_mod.set_pairwise_cutoff("Zr","Zr",4.476)
#### Al-Cu-Zr
# create_bonds_mod.set_pairwise_cutoff("Al","Al",3.843)
# create_bonds_mod.set_pairwise_cutoff("Al","Cu",3.664)
# create_bonds_mod.set_pairwise_cutoff("Al","Zr",4.081)
# create_bonds_mod.set_pairwise_cutoff("Cu","Cu",3.486)
# create_bonds_mod.set_pairwise_cutoff("Cu","Zr",4.022)
# create_bonds_mod.set_pairwise_cutoff("Zr","Zr",4.618)
pipeline.modifiers.append(create_bonds_mod)

# Compute CNA indices on the basis of the created bonds.
pipeline.modifiers.append(CommonNeighborAnalysisModifier(mode = CommonNeighborAnalysisModifier.Mode.BondBased))

def row_histogram(a):
    ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
    unique, indices, inverse = numpy.unique(ca, return_index=True, return_inverse=True)
    counts = numpy.bincount(inverse)
    sort_indices = numpy.argsort(counts)[::-1]
    return (a[indices[sort_indices]], counts[sort_indices])

# Let OVITO's data pipeline do the heavy work.
CNA_indices = numpy.array([], dtype= int)
for frame in range(pipeline.source.num_frames):
    data = pipeline.compute(frame)
    cna_indices = data.particles.bonds['CNA Indices']
    bond_enumerator = BondsEnumerator(data.particles.bonds)
    for particle_index in range(data.particles.count):
        #### Al-Cu alloy        
        # if particle_index < 7778:            ### Al-CNA
        # if 7777 < particle_index < 10000:      ### Cu-CNA 

        #### Al-Zr alloy
        # if particle_index < 8750:            ### Al-CNA
        if 8749 < particle_index < 10000:      ### Zr-CNA  
          
        #### Al-Cu-Zr alloy        
        # if particle_index < 7000:            ### Al-CNA
        # if 6999 < particle_index < 9000:      ### Cu-CNA
        # if 8999 < particle_index < 10000:      ### Zr-CNA

        ### all-CNA
        # if particle_index < 10000:              
            bond_index_list = list(bond_enumerator.bonds_of_particle(particle_index))
            local_cna_indices = cna_indices[bond_index_list]
            CNA_indices = numpy.append(CNA_indices, local_cna_indices)

# print(len(numpy.array(CNA_indices).reshape(-1,3)))
unique_indices, counts = row_histogram(CNA_indices.reshape(-1, 3))

# Print the ten most frequent histogram entries. 
for i in range(15):
    print("%s\t%i\t(%.1f %%)" % (tuple(unique_indices[i]),
                                counts[i],
                                100.0*float(counts[i])/len(numpy.array(CNA_indices).reshape(-1,3))))

