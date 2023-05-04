# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import VoronoiAnalysisModifier
from ovito.data import BondsEnumerator

# Import standard NumPy modules.
import numpy

# Load the simulation dataset obtained from RMC simulation to be analyzed
# pipeline = import_file('D:/INP_work/RMC/Result/cfg-Al78Cu/869K/al_pot.869K.xyz')
pipeline = import_file('D:/INP_work/RMC/Result/cfg-Al87Zr/1821K/al_pot.1821K.xyz')
# pipeline = import_file('D:/INP_work/RMC/Result/cfg-Al70Cu20Zr/1834K/al_pot.3.xyz')

# atomic radii of  atoms now (required for polydisperse Voronoi tessellation).
def assign_particle_radii(frame, data):
    atom_types = data.particles_.particle_types_
    ### Al-Zr
    atom_types.type_by_id_(1).radius = 1.43   # atomic radius assigned to atom type (Al)
    atom_types.type_by_id_(2).radius = 1.60   # atomic radius assigned to atom type (Zr)

pipeline.modifiers.append(assign_particle_radii)

# Set up the Voronoi analysis modifier.
voro = VoronoiAnalysisModifier(compute_indices = True, use_radii = True, edge_threshold = 0)
pipeline.modifiers.append(voro)

# Let OVITO compute the results.
data = pipeline.compute()

# Access computed Voronoi indices.
# This is an (N) x (M) array, where M is the maximum face order.
voro_indices = data.particles['Voronoi Index']

def row_histogram(a):
    ca = numpy.ascontiguousarray(a).view([('', a.dtype)] * a.shape[1])
    unique, indices, inverse = numpy.unique(ca, return_index=True, return_inverse=True)
    counts = numpy.bincount(inverse)
    sort_indices = numpy.argsort(counts)[::-1]
    return (a[indices[sort_indices]], counts[sort_indices])

# Let OVITO's data pipeline do the heavy work.
Voro_indices = numpy.array([], dtype= int)
F5 = numpy.array([], dtype= float)
for frame in range(pipeline.source.num_frames):
    data = pipeline.compute(frame)
    voro_indices = data.particles['Voronoi Index']
    for particle_index in range(data.particles.count):

        #### Al-Zr alloy
        # if particle_index < 8750:            ### Al-CNA
        # if 8749 < particle_index < 10000:      ### Zr-CNA  
          
        #### all
        if particle_index < 10000:
            local_voro_indices_n10 = voro_indices[particle_index]
            local_voro_indices = local_voro_indices_n10[2:6] #<n3 n4 n5 n6>
            local_f5 = local_voro_indices[2]/(local_voro_indices[0]+local_voro_indices[1]+local_voro_indices[2]+local_voro_indices[3])
            Voro_indices = numpy.append(Voro_indices, local_voro_indices)
            F5 = numpy.append(F5, local_f5)

# Compute frequency histogram and average five-fold symmetry.
unique_indices, counts = row_histogram(Voro_indices.reshape(-1,4))
Average_F5 = numpy.average(F5)
print("平均五重对称性:" + str(Average_F5))

### Compute the fraction of VP with five-fold symmetry larger than 0.5.
n = 0
for i in range(len(F5)):
    if F5[i] > 0.5:
        n += 1
Bigger_than_half = n/len(Voro_indices.reshape(-1,4))
print("fraction of f5 > 0.5:" +  str(Bigger_than_half))

### number of VP types
print("团簇类型的数量是：" + str(len(unique_indices)))

### varience of cluster fraction
counts_fraction = counts/len(Voro_indices.reshape(-1,4))*100
print("团簇分数的变量是：" + str(numpy.var(counts_fraction)))

# # Print the ten most frequent histogram entries. 
for i in range(10):
    print("%s\t%i\t(%.2f %%)" % (tuple(unique_indices[i]),
                                counts[i],
                                100.0*float(counts[i])/len(Voro_indices.reshape(-1,4))))