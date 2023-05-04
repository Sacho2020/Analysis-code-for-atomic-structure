import numpy as np
from Spitial_Correlation import SpitialCorrelation, draw_heatmap

### calculate the correlation for all configurations
cij = []
for i in range(1,4):
    cij.append(SpitialCorrelation(i))
    print()
cij_average = np.mean(cij,axis=0)
print(np.around((cij_average),3))

### draw the heatmap
# xlabels=['0','0.5','1']
# ylabels=['0','0.5','1']
# draw_heatmap(cij_average,xlabels,ylabels)

# ### draw the heatmap
xlabels=['0','0.2','0.4','0.6','0.8','0.1']
ylabels=['0','0.2','0.4','0.6','0.8','0.1']
draw_heatmap(cij_average,xlabels,ylabels)

# ### draw the heatmap
# xlabels=['0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1']
# ylabels=['0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1']
# draw_heatmap(cij_average,xlabels,ylabels)