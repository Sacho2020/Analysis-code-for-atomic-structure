import numpy as np

def SpitialCorrelation(Num_Config):

    ## inputdata----atom number and f5
    filepath='D:/INP_work/Analysis code/Spitial Correlation'
    Atom_Five = np.loadtxt(filepath+'/Atom_Five.'+str(Num_Config)+'.dat')
    f5 = Atom_Five[:,1]
    Atom_Num = Atom_Five[:,0]
    pos = Atom_Five[:,:]
    ## inputdata----bond pair
    Bond = np.loadtxt(filepath+'/Bond.'+str(Num_Config)+'.dat')
    Bond_pairs = Bond[:,:].tolist()

    #################################### 根据f5进行分类
    ##创建含有10个子列表的列表（10组不同的结构）
    Group = []
    Group_Atom_Num= []
    for i in range(11):
        Group.append([])
        Group_Atom_Num.append([])
    ##分3类
    # for i in range(len(f5)):
    #     if 0 <= f5[i] < 0.5:
    #         Group[0].append(Atom_Num[i])
    #     if f5[i] == 0.5:
    #         Group[1].append(Atom_Num[i])
    #     if 0.5 < f5[i] <= 1:
    #         Group[2].append(Atom_Num[i])
    # ##分5类
    # for i in range(len(f5)):
    #     if f5[i] == 0:
    #         Group[0].append(Atom_Num[i])
    #     if 0 < f5[i] <= 0.2:
    #         Group[1].append(Atom_Num[i])
    #     if 0.2 < f5[i] <= 0.4:
    #         Group[2].append(Atom_Num[i])
    #     if 0.4 < f5[i] <= 0.6:
    #         Group[3].append(Atom_Num[i])
    #     if 0.6 < f5[i] <= 0.8:
    #         Group[4].append(Atom_Num[i])
    #     if 0.8 < f5[i] <= 1:
    #         Group[5].append(Atom_Num[i])
    ##分10类
    for i in range(len(f5)):
        if f5[i] == 0:
            Group[0].append(Atom_Num[i])
        if 0 < f5[i] <= 0.1:
            Group[1].append(Atom_Num[i])
        if 0.1 < f5[i] <= 0.2:
            Group[2].append(Atom_Num[i])
        if 0.2 < f5[i] <= 0.3:
            Group[3].append(Atom_Num[i])
        if 0.3 < f5[i] <= 0.4:
            Group[4].append(Atom_Num[i])
        if 0.4 < f5[i] <= 0.5:
            Group[5].append(Atom_Num[i])
        if 0.5 < f5[i] <= 0.6:
            Group[6].append(Atom_Num[i])
        if 0.6 < f5[i] <= 0.7:
            Group[7].append(Atom_Num[i])
        if 0.7 < f5[i] <= 0.8:
            Group[8].append(Atom_Num[i])
        if 0.8 < f5[i] <= 0.9:
            Group[9].append(Atom_Num[i])
        if 0.9 < f5[i] <= 1:
            Group[10].append(Atom_Num[i])

    ## 计算每组中的原子个数
    for i in range(len(Group)):
        Group_Atom_Num[i] = len(Group[i])
    print(Group_Atom_Num)
    
    def Bond_Number(a,b,c):
        G = set([tuple(t) for t in c])
        # print(len(G))
        bond_num = 0
        h = []

        if a == b:
            for i in a:
                for j in b:
                    if (i,j) in G:
                        h.append([i,j])
                        bond_num +=1
        else:
            for i in a:
                for j in b:
                    if (i,j) in G or (j,i) in G:
                        h.append([i,j])
                        bond_num +=1
        return(bond_num)
    #################################### 调用Bond_Number函数计算各组之间的成键数量
    ## 创建10*10列表
    Pij_real = []
    Pij_theo = []
    Cij = []
    for i in range(len(Group_Atom_Num)):
        Pij_real.append([])
        Pij_theo.append([])
        Cij.append([])
        for j in range(len(Group_Atom_Num)):
            Pij_real[i].append([])
            Pij_theo[i].append([])
            Cij[i].append([])
    ## 计算对应位置的键的数量及计算成键的真实概率Pij
    AllBond_number = int(len(Bond))
    for i in range(len(Group_Atom_Num)):
        for j in range(len(Group_Atom_Num)):
            Pij_real[i][j] = Bond_Number(Group[i],Group[j],Bond_pairs)/AllBond_number

    ###################################计算成键的理论值Pij0
    N = int(len(Atom_Five))
    for i in range(len(Group_Atom_Num)):
        for j in range(len(Group_Atom_Num)):
            if i == j:
                Pij_theo[i][j] = Group_Atom_Num[i]*(Group_Atom_Num[i]-1)/(N*(N-1))
            else:
                Pij_theo[i][j] = 2*Group_Atom_Num[i]*Group_Atom_Num[j]/(N*(N-1))            

    #####################################计算相关性Cij
    for i in range(len(Group_Atom_Num)):
        for j in range(len(Group_Atom_Num)):
            Cij[i][j] = Pij_real[i][j]/Pij_theo[i][j]-1
    
    for i in range(len(Group_Atom_Num)):
        print(Cij[i])

    return(Cij)

### draw the heatmap
from matplotlib import pyplot as plt
from matplotlib import cm 
from matplotlib import axes
def draw_heatmap(data,xlabels,ylabels):
    cmap = cm.get_cmap('rainbow',1000)    
    figure=plt.figure(facecolor='w')
    ax=figure.add_subplot(1,1,1,position=[0.1,0.15,0.8,0.8])
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels)
    vmax=data[0][0]
    vmin=data[0][0]
    for i in data:
        for j in i:
            if j>vmax:
                vmax=j
            if j<vmin:
                vmin=j
    map=ax.imshow(data,interpolation='nearest',cmap=cmap,aspect='auto',vmin=vmin,vmax=vmax)
    plt.colorbar(mappable=map,cax=None,ax=None,shrink=1)
    plt.show()

