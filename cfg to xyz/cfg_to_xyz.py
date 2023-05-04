################################################################################  
#FILENAME:cfg_to_xyz.py
#Author:Cao Saichao
#DESCRIPTION:a batch processingtool written by python language for converting .cfg to .xyz format files for RMC simulation
#Notes:
################################################################################  

import numpy as np
import os
import linecache

filepath = 'E:\\RMC_processing\\cfg_file\\'

pathdir = os.listdir(filepath)

p = []
pathdir_new = []
for filename in pathdir:
    if os.path.splitext(filename)[1] == ".cfg":
        p.append(np.loadtxt(filepath+filename, skiprows=30))  ###三元合金
        # p.append(np.loadtxt(filepath+filename, skiprows=26))   ###二元合金
        pathdir_new.append(filename)

coordinate = []
vector = float(np.array(linecache.getline(filepath+pathdir_new[0], 15).strip().split(" "))[0])
for item in p:
    a = (item*vector).tolist()
    coordinate.append(a)

#### 三元合金
atom_number = []
for item in [19, 23, 27]:
    number = np.array(linecache.getline(filepath+pathdir_new[0], item).strip().split(" "))[0]
    atom_number.append(int(number))

#### 二元合金
# atom_number = []
# for item in [19, 23]:
#     number = np.array(linecache.getline(filepath+pathdir_new[0], item).strip().split(" "))[0]
#     atom_number.append(int(number))

atom_type = int(input("请输入原子种类数量: "))

name = []
if atom_type == 1:
    atom_name_1 = input("请输入元素符号(如Cr): ")
    for i in range(atom_number[0]):
        name.append(atom_name_1)
elif atom_type == 2:
    atom_name_1 = input("请输入元素符号(如Cr): ")
    atom_name_2 = input("请输入元素符号(如Cr): ")
    for i in range(atom_number[0]):        
        name.append(atom_name_1)
    for i in range(atom_number[1]):            
        name.append(atom_name_2)
elif atom_type == 3:
    atom_name_1 = input("请输入元素符号(如Cr): ")
    atom_name_2 = input("请输入元素符号(如Cr): ")
    atom_name_3 = input("请输入元素符号(如Cr): ")
    for i in range(atom_number[0]):        
        name.append(atom_name_1)
    for i in range(atom_number[1]):            
        name.append(atom_name_2)
    for i in range(atom_number[2]):        
        name.append(atom_name_3)

for k in range(len(pathdir_new)):
    filename = pathdir_new[k]
    file = open (filepath+os.path.splitext(filename)[0]+".xyz", "w")
    file.write ("{}\n" .format(len(name)))
    file.write ("\n")
    for i in range(len(name)):
            file.write ("{0}\t{1:.7f}\t{2:.7f}\t{3:.7f}\n".format(name[i], coordinate[k][i][0], coordinate[k][i][1], coordinate[k][i][2]))
    file.close()