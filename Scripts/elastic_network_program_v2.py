#!/usr/bin/env python


import pandas as pd
import math
df = pd.DataFrame(columns=["Bead_no","Bead_name","Residue","Chain_id","Residue_num","x","y","z",
                           "mass","charge","isOverhang"])


Bead_no = []
Bead_name = []
Residue = []
Chain_id = []
Residue_num = []
x = []
y = []
z = []
mass = []
charge = []
path="PATHNAME/"+"CG.pdb"
file = open("PATHNAME/"+"spine_defn.txt","w")
with open(path) as f:
    for line in f.read().split("\n"):
        if(line.split()[0] == "TER"):
            break
        Bead_no.append(line.split()[1])
        Bead_name.append(line.split()[2])
        Residue.append(line.split()[3])
        Chain_id.append(line.split()[4])
        Residue_num.append(line.split()[5])
        x.append(line.split()[6])
        y.append(line.split()[7])
        z.append(line.split()[8])
        mass.append(line.split()[9])
        charge.append(line.split()[10])
    df["Bead_no"] = Bead_no
    df["Bead_name"] = Bead_name
    df["Residue"] = Residue
    df["Chain_id"] = Chain_id
    df["Residue_num"] = Residue_num
    df["x"] = x
    df["y"] = y
    df["z"] = z
    df["charge"] = charge
    df["mass"] = mass
    df["isOverhang"] = "No" 
    df["isHairpinLoop"] = "No"



#Spine definitions
print(";with the spine and the molecule" +"\n")
print("#define RUBBER_FC_SPINE_constraint 800.0 \n")
file.write(";with the spine and the molecule \n")
file.write("#define RUBBER_FC_SPINE_constraint 800.0 \n")
spine_coords = []
molecular_coords = []
for index,row in df.iterrows():
    if(row["Bead_name"]=="CENT"):
        spine_coords.append({"atom_no":row["Bead_no"],"x":row["x"],"y":row["y"],"z":row["z"]})
    else:
        molecular_coords.append({"atom_no":row["Bead_no"],"x":row["x"],"y":row["y"],"z":row["z"]})

#radial distance to try
max_radial_dis = 1.8
min_radial_dis = 0.2
# elastic network definitions
# print(centroid)
for centroid in spine_coords:
    for atom in molecular_coords:
        dis = distance((centroid['atom_no'],centroid['x'],centroid['y'],centroid['z']),(atom['atom_no'],atom['x'],atom['y'],atom['z']))
        dis = dis/10.0
        if(dis<max_radial_dis):
             print(centroid["atom_no"]+" "+atom["atom_no"]+" "+"6"+" "+str(round((dis),2))+" "+"RUBBER_FC_SPINE_constraint*1.000000"+" \n")
             file.write(centroid["atom_no"]+" "+atom["atom_no"]+" "+"6"+" "+str(round((dis),2))+" "+"RUBBER_FC_SPINE_constraint*1.000000"+" \n\n")

file.close()

left_chain=[1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
right_chain=[40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21]

left_chain=[9, 10, 11, 12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
right_chain=[40,39,38,37,36,35,34,33,32,31,30,8,7,6,5,4,3,2,1]

residues_strand_C = []
residues_strand_D = []

print("#define RUBBER_FC 300.000000")
print(";HBond elastic network")

for i in left_chain:
    for index,row in df.iterrows():
        if(int(df['Residue_num'][index]) == i):
            if(row['Bead_name'] == "ASC2" or row['Bead_name'] == 'ASC3' 
               or row['Bead_name'] == "USC2" or row['Bead_name'] == 'USC3' 
              or row['Bead_name'] == "GSC2" or row['Bead_name'] == "GSC3"
              or row['Bead_name'] == "CSC2" or row['Bead_name'] == "CSC3"):
                residues_strand_C.append((df['Bead_name'][index],df['Bead_no'][index],
                                      float(df['x'][index]),float(df['y'][index]),float(df['z'][index])))
for j in right_chain:
    for index,row in df.iterrows():
        if(int(df['Residue_num'][index]) == j):
            if(row['Bead_name'] == "ASC2" or row['Bead_name'] == 'ASC3' 
               or row['Bead_name'] == "USC2" or row['Bead_name'] == 'USC3' 
              or row['Bead_name'] == "GSC2" or row['Bead_name'] == "GSC3"
              or row['Bead_name'] == "CSC2" or row['Bead_name'] == "CSC3"):
                residues_strand_D.append((df['Bead_name'][index],df['Bead_no'][index],float(df['x'][index]),
                              float(df['y'][index]),float(df['z'][index])))

for k in range(0,len(residues_strand_C)):
    dis = distance((residues_strand_C[k][1],residues_strand_C[k][2],residues_strand_C[k][3],
                    residues_strand_C[k][4]),
                  (residues_strand_D[k][1],residues_strand_D[k][2],residues_strand_D[k][3],
                   residues_strand_D[k][4]))
    dis = dis/10.0
    print(residues_strand_C[k][1]+" "+residues_strand_D[k][1]+" "+"6"+" "+str(round(dis,2))+" "
          +"RUBBER_FC*1.000000"+" "+'\n')


# In[47]:


len(residues_strand_D)


# #### dsRNA old H-bond implementation

# In[8]:


## HBond definitions
residues_strand_C = []
residues_strand_D = []
for index,row in df.iterrows():
    if(row["Chain_id"] == "C" and row["isOverhang"] == "No"):
        if(row["Bead_name"] == "ASC2"):
            residues_strand_C.append(("ASC2",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "ASC3"):
            residues_strand_C.append(("ASC3",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "USC2"):
            residues_strand_C.append(("USC2",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "USC3"):
            residues_strand_C.append(("USC3",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "GSC2"):
            residues_strand_C.append(("GSC2",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "GSC3"):
            residues_strand_C.append(("GSC3",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "CSC2"):
            residues_strand_C.append(("CSC2",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "CSC3"):
            residues_strand_C.append(("CSC3",row["Bead_no"],row["x"],row["y"],row["z"]))
    
    elif(row["Chain_id"] == "D" and row["isOverhang"] == "No"):
        if(row["Bead_name"] == "ASC2"):
            residues_strand_D.append(("ASC2",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "ASC3"):
            residues_strand_D.append(("ASC3",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "USC2"):
            residues_strand_D.append(("USC2",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "USC3"):
            residues_strand_D.append(("USC3",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "GSC2"):
            residues_strand_D.append(("GSC2",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "GSC3"):
            residues_strand_D.append(("GSC3",row["Bead_no"],row["x"],row["y"],row["z"]))
        elif(row["Bead_name"] == "CSC2"):
            residues_strand_D.append(("CSC2",row["Bead_no"],row["x"],row["y"],row["z"])) 
        elif(row["Bead_name"] == "CSC3"):
            residues_strand_D.append(("CSC3",row["Bead_no"],row["x"],row["y"],row["z"]))        


# In[9]:


residues_strand_C_invert = []
for i in range(0,len(residues_strand_C),2):
    residues_strand_C_invert.append(residues_strand_C[i+1])
    residues_strand_C_invert.append(residues_strand_C[i])


# In[10]:


i = len(residues_strand_C)-1
j = 0
print("#define RUBBER_FC 300.000000")
print(";HBond elastic network")
for k in range(0,len(residues_strand_C)):
#     print(residues_strand_C_invert[i][0],residues_strand_D[j][0])
    dis = distance((residues_strand_C_invert[i][1],residues_strand_C_invert[i][2],residues_strand_C_invert[i][3],
                   residues_strand_C_invert[i][4]),(residues_strand_D[j][1],residues_strand_D[j][2],
                   residues_strand_D[j][3],residues_strand_D[j][4]))
    dis = dis/10.0
    print(residues_strand_C_invert[i][1]+" "+residues_strand_D[j][1]+" "+"6"+" "+str(round(dis,2))+" "
          +"RUBBER_FC*1.000000"+" "+'\n')
    i = i-1
    j = j+1


# In[6]:


def distance(point_A,point_B):
#     print(point_A[0],point_B[0])
    x1 = float(point_A[1])
    y1 = float(point_A[2])
    z1 = float(point_A[3])
    x2 = float(point_B[1])
    y2 = float(point_B[2])
    z2 = float(point_B[3])
    dis = math.sqrt(math.pow((x1-x2),2) + math.pow((y1-y2),2) + math.pow((z1-z2),2) )
    return dis


# In[12]:


#grooves height code
groove_beads_C_BB1 = []
groove_beads_D_BB1 = []

groove_beads_C_BB2 = []
groove_beads_D_BB2 = []

groove_beads_C_BB3 = []
groove_beads_D_BB3 = []

for index,row in df.iterrows():
    if(row['Bead_name'] == "BB1" and row['Chain_id'] == 'B'):
        groove_beads_C_BB1.append((row["Bead_no"],row["x"],row["y"],row["z"]))
#     if(row['Bead_name'] == "BB1" and row['Chain_id'] == 'D'):
#         groove_beads_D_BB1.append((row["Bead_no"],row["x"],row["y"],row["z"])) 
        
    if(row['Bead_name'] == "BB2" and row['Chain_id'] == 'B'):
        groove_beads_C_BB2.append((row["Bead_no"],row["x"],row["y"],row["z"]))
#     if(row['Bead_name'] == "BB2" and row['Chain_id'] == 'D'):
#         groove_beads_D_BB2.append((row["Bead_no"],row["x"],row["y"],row["z"])) 

    if(row['Bead_name'] == "BB3" and row['Chain_id'] == 'B'):
        groove_beads_C_BB3.append((row["Bead_no"],row["x"],row["y"],row["z"]))

#groove width code
print("#define RUBBER_FC 300.0")
print(";Groove width elastic network")
residues_strand_C_bb = []
residues_strand_D_bb = []
for index,row in df.iterrows():
    if(row["Chain_id"] == "C" and row["isOverhang"] == "No"):    
        if(row["Bead_name"] == "BB1"):
            residues_strand_C_bb.append(("BB1",row["Bead_no"],row["x"],row["y"],row["z"]))
        if(row["Bead_name"] == "BB2"):
            residues_strand_C_bb.append(("BB2",row["Bead_no"],row["x"],row["y"],row["z"]))
#         if(row["Bead_name"] == "BB3"):
#             residues_strand_C_bb.append(("BB3",row["Bead_no"],row["x"],row["y"],row["z"]))
        ####################################################################
    elif(row["Chain_id"] == "D" and row["isOverhang"] == "No"):
        if(row["Bead_name"] == "BB1"):
            residues_strand_D_bb.append(("BB1",row["Bead_no"],row["x"],row["y"],row["z"]))
        if(row["Bead_name"] == "BB2"):
            residues_strand_D_bb.append(("BB2",row["Bead_no"],row["x"],row["y"],row["z"]))
#         if(row["Bead_name"] == "BB3"):
#             residues_strand_D_bb.append(("BB3",row["Bead_no"],row["x"],row["y"],row["z"]))
        #######################################################################
residues_strand_C_invert_bb = []
for i in range(0,len(residues_strand_C_bb),2):
#     residues_strand_C_invert_bb.append(residues_strand_C_bb[i+2])
    residues_strand_C_invert_bb.append(residues_strand_C_bb[i+1])
    residues_strand_C_invert_bb.append(residues_strand_C_bb[i])

i = len(residues_strand_C_invert_bb)-1
j = 0
for k in range(0,len(residues_strand_C_invert_bb)):
    dis = distance((residues_strand_C_invert_bb[i][1],residues_strand_C_invert_bb[i][2],residues_strand_C_invert_bb[i][3],
                   residues_strand_C_invert_bb[i][4]),(residues_strand_D_bb[j][1],residues_strand_D_bb[j][2],
                   residues_strand_D_bb[j][3],residues_strand_D_bb[j][4]))
    print(residues_strand_C_invert_bb[i][1]+" "+residues_strand_D_bb[j][1]+" "+"6"+" "+
          str(round((dis/10.0),2))+" "
          +"RUBBER_FC_width*1.000000"+" \n")
    i = i-1
    j = j+1

