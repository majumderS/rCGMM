#!/usr/bin/env python
# coding: utf-8


#import modules
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
import pandas as pd
import math


###### ------------------------------------------------ Function definitions --------------------------------------------

def find_neighbour_density(point_A,set_of_points):
    neighbours = []
    for point_B in set_of_points:
        if(point_A != point_B):
            dis = distance(point_A,point_B)
            if(dis < 7):
                neighbours.append(set_of_points[i])
    return {"Point":point_A,"count":len(neighbours)}


def distance(point_A,point_B):
    x1 = float(point_A[0])
    y1 = float(point_A[1])
    z1 = float(point_A[2])
    x2 = float(point_B[0])
    y2 = float(point_B[1])
    z2 = float(point_B[2])
    dis = math.sqrt(math.pow((x1-x2),2) + math.pow((y1-y2),2) + math.pow((z1-z2),2) )
    return dis


def centroid(points):
    sum_x = 0.0
    sum_y = 0.0
    sum_z = 0.0
    l = float(len(points))
    for point in points:
        sum_x = sum_x + point[0]
        sum_y = sum_y + point[1]
        sum_z = sum_z + point[2]
    return [round((sum_x/l),3),round((sum_y/l),3),round((sum_z/l),3)]


def return_grid_coordinates(df):
    x_max = df['x'][0]
    y_max = df['y'][0]
    z_max = df['z'][0]

    x_min = df['x'][0]
    y_min = df['y'][0]
    z_min = df['z'][0]

    for value in df['x']:
        x_coordinate = float(value)
        if (x_coordinate > x_max):
            x_max = x_coordinate
    for value in df['y']:
        y_coordinate = float(value)
        if (y_coordinate > y_max):
            y_max = y_coordinate
    for value in df['z']:
        z_coordinate = float(value)
        if (z_coordinate > z_max):
            z_max = z_coordinate
    for value in df['x']:
        x_coordinate = float(value)
        if (x_coordinate < x_min):
            x_min = x_coordinate
    for value in df['y']:
        y_coordinate = float(value)
        if (y_coordinate < y_min):
            y_min = y_coordinate
    for value in df['z']:
        z_coordinate = float(value)
        if (z_coordinate < z_min):
            z_min = z_coordinate
    return ([x_min,x_max,y_min,y_max,z_min,z_max])

def get_midpoint(A,B):
    mdpt = []
    i = (A[0]+B[0])/2
    j = (A[1]+B[1])/2
    k = (A[2]+B[2])/2
    mdpt.append(round(i,2))
    mdpt.append(round(j,2))
    mdpt.append(round(k,2))
    return mdpt


#parse PDB
df = pd.DataFrame(columns=["Bead_no","Bead_name","Residue","Chain_id","Residue_num","x","y","z",
                           "mass","charge","isOverhang"])
Bead_no = []
Bead_name = []
Residue = []
Chain_id = []
charge = []
mass = []
Residue_num = []
x = []
y = []
z = []
path="/home/subhasree/IISC_Bangalore/My_PhD_Studies/Lab_work/oralcancer_mirna_3d/pdb/miR1273p/miR1273p_CG.pdb"
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


#some initizlization for plotting
#get the x,y,z coordinates, molecule data
x = np.array([float(item) for item in df['x']])
y = np.array([float(item) for item in df['y']])
z = np.array([float(item) for item in df['z']])
#initialize empty arrays to plot boundary data
bx = []
by = []
bz = []

#get the boundary points
points = []
for i in range(0,len(df['x'])):
    points.append([float(df['x'][i]),float(df['y'][i]),float(df['z'][i])])
boundary_points = []
for pt in points:
    density_of_point = find_neighbour_density(pt,points)
    print(density_of_point)
    if(density_of_point['count']<15):
        boundary_points.append(pt)

#plot molecule with its boundary points
for i in range(0,len(boundary_points)):
    np.array(bx.append(boundary_points[i][0]))
    np.array(by.append(boundary_points[i][1]))
    np.array(bz.append(boundary_points[i][2]))

fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
# Creating plot
ax.scatter3D(x,y,z, color = "green")
ax.scatter3D(bx,by,bz, color = "red")
plt.title("simple 3D scatter plot")
plt.show()



#plot molecule with its boundary points
for i in range(0,len(boundary_points)):
    np.array(bx.append(boundary_points[i][0]))
    np.array(by.append(boundary_points[i][1]))
    np.array(bz.append(boundary_points[i][2]))

fig = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
# Creating plot
ax.scatter3D(x,y,z, color = "green")
ax.scatter3D(bx,by,bz, color = "red")
plt.title("simple 3D scatter plot")
plt.show()


# Just the boundary plot 
fig2 = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
# ax.scatter3D(x,y,z, color = "green")
ax.scatter3D(bx,by,bz, color = "red",s=2)
plt.title("simple 3D scatter plot")
ax.plot(bx,by,bz,color='red')
plt.xlabel('x')
plt.ylabel('y')
plt.show()




#get spline data
spline = []
bndry_pts = boundary_points[:]
total_points_length = len(bndry_pts)
# for i in range(0,total_points_length):
i = 0
while i<len(bndry_pts):
    neighbours_array = []
    #thee is thy's first neighbour
    neighbours_array.append(bndry_pts[i])
#     print(i)
#     for j,val in enumerate(bndry_pts[i+1:]):
    j = i+1
    while(j<len(bndry_pts)):
        dis = distance(bndry_pts[i],bndry_pts[j])
#         print(dis)
        if(dis<15):
            neighbours_array.append(bndry_pts[j])
            del bndry_pts[j]
            total_points_length = total_points_length - 1
        j = j+1
#     if(len(centroid_array)>0):
    if(len(neighbours_array)>=3):
        print(neighbours_array)
        #calculate centroid for the neighbours
        centre = centroid(neighbours_array)
#         print(centre)
        spline.append(centroid(neighbours_array))
#         print(spline)
    i = i+1


#Plot spine of molecule
#initialize arrays to plot the spine coordinates - centroidal axis
spline_x = []
spline_y = []
spline_z = []

print(spline)
for i in range(0,len(spline)):
    np.array(spline_x.append(spline[i][0]))
    np.array(spline_y.append(spline[i][1]))
    np.array(spline_z.append(spline[i][2]))
fig3 = plt.figure(figsize = (10, 7))
ax = plt.axes(projection ="3d")
# Creating plot
# ax.scatter3D(x,y,z, color = "green")
ax.scatter3D(bx,by,bz, color = "red",s=5)
ax.scatter3D(spline_x,spline_y,spline_z, color = "black",linestyle='-',s=5)
# ax.plot(spline_x,spline_y,spline_z, color='black')

plt.title("simple 3D scatter plot")
plt.xlabel('x')
plt.ylabel('y')
plt.show()


# In[19]:


#print
# spline
i = 140
j = 23
for line in spline:
    print("ATOM    "+str(i)+" CENT P   G  "+str(j)+"       "+str(line[0])+"  "+str(line[1])+" "+str(line[2])+"  1.00  0.00")
    i=i+1
    j=j+1


# #### Take only end points from spine and make a straight spine

#not necessary the end points in the spline list are the terminal points.
#add logic to detect end points == terminal points
p0 = spline[4]
pn = spline[9]
p1 = []
p2 = [] 
p3 = []
p4 = []
p5 = []
p6 = []
p7 = []


newspline = []
for i in range(0,3):
    p1.append(round (((p0[i]+pn[i])/2),2))
for i in range(0,3):
    p2.append(round(((p0[i]+p1[i])/2),2))
for i in range(0,3):
    p3.append(round(((p1[i]+pn[i])/2),2))
for i in range(0,3):
    p4.append(round(((p0[i]+p2[i])/2),2))
for i in range(0,3):
    p5.append(round(((p2[i]+p1[i])/2),2))
for i in range(0,3):
    p6.append(round(((p1[i]+p3[i])/2),2))
for i in range(0,3):
    p7.append(round(((p3[i]+pn[i])/2),2))
newspline.append(p0)
newspline.append(p1)
newspline.append(p2)
newspline.append(p3)
newspline.append(p4)
newspline.append(p5)
newspline.append(p6)
newspline.append(p7)
newspline.append(pn)



# spline
i = 140
j = 23
for line in newspline:
    print("ATOM    "+str(i)+" CENT P   G  "+str(j)+"     "+str(line[0])+"  "+str(line[1])+" "+str(line[2])+"  1.00  0.00")
    i=i+1
    j=j+1



# Add these lines to molecule_0.itp

i = 104
j = 17
for line in newspline:
    print(str(i)+" "+"SPN  "+str(j)+" "+"P CENT "+str(i)+"  0.0")
    i=i+1
    j=j+1


#increase points in spine
newspline = []
for i in range(0,len(spline)-1):
    mdpt = get_midpoint(spline[i],spline[i+1])
    newspline.append(spline[i])
    newspline.append(mdpt)
    newspline.append(spline[i+1])



j = 19
for line in newspline:
    print("ATOM    "+str(i)+" CENT P   G  "+str(j)+"     "+str(line[0])+"  "+str(line[1])+" "+str(line[2])+"  1.00  0.00")
    i=i+1
    j=j+1




