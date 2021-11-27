# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 14:06:58 2021

@author: bread
CURVFAM
For Almighty Andrej. The Great and Powerful. 
Mightyest wizard in all the land of Kickassia.
This program prescribes a curvature to the universe and measures the obseravbles 
of uniform population in curved space.
2021 Danny Jensen
"""

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from astropy.coordinates import cartesian_to_spherical as cts
import itertools

###############################################################################
" Data Import and Processing "

data=[i for i in range(-10,0,2)]+[0]+[i for i in range(1,11,2)]

"Azimuthal Angle, Polar Angle, and Radial Distance"
location = [[[(i,j,k) for k in data] for j in data] for i in data]

x=list(i for i in data)
y=list(i for i in data)
z=list(i for i in data)
points=list(itertools.product(x,y,z))

###############################################################################
" Curvature Functions and Prime Location Function "

def kx(position):
    """
    Phi
    Position must be given as a vector
    """
    a=position[0]
    b=position[1]
    c=position[2]
    return a**2+c*b
def ky(position):
    """
    Theta
    Position must be given as a vector
    """
    a=position[0]
    b=position[1]
    c=position[2]
    return a**2-a*b
def kz(position):
    """
    Distance (time)
    Position must be given as a vector
    """
    a=position[0]
    b=position[1]
    c=position[2]
    return a**2-a-b
def primer(position):
    """
    Returns primed position.
    Position must be given as a vector.
    """
    x_prime=kx(position)
    y_prime=ky(position)
    z_prime=kz(position)
    prime=(x_prime,y_prime,z_prime)
    return prime

###############################################################################
" Length and Angle Functions "

def length(vector):
    "Returns the length of a vector"
    return np.linalg.norm(vector)
def unit(vector):
    "Returns the unit vector of the vector"
    if vector==(0,0,0):
        return (0,0,0)
    else:
        return vector/np.linalg.norm(vector)
def angle(v1, v2):
    "Returns the angle in degrees between 2 vectors"
    v1_u=unit(v1)
    v2_u=unit(v2)
    return np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

###############################################################################
" Rectangular Spherical Conversion Functions "

def rectangularize(position):
    "Given position vector in spherical returns cartesian"
    xx=position[0]
    yy=position[1]
    zz=position[2]
    x=zz*np.cos(xx)*np.sin(yy)
    y=zz*np.sin(xx)*np.sin(yy)
    z=zz*np.cos(yy)
    return x,y,z

def sphericalize(position):
    "Given position vector in cartesian returns spherical (azi,polar,dist)"
    xx=position[0]
    yy=position[1]
    zz=-position[2]
    r,t,p=cts(xx,yy,zz)
    r=r.value
    t=np.pi/2+t.value
    p=p.value
    return p,t,r

###############################################################################
" Find Adjacent Pairs "
    
adjacents=[]        # This is the list of adjacent pairs (tuple)
for i in range(0,len(location)-1):
    for j in range(0,len(location)-1):
        for k in range(0,len(location)-1):
            adj=(location[i][j][k],location[i+1][j][k])
            adjacents.append(adj)
            adj=(location[i][j][k],location[i][j+1][k])
            adjacents.append(adj)
            adj=(location[i][j][k],location[i][j][k+1])   
            adjacents.append(adj)
            
###############################################################################
" Find Angle Between Adjacent Pairs "

angles=[]       # This is the list of the flat space angles
for i in range(0,len(adjacents)):
    angle_value=angle(adjacents[i][0],adjacents[i][1])
    angles.append(angle_value)
    
plt.figure("Flat Space Adjacent Angles")
plt.hist(angles,bins=100,color="fuchsia")
plt.xlabel("Angle (deg)")
plt.ylabel("Number")
plt.title("Distribution of Angles Between Adjacent Points in Flat Space")
flat_y,flat_x= np.histogram(angles,bins=100)
flat_mean=np.mean(angles)
plt.vlines(flat_mean,ymin=0,ymax=max(flat_y),label="Mean Angle Value")
plt.legend()
print("Average angle value for flat space:", flat_mean, "degrees")

###############################################################################
" Apply Curvature to Adjacent Pairs "

primed_adjacents=[]     # This is the list of the adjacent pairs once curved (tuple)
spherical_adjacents=[]
for i in range(0,len(adjacents)):       # Make adjacents spherical
    value1=sphericalize(adjacents[i][0])
    value2=sphericalize(adjacents[i][1])
    value=(value1,value2)
    spherical_adjacents.append(value)
    
for i in range(0,len(spherical_adjacents)):     # Prime adjacents
    primed_value1=primer(spherical_adjacents[i][0])
    primed_value2=primer(spherical_adjacents[i][1])
    value=(primed_value1,primed_value2)
    primed_adjacents.append(value)
    
rect_primed_adjacents=[]    # Primed adjacents converted to rectangular
for i in range(0,len(primed_adjacents)):
    adj=(rectangularize(primed_adjacents[i][0]),rectangularize(primed_adjacents[i][1]))
    rect_primed_adjacents.append(adj)
    
primed_adjacents=rect_primed_adjacents
###############################################################################
" Find Angle Between Primed Adjacent Pairs "

primed_angles=[]       # This is the list of the curved space angles
for i in range(0,len(primed_adjacents)):
    angle_value=angle(primed_adjacents[i][0],primed_adjacents[i][1])
    primed_angles.append(angle_value)

plt.figure("Curved Space Adjacent Angles")
plt.hist(primed_angles,bins=100,color="fuchsia")
plt.xlabel("Angle (deg)")
plt.ylabel("Number")
plt.title("Distribution of Angles Between Adjacent Points in Curved Space")
primed_y,primed_x= np.histogram(primed_angles,bins=100)
primed_mean=np.mean(primed_angles)
plt.vlines(primed_mean,ymin=0,ymax=max(primed_y),label="Mean Angle Value")
plt.legend()
print("Average angle value for curved space:", primed_mean, "degrees")

###############################################################################
" 3-D Plot "

xx=[]
yy=[]
zz=[]
for i in points:
    xx.append(i[0])
    yy.append(i[1])
    zz.append(i[2])
xxx=[]
yyy=[]
zzz=[]
for i in range(0,len(xx)):
    xxx.append(zz[i]*np.cos(xx[i])*np.sin(yy[i]))
    yyy.append(zz[i]*np.sin(xx[i])*np.sin(yy[i]))
    zzz.append(zz[i]*np.cos(yy[i]))

plt.figure("3D")
ax = plt.axes(projection='3d')
ax.scatter3D(xx, yy, zz, color="fuchsia",s=5)
plt.title("Position in Flat Space")

###############################################################################
" Primed 3-D Plot "

x_prime,y_prime,z_prime=[],[],[]

spherical_points=[]

for i in points:
    spherical_points.append(sphericalize(i))

for i in range(0,len(spherical_points)):
    x,y,z=primer(spherical_points[i])
    x_prime.append(x)       # Phi
    y_prime.append(y)       # Theta
    z_prime.append(z)       # r
    
xx_prime,yy_prime,zz_prime=[],[],[]

for i in range(0,len(x_prime)):
    xx_prime.append(z_prime[i]*np.cos(x_prime[i])*np.sin(y_prime[i]))   # x
    yy_prime.append(z_prime[i]*np.sin(x_prime[i])*np.sin(y_prime[i]))   # y
    zz_prime.append(z_prime[i]*np.cos(y_prime[i]))                      # z

plt.figure("3D Primed")
ax = plt.axes(projection='3d')
ax.scatter3D(xx_prime, yy_prime, zz_prime, color="fuchsia",s=5)
plt.title("Position in Curved Space")

###############################################################################
" Distance Histograms "

distances=[]
for i in range(0,len(points)):
    distances.append(length(points[i]))

plt.figure("Flat Space Distances")
plt.hist(distances,bins=100,color="fuchsia")
plt.xlabel("Distance (units)")
plt.ylabel("Number")
plt.title("Distribution of Distances From Reference Point in Flat Space")
dist_x,dist_y= np.histogram(distances,bins=100)
dist_mean=np.mean(distances)
#plt.vlines(dist_mean,ymin=0,ymax=max(dist_x),label="Mean Distance Value")
#plt.legend()
print("Average distance value for flat space:", dist_mean, "units")

primed_distances=[]
for i in range(0,len(points)):
    primed_distances.append(length(rectangularize(primer(spherical_points[i]))))

plt.figure("Curved Space Distances")
plt.hist(primed_distances,bins=100,color="fuchsia")
plt.xlabel("Distance (units)")
plt.ylabel("Number")
plt.title("Distribution of Distances From Reference Point in Curved Space")
primed_dist_x,primed_dist_y= np.histogram(primed_distances,bins=100)
primed_dist_mean=np.mean(primed_distances)
#plt.vlines(primed_dist_mean,ymin=0,ymax=max(primed_dist_x),label="Mean Distance Value")
#plt.legend()
print("Average distance value for curved space:", primed_dist_mean, "units")

###############################################################################

# Code is currently a work in progess
# Forgive any clutter or complications
