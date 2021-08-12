#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 15:08:44 2021

@author: Aron Dahlberg
"""

import numpy as np
import matplotlib.pyplot as plt




## FOR USER TO DEFINE: 


## Closest approach distance given in kilometers.
d=5000


##Velocity in km/s
v=70

## Fly-by time interval of the detailed imaging interval and outer interval. If the detailed imaging interval is from -120 to 120 seconds, give the parameter
## as 120. Seconds_per parameters tell how many seconds are between two images within a given interval.

detailed_flyby_time=120
seconds_per_image_detailed=10

outer_flyby_time=300
seconds_per_image_outer=60


## Sun angle in relation to the closest approach vector, shown as +- 45 degrees in figure (number) of the documentation.
sun_angle=0

##Angle of inclination in degrees. Negative angle means that the trajectory is below the sun and positive above. 
gamma=0

##Restrictions of what coordinates are taken into consideration during the flyby. First term indicates the angle "delta" on the left side of closest approach
##and the second term the angle on the right side (Figure (referoi)). If you do not want any restrictions, insert values -0,0. Note: keep the second term 
##negative to indicate the side.
delta_restrictions=[0,-0]

## Other things user might want to change, if necessary:

## Sun Irradiance and rotation of the asteroid, used in Mario's code.

sun_irr=250

rot=60




## HERE BEGIN'S THE CODE, WHICH USER DOES NOT HAVE TO CHANGE.


## Closest-approach vector p. Keep this as 0,1,0 if you use spherical coordinates in Mario's code.
closest_approach=np.array([0,1,0])*d


## These are used in calculating the perpendicular vector, which is used to construct the trajectory. 
p_x=closest_approach[0]
p_y=closest_approach[1]
p_z=closest_approach[2]



    
## Calculating a vector which is perpendicular to the closest approach vector. Derivation of these equations can be found in the documentation.
if (p_y!=0):
    a=np.sqrt(1/(1+p_x**2/p_y**2))

    b=-(p_x/p_y)*a

else:
    a=0
    b=d


c=0



## Unit-length size perpendicular vector w which is then scaled,
## so that the scaled vector + the closest approach vector determines the satellite's initial position.
unit_w=np.array([a,b,c])



## Sun vector, which is set to be pointing to y-direction (scale of 15000 has no significance, could be anything. Just looked good in the plots). 
sun_vector=(closest_approach/d)*15000



## Sun angle in radians and the rotation matrix, which rotates the whole system by the angle specified with the parameter sun_angle. Explanation in
## the documentation.
sun_angle=np.deg2rad(sun_angle)
rotation_matrix=np.array([[np.cos(sun_angle),-np.sin(sun_angle)],[np.sin(sun_angle),np.cos(sun_angle)]])


## Function which rotates the whole system.
def RotateSystem(closest_app,sat_vec, rotation_matrix):
   
    closest_app=np.append(rotation_matrix.dot(np.array([closest_app[0],closest_app[1]])),closest_app[2])
    sat_vec=np.append(rotation_matrix.dot(np.array([sat_vec[0],sat_vec[1]])),sat_vec[2])
    return closest_app, sat_vec


## Rotated system parameters.
    
closest_approach, unit_w = RotateSystem(closest_approach,unit_w,rotation_matrix)



def RotationVertical(closest_approach, gamma):
    gamma=np.deg2rad(gamma)
    help_vector=np.cross(closest_approach,unit_w)
    closest_approach_rotated=closest_approach*np.cos(-gamma) + help_vector*np.sin(-gamma)
    return closest_approach_rotated


closest_approach=RotationVertical(closest_approach,gamma)
    
##Allocating a list which contains time points of fly-by images. The while loop creates the coordinates to be 
## rendered according to velocity and the time parameters.

t=-(outer_flyby_time + detailed_flyby_time)

time=[[t]]

rendering_coordinates=[]

while t<=(outer_flyby_time + detailed_flyby_time):
    
    if(t<=0):
        delta=np.rad2deg(np.arcsin(d/np.linalg.norm(closest_approach + unit_w*v*t)))
        statement=delta>= delta_restrictions[0]
        
        
     
    else:
        delta=-np.rad2deg(np.arcsin(d/np.linalg.norm(closest_approach + unit_w*v*t)))
        statement= delta<=delta_restrictions[1]
       
    
    
    if (t<-1*(detailed_flyby_time) or detailed_flyby_time<=t) and (statement==True):
        rendering_coordinates.append(closest_approach + unit_w*v*t)
        time[0].append(t)
        
        t=t+seconds_per_image_outer
        
        
    elif (t<-1*(detailed_flyby_time) or detailed_flyby_time<=t) and (statement==False):
         
         t=t+seconds_per_image_outer
         
        
    
    
    
    elif (-detailed_flyby_time<=t<detailed_flyby_time) and (statement==True):
        rendering_coordinates.append(closest_approach + unit_w*v*t)
        time[0].append(t)
        
        t=t+seconds_per_image_detailed
        
          
    elif (-detailed_flyby_time<=t<detailed_flyby_time) and (statement==False):
        
        t=t+seconds_per_image_detailed
   
    
rendering_coordinates=np.array(rendering_coordinates)

time=np.append(np.array([["time (s)"]]),np.transpose(np.array([time[0][1:]])),axis=0)




##Phase angle and distance from satellite to asteroid.
phi_and_distance=[]

np.rad2deg(np.arccos(np.dot(closest_approach,sun_vector)/(np.linalg.norm(closest_approach)*np.linalg.norm(sun_vector))))

##This is the time value, when phase angle is zero. Crucial to defining when phi changes signs.
time_of_index=np.tan(sun_angle)*d/(np.linalg.norm(unit_w*v))

for k in range(len(rendering_coordinates)):
    phi=np.rad2deg(np.arccos(np.dot(rendering_coordinates[k],sun_vector)/(np.linalg.norm(rendering_coordinates[k])*np.linalg.norm(sun_vector)))) 
    
    if(float(time[k+1])<=time_of_index):
        phi=phi
    else:
        phi=-phi
      
    distance=np.linalg.norm(rendering_coordinates[k])
    
    phi_and_distance.append([phi,distance])
    

 
phi_and_distance=np.append(np.array([["phi (deg)","distance(km)"]]),np.array(phi_and_distance),axis=0)




## These next 30 lines make the data files to look as wanted. All of the divisions by 1000 are because Mario's code takes inputs of 
## coordinates in 1:1000 ratio. Meaning that for example value camera_x = 10 means 10 000 km.
camera_lon=np.zeros((len(rendering_coordinates),1))
camera_spherical=np.concatenate((camera_lon,phi_and_distance[1:]),axis=1)

##Allocate the postionalData files in cartesian and spherical coordinates, that Mario's code can read.
positionalData_cart=np.array([["camera_x", "camera_y", "camera_z", "sun_x", "sun_y", "sun_z", "sun_irr", "rot"]])
positionalData_spherical=np.array([["camera_lat", "camera_lon", "camera_dist", "sun_lat", "sun_lon", "sun_irr", "rot"]])


##dont worry about this. This is just to make a list of length len(renering_coordinates) of these variables, that do not change.
##Mario's code needs these for each set of satellite points, thats why.
inputs_cart=np.array([[sun_vector[0]/1000,sun_vector[1]/1000,sun_vector[2]/1000, sun_irr, rot]])
inputs_spherical=np.array([[0,0,sun_irr,rot]])

for l in range(len(rendering_coordinates)-1):
    
    inputs_cart=np.append(inputs_cart,np.array([[sun_vector[0]/1000,sun_vector[1]/1000,sun_vector[2]/1000, sun_irr, rot]]),axis=0)
    inputs_spherical=np.append(inputs_spherical,np.array([[0,0,sun_irr,rot]]),axis=0)

##Fly-by data for Mario's code to read in cartesian coordinates.
positionalData_cart=np.append(positionalData_cart,(np.append(rendering_coordinates/1000,inputs_cart,axis=1)),axis=0)

##Fly-by data for Mario's code to read in spherical coordinates.
positionalData_spherical=np.append(positionalData_spherical,(np.append(camera_spherical,inputs_spherical,axis=1)),axis=0)


phi_and_distance=np.concatenate((time,phi_and_distance),axis=1)








##Writing the data to text-files.

#Cartesian coordinate data
try:
    textfile=open('fly_by_coordinates_cart.txt','w')
    for row in positionalData_cart:
       row=np.reshape(row,(1,positionalData_cart.shape[1])).tolist()
       string="         ".join(str(item) for item in row[0])  
       textfile.write(string + '\n')
    textfile.close()
    print("Text file containing fly-by path in cartesian coordinates has been succesfully created with the name fly_by_coordinates_cart.txt")
    
except:
    print("Writing the text file fly_by_coordinates_cart has failed.")

print("\n")   

##Spherical coordinate data
try:
    textfile=open('fly_by_coordinates_spherical.txt','w')
    for row in positionalData_spherical:
       row=np.reshape(row,(1,positionalData_spherical.shape[1])).tolist()
       string="         ".join(str(item) for item in row[0])  
       textfile.write(string + '\n')
    textfile.close()
    print("Text file containing fly-by path in spherical coordinates has been succesfully created with the name fly_by_coordinates_spherical.txt")
    
except:
    print("Writing the text file fly_by_coordinates_spherical has failed.")
    
print("\n")   
##Phase angle and distance from satellite to asteroid
try:
    textfile=open('phi_and_distance.txt','w')
    for row in phi_and_distance:
       row=np.reshape(row,(1,3)).tolist()
       string="         ".join(str(item) for item in row[0])  
       textfile.write(string + '\n')
    textfile.close()
    print("Text file containing the phase angle phi and distance from satellite to asteroid has been succesfully created with the name phi_and_distance.txt")

except:
    print("Writing the text file delta_and_distance has failed.")




##Plotting the satellites initial vector, closest approach vector and the fly-by trajectory.
##NOTE: only for visualization purposes, so that the user gets the idea behind these calculations. I will also make a diagram of the 
##why this calculates the trajectory.

##origin, satellites initial position and satellites last position for plotting vectors.
X,Y,Z=0,0,0


sat_ini_x,sat_ini_y,sat_ini_z=rendering_coordinates[0][0],rendering_coordinates[0][1],rendering_coordinates[0][2]
sat_last_x,sat_last_y,sat_last_z=rendering_coordinates[-1][0],rendering_coordinates[-1][1],rendering_coordinates[-1][2]
## The scaling factor is the distance which the satellite travels at said velocity in the given time interval to both sides around t=0. Unit is in kilometers.
scaling_factor=v*(outer_flyby_time+detailed_flyby_time)
w=unit_w*scaling_factor*2

fig=plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_xlim(-scaling_factor-1000,scaling_factor+1000)
ax.set_ylim(-scaling_factor-1000,scaling_factor+1000)
ax.set_zlim(-scaling_factor-1000,scaling_factor+1000)

ax.quiver(X,Y,Z,sat_ini_x,sat_ini_y,sat_ini_z,color='k',arrow_length_ratio=0.3,label="satellite's initial position vector")
ax.quiver(sat_ini_x,sat_ini_y,sat_ini_z,w[0],w[1],w[2], color='r', arrow_length_ratio=0.1,label="satellite's trajectory")
ax.quiver(X,Y,Z,closest_approach[0],closest_approach[1],closest_approach[2],color='g', arrow_length_ratio=0.4, label="closest approach vector")


ax.quiver(X,Y,Z, sun_vector[0],sun_vector[1],sun_vector[2],color='y',arrow_length_ratio=0.3, label="sun vector")

ax.set_title("Fly-by trajectory with trajectory lifted by "+ str(gamma)+" degrees from xy-plane",y=1.1)
plt.legend(loc = 'best')
plt.show()




