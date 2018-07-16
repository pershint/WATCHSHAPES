#We will populate the sphere using the Fibonacci Sphere algorithm, which gets
#Close to 

import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class GeometryPoints(object):
    def __init__(self, numpoints = 10):
        self.numpoints = numpoints
        self.positions = []
        self.geometry = "N/A"

    def ShowPositions(self):
        if len(self.positions)<=0:
            print("No positions have been filled yet.")
            return
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            x, y, z = [], [], []
            for i in xrange(len(self.positions)):
                x.append(self.positions[i][0])
                y.append(self.positions[i][1])
                z.append(self.positions[i][2])
            #X,Y = np.meshgrid(x, y)
            ax.scatter(x,y, z)
            ax.set_xlabel("X position (mm)", fontsize=22)
            ax.set_ylabel("Y position (mm)", fontsize=22)
            ax.set_zlabel("Z position (mm)", fontsize=22)
            for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(20)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
            plt.title("%s points distributed across a %s"%(\
                    str(self.numpoints),str(self.geometry)),fontsize=34)
            plt.show()

class SpherePoints(GeometryPoints):
    '''Given a radius (mm) and number of points, can generate
    positions equidistant on the sphere's surface using Fibonacci Sphere algorithm'''
    
    def __init__(self, numpoints = 10, radius=8000.0):
        super(SpherePoints, self).__init__(numpoints=numpoints)
        self.radius = radius
        self.geometry = "SPHERE"
    #source: https://stackoverflow.com/questions/9600801/evenly-distributing-
    #n-points-on-a-sphere
    def PopulatePositions(self, randomize=False):
        rnd = 1.
        if randomize is True:
            rnd = random.random() * self.numpoints

        points = []
        offset = 2./(self.numpoints)
        increment = np.pi * (3. - np.sqrt(5.)) #is the golden angle

        for i in range(self.numpoints):
            y = ((i * offset) - 1) + (offset/2)
            r = np.sqrt(1 - y**2)
            phi = ((i+rnd) % self.numpoints) * increment

            x = np.cos(phi) * r
            z = np.sin(phi) * r

            posx = x * self.radius
            posy = y * self.radius
            posz = z * self.radius
            points.append([posx,posy,posz])
        
        self.positions= points
        self.numpoints = len(self.positions)

class CylinderPoints(GeometryPoints):
    '''Given a radius (mm), height and number of points, can generate
    positions equidistant on a cylinder surface'''
    
    def __init__(self, numpoints = 10, radius=8000.0, height=16000.0):
        super(CylinderPoints, self).__init__(numpoints=numpoints)
        self.radius = radius
        self.height = height
        self.geometry = "CYLINDER"
        self.point_dict = {}

    def PopulatePositions(self):
        #First, we'll subtract a bit of the radius on the top, and the height
        #off to ensure there's no points overlapping/ on the edges
        topbot_rad = self.radius - (self.radius/20.0)
        strip_height = self.height - (self.height/20.0)

        #Find the fraction of points that should go on the cylinder top & bottom
        SA_side = 2.0*np.pi*self.radius*strip_height
        SA_tops = 2.0*np.pi*topbot_rad**2
        frac_side = SA_side/(SA_side + SA_tops)
        numpoints_side = int(self.numpoints * frac_side) + 1
        numpoints_topbot = int(self.numpoints * (1.0-frac_side))
        #For the cylinder side, we place the points at the center of
        #Squares distributed across the surface of the side
        side_points = []
        len_side = np.sqrt(SA_side/numpoints_side)
        #We're battling rounding errors here... I know it's silly
        phi, z = [], [] 
        phi = np.arange(0,2*np.pi,2*np.pi*(len_side/(2*np.pi*self.radius)))
        #Now, we have rounding errors that bring all values too low and give one extra row.
        #We need to shift all phis&zs by this rounding error, then delete the last point.
        corr = (2*np.pi - phi[len(phi)-1])
        for j,phival in enumerate(phi): 
            phi[j] = phi[j] + (float(j)/(len(phi)-1)*corr)
        phi = list(phi)
        del phi[len(phi)-1]
        phi = np.array(phi)
        z= np.arange(-strip_height/2.0,(strip_height/2.0), len_side)
        z_corr = ((strip_height/2.0) - z[len(z)-1])
        for j,zval in enumerate(z): 
            z[j] = z[j] + (float(j)/(len(z)-1)*z_corr)
        z = list(z)
        del z[len(z)-1]
        z = np.array(z)
        z = z + (len_side/2.0) 
        posx = self.radius * np.cos(phi)
        posy = self.radius * np.sin(phi)
        posz = z
        for j,zind  in enumerate(posz):
            for i,phi_ind in enumerate(posy):
                side_points.append([posx[i],posy[i],posz[j]])

        #Now, for the disk face, we're going to use Vogel's method to
        #populate each disk.  The reference for this source code is here:
        #blog.marmakoide.org/?p=1
        top_points, bot_points = [], []
        increment = np.pi * (3.0 - np.sqrt(5))  #The golden angle again
        #FIXME: Add a check in for if numpoints_topbot is an integer. if not,
        #Need to add one tube to either side for rounding correction
        p = int(numpoints_topbot/2.0)
        if numpoints_topbot%2 == 1:
            need_extratube = True
        for i in xrange(p):
            phi = i * increment
            r = topbot_rad*(np.sqrt(i) / np.sqrt(p))
            posx = r * np.cos(phi)
            posy = r * np.sin(phi)
            posz = self.height/2.0
            top_points.append([posx, posy, posz])
            bot_points.append([posx, posy, -posz])
        #add them all together
        points = side_points + top_points + bot_points
        self.positions = points
        self.numpoints = len(self.positions)
        self.point_dict = {"top": top_points, "bottom": bot_points,
                "side": side_points}

