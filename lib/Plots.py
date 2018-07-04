import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

def ShowPositions(positions):
    if len(positions)<=0:
        print("No positions have been filled yet.")
        return
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x, y, z = [], [], []
        for i in xrange(len(positions)):
            x.append(positions[i][0])
            y.append(positions[i][1])
            z.append(positions[i][2])
        #X,Y = np.meshgrid(x, y)
        ax.scatter(x,y, z,label='PMT positions')
        ax.set_xlabel("X position (mm)")
        ax.set_ylabel("Y position (mm)")
        ax.set_zlabel("Z position (mm)")
        plt.title("%s points distributed across a %s"%(\
                str(self.numpoints),str(self.geometry)))
        plt.legend() 
        plt.show()

def ContourMap_XYSlice(positions,light_factors,zrange=[-1000.0,1000.0],pmt_positions=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x, y, lf = [], [], []
    for i in xrange(len(positions)):
        if zrange[0] < positions[i][2] < zrange[1]:
            x.append(positions[i][0])
            y.append(positions[i][1])
            lf.append(light_factors[i])
    if pmt_positions is not None: 
        for i in xrange(len(pmt_positions)):    
            px, py, pz = [], [], []
            for i in xrange(len(positions)):
                px.append(positions[i][0])
                py.append(positions[i][1])
                pz.append(positions[i][2])
            #X,Y = np.meshgrid(x, y)
            ax.scatter(px,py, pz,label='PMT positions')
    ax.plot_trisurf(x,y,light_factors,cmap=plt.cm.jet, linewidth=0.2,label='LC factor')

    ax.set_xlabel("X position (mm)")
    ax.set_ylabel("Y position (mm)")
    ax.set_zlabel("LC Factor")
    plt.title("Light collection factor through WATCHMAN fiducial volume\n"+\
            "%s points in slice; avg. LC factor: %s\n"%(str(len(lf)),\
            str(np.average(lf)))+"zrange (mm): %s"%(str(zrange)))
    plt.legend() 
    plt.show()

def ColorMap(positions,light_factors):
    if len(positions)<=0:
        print("No positions have been filled yet.")
        return
    else:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x, y, z = [], [], []
        for i in xrange(len(positions)):
            x.append(positions[i][0])
            y.append(positions[i][1])
            z.append(positions[i][2])
        ax.scatter(x,y,z,np.array(light_factors)*1000.0,cmap=plt.cm.spring)
        ax.set_xlabel("X position (mm)")
        ax.set_ylabel("Y position (mm)")
        ax.set_zlabel("Z position (mm)")
        plt.title("Light collection factor through WATCHMAN fiducial volume\n"+\
                "%s points presented"%(str(len(light_factors))))
        plt.show()
