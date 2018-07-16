import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import scipy.optimize as spc

def OOR2(x, p1):
    return p1/(x**2)

def PlotLightCollection(positions, LC_factors, axis=None):
    '''For a given axis, plot the light collection as a function of
    that axis (x,y,or z supported)'''
    sns.set_style("whitegrid")
    xkcd_colors = ['slate blue','black']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    axis_dict = {'x':0, 'y':1, 'z':2} 
    for a in axis_dict: 
        if axis==a:
            theaxis = axis_dict[a] 
    the_positions = [] 
    for i in xrange(len(positions)):
        the_positions.append(positions[i][theaxis])
    the_positions = np.array(the_positions)
    LC_factors = np.array(LC_factors)
    plt.plot(the_positions, LC_factors, linestyle='none',marker='o',markersize=7,
            label="LC factor")
    plt.legend(loc=1) 
    plt.show()


def PlotLightCollection_OnePMT(positions, LC_factors, axis=None):
    '''For a given axis, plot the light collection as a function of
    that axis (x,y,or z supported)'''
    sns.set_style("whitegrid")
    xkcd_colors = ['slate blue','black']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))
    axis_dict = {'x':0, 'y':1, 'z':2} 
    for a in axis_dict: 
        if axis==a:
            theaxis = axis_dict[a] 
    the_positions = [] 
    for i in xrange(len(positions)):
        the_positions.append(positions[i][theaxis])
    the_positions = np.array(the_positions)
    the_distance = 8203.0 - the_positions
    LC_factors = np.array(LC_factors)
    plt.plot(the_distance, LC_factors, linestyle='none',marker='o',markersize=7,
            label="LC factor")
    popt, pcov = spc.curve_fit(OOR2, the_distance, LC_factors, p0=[0.01])
    print("BEST FIT VALUES: " + str(popt))
    print("PCOVARIANCE: " + str(pcov))
    x = np.arange(min(the_distance),max(the_distance),
            (max(the_distance)-min(the_distance))/100.0)
    bfvals = OOR2(x, popt[0]) 
    plt.plot(x, bfvals, linewidth=4,label=r'$A/r^{2}$ fit')
    plt.legend(loc=1) 
    plt.show()

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
        plt.title("Distribution of positions in input array")
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
