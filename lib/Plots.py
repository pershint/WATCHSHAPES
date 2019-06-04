import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
sns.set(font_scale=2)
import pandas
import scipy.optimize as spc

def OOR2(x, p1):
    return p1/(x**2)

def PlotFVLC_sphere(positions, LC_factors,detector_specs):
    '''For a given geometry, plot the light collection as a function of
    that axis (x,y,or z supported)'''
    sns.set_style("darkgrid")
    xkcd_colors = ['light eggplant', 'clay', 'leaf',
            'aqua blue','vomit', 'red','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
    R3_FV = detector_specs['radius']**3 
    data_dict = {'x': [], 'y': [], 'z': [], 'R3': [], 'LC': []} 
    axis_dict = {'x':0, 'y':1, 'z':2} 
    for a in axis_dict: 
        for j,lc in enumerate(LC_factors[0:10000]): 
            data_dict[a].append(int(positions[j][axis_dict[a]])) 
    for a in axis_dict:
        data_dict[a] = np.array(data_dict[a]) 
    for j, pos in enumerate(LC_factors[0:10000]):
        data_dict['R3'].append((np.sqrt(data_dict['x'][j]**2 + data_dict['y'][j]**2 + \
                data_dict['z'][j]**2)**3)/
                R3_FV)
        data_dict['LC'] = 100*np.array(LC_factors[0:10000])/(4.0*np.pi)
    plt.plot(data_dict['R3'],data_dict['LC'],marker='o',linewidth=0,markersize=8)
    plt.xlabel(r'$(R/R_{PMT})^{3}$',fontsize=24)
    plt.ylabel('PC%',fontsize=24)
    plt.title("Effective photocoverage in WATCHMAN",fontsize=30)
    plt.legend(fontsize=22)
    plt.ion()
    plt.show()

def PlotFVLC_cylinder_wdata(therealdeal,griddims=[20,20],hbounds=None):
    '''Takes the return from PlotFVLC_cylinder() and re-plots everything without
    having to send the data to a pandas dataframe again'''
    import seaborn as sns
    sns.set(font_scale=2)
    hm = therealdeal.pivot(index='zavg',columns='rho2', values='LC')
    #print(hm)
    if hbounds is None:
        ax = sns.heatmap(hm,cmap=plt.get_cmap('RdYlGn_r'),\
                cbar_kws={'label':r'$\Delta$PC'})
    else:
        ax = sns.heatmap(hm,cmap=plt.get_cmap('jet_r'),\
                cbar_kws={'label':'DeltaPC%'},vmin=hbounds[0],
                vmax=hbounds[1])
    plt.xlabel(r'$(\rho/\rho_{PMT})^{2}$',fontsize=24)
    plt.ylabel(r'$Z \, (mm)$',fontsize=24)
    plt.title("Effective photocoverage in WATCHMAN",fontsize=30)
    plt.legend(fontsize=22)
    plt.ion()
    plt.show()

def PlotFVLC_cylinder(datdict,griddims=[20,20],hbounds=None):
    '''For a given geometry, plot the light collection as a function of
    that axis (x,y,or z supported)'''
    positions = datdict["FV_positions"]
    LC_factors = datdict["light_factors"]
    detector_specs = datdict["detector_specs"]
    FV_specs = datdict["FV_specs"]
    rho2_FullVol = detector_specs['radius']**2 
    rho2_FV = FV_specs["radius"]**2
    data_dict = {'x': [], 'y': [], 'z': [], 'rho2': [], 'LC': []} 
    axis_dict = {'x':0, 'y':1, 'z':2} 
    for a in axis_dict: 
        for j,lc in enumerate(LC_factors): 
            data_dict[a].append(int(positions[j][axis_dict[a]])) 
    for j, pos in enumerate(LC_factors):
        data_dict['rho2'].append((data_dict['x'][j]**2 + data_dict['y'][j]**2)/
                rho2_FullVol)
    data_dict['LC'] = 100*np.array(LC_factors)/(4.0*np.pi)
    data_pd = pandas.DataFrame(data=data_dict)
    #We have to re-bin our data according to the input griddims
    rhobins = np.linspace(0, rho2_FV/rho2_FullVol, griddims[1]+1)
    print("RHOBINS: " + str(rhobins))
    rhobins_center = []
    for j,b in enumerate(rhobins):
        if j == 0:
            continue
        rbc = rhobins[j] - ((rhobins[j] - rhobins[j-1])/2)
        rhobins_center.append(rbc)
    zbins = np.linspace(data_pd.z.min(), data_pd.z.max(), griddims[0]+1)
    all_bins = []
    for j in xrange(len(zbins)):
        if j == 0: continue
        data_thisz = data_pd[data_pd.z > zbins[j-1]]
        data_thisz = data_thisz[data_pd.z < zbins[j]]
        data_thisz['zavg'] = pandas.Series(np.ones(len(data_thisz.z))*\
                np.mean(data_thisz.z),index=data_thisz.index)
        rho2_binned = data_thisz.groupby(pandas.cut(data_thisz.rho2, rhobins))
        rho2_binned = rho2_binned.aggregate(np.mean)
        #Now, set the rho2 value to the bin center, not the average
        rho2_binned.rho2 = rhobins_center
        rho2_binned.z = rho2_binned.z.astype(int)
        rho2_binned.zavg = rho2_binned.zavg.astype(int)
        rho2_binned.rho2 = np.round(rho2_binned.rho2,2)
        rho2_binned.LC = np.round(rho2_binned.LC,2)
        all_bins.append(rho2_binned) 
    therealdeal = None
    for j,b in enumerate(all_bins):
        if j==0: 
            therealdeal = b
        else:
            therealdeal=pandas.concat([therealdeal,b],ignore_index=True)
    print(therealdeal)
    hm = therealdeal.pivot(index='zavg',columns='rho2', values='LC')
    #print(hm)
    if hbounds is None:
        ax = sns.heatmap(hm,cmap=plt.get_cmap('ocean'),\
                cbar_kws={'label':'PC%'})
    else:
        ax = sns.heatmap(hm,cmap=plt.get_cmap('ocean'),\
                cbar_kws={'label':'PC%'},vmin=hbounds[0],
                vmax=hbounds[1])
    plt.xlabel(r'$(\rho/\rho_{PMT})^{2}$',fontsize=24)
    plt.ylabel(r'$Z \, (mm)$',fontsize=24)
    plt.title("Effective photocoverage in WATCHMAN",fontsize=30)
    plt.legend(fontsize=22)
    plt.ion()
    plt.show()
    return therealdeal
 

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
