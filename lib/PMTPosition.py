from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import PointGeometries as pg

class PMTGen(object):
    '''Class takes in a geometry of points and can build PMT direction and
    radius info based on the point geometry'''
    def __init__(self, PointGeometry=None,PMTRadius=254.0):
        self.PointGeometry = PointGeometry 
        self.positions = PointGeometry.positions
        self.directions = []
        self.geometry = PointGeometry.geometry
        print("GEO: %s"%(self.geometry))
        self.numPMTs = PointGeometry.numpoints
        self.PMTRadius = PMTRadius
        self.PMT_center_facing = None

    def SetPMTRadius(self, rad):
        '''Set the PMT radius in millimeters'''
        self.PMTRadius = rad

    def ShowPMTPositions(self):
        #FIXME: Memory error at meshgrid, supes broke
        if len(self.positions)<=0 or len(self.directions)<=0:
            print("No positions or directions have been filled yet.")
            return
        else:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            x, y, z = [], [], []
            for i in np.arange(0,len(self.positions), len(self.positions)/50):
                x.append(self.positions[i][0])
                y.append(self.positions[i][1])
                z.append(self.positions[i][2])
            x,y,z = np.meshgrid(np.array(x),np.array(y),np.array(z))
            u, v, w = [], [], []
            for i in np.arange(0,len(self.positions),len(self.positions)/50):
                u.append(self.positions[i][0])
                v.append(self.positions[i][1])
                w.append(self.positions[i][2])
            ax.quiver(x,y, z, u, v, w)
            ax.set_xlabel("X position (mm)")
            ax.set_ylabel("Y position (mm)")
            ax.set_zlabel("Z position (mm)")
            plt.title("PMT positions and directions for given geometry")
            plt.show()

    def BuildPMTDirections(self, face_center=True):
        '''Fill self.directions with the direction of each PMT.
        if face_center==False, PMTs face normal from wall of geometry'''
        if self.geometry == "SPHERE" and face_center is False:
            print("spheres always face the center...")
            face_center = True
        self.PMT_center_facing = face_center 
        if face_center is True:
            for PMTPos in self.positions:
                norm = np.sqrt(PMTPos[0]**2 + PMTPos[1]**2 + PMTPos[2]**2)
                u = -PMTPos[0] / norm
                v = -PMTPos[1] / norm
                w = -PMTPos[2] / norm
                self.directions.append([u,v,w])
        else:
            if self.geometry == "CYLINDER":
                #Use the point_dict to determine facing vectors
                self.positions = []
                pd = self.PointGeometry.point_dict
                for face in pd:
                    for pos in pd[face]:
                        norm = np.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
                        xy_norm = np.sqrt(pos[0]**2 + pos[1]**2)
                        self.positions.append(pos)
                        if face == "top": 
                            self.directions.append([0,0,-1])
                        if face == "bottom": 
                            self.directions.append([0,0,1])
                        if face == "side":
                            self.directions.append(\
                                    [-pos[0]/xy_norm,-pos[1]/xy_norm,0])

class SegSurfaceGen(PMTGen):
    def __init__(self,PointGeometry=None, PMTRadius=254.0, frac_exposed=(2.0-np.sqrt(2))):
        super(SegSurfaceGen,self).__init__(PointGeometry=PointGeometry,PMTRadius=PMTRadius)
        self.frac_exposed = frac_exposed
        self.SA = 4*np.pi*(self.PMTRadius**2)*self.frac_exposed
        self.PMTFacePosns = []
        self.PMTFaceDirs = []
        self.PMTPointRadiuses = []

    def setFracExposed(self):
        '''Set what fraction of the PMT 'sphere' is actually photocathode.
        default is currently at frac. of sphere exposed for 45 degrees'''

    def ShowPMTFacePosns(self):
        if len(self.positions)<=0 or len(self.directions)<=0:
            print("No positions or directions have been filled yet.")
            return
        else:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            x, y, z = [], [], []
            for i in np.arange(0,len(self.PMTFacePosns), len(self.positions)/50):
                for j in xrange(len(self.PMTFacePosns[i])):
                    x.append(self.PMTFacePosns[i][j][0])
                    y.append(self.PMTFacePosns[i][j][1])
                    z.append(self.PMTFacePosns[i][j][2])
            ax.scatter(x,y, z)
            ax.set_xlabel("X position (mm)", fontsize=22)
            ax.set_ylabel("Y position (mm)", fontsize=22)
            ax.set_zlabel("Z position (mm)", fontsize=22)
            for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(20)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
            plt.title("Points making up some of the PMT faces", fontsize=34)
            plt.show()
    
    def ShowPMTFaceDirs(self):
        if len(self.positions)<=0 or len(self.directions)<=0:
            print("No positions or directions have been filled yet.")
            return
        else:
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            x, y, z = [], [], []
            for i in np.arange(0,len(self.PMTFacePosns), len(self.PMTFacePosns)/3):
                for j in xrange(len(self.PMTFacePosns[i])):
                    x.append(self.PMTFacePosns[i][j][0])
                    y.append(self.PMTFacePosns[i][j][1])
                    z.append(self.PMTFacePosns[i][j][2])
            ax.set_xlabel("X position (mm)", fontsize=22)
            ax.set_ylabel("Y position (mm)", fontsize=22)
            ax.set_zlabel("Z position (mm)", fontsize=22)
            for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(20)
            for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
            for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
            x,y,z = np.meshgrid(np.array(x),np.array(y),np.array(z))
            u, v, w = [], [], []
            for i in np.arange(0,len(self.PMTFaceDirs),len(self.PMTFaceDirs)/3):
                u.append(self.PMTFaceDirs[i][j][0])
                v.append(self.PMTFaceDirs[i][j][1])
                w.append(self.PMTFaceDirs[i][j][2])
            ax.quiver(x,y, z, u, v, w)
            plt.title("PMT positions and directions for given geometry")
            plt.show()

    def BuildPMTSurfacePoints(self,pointresolution=1000):
        '''Function builds an array of points on a spherical surface of PMT radius
        covering only the fraction of PMT surface chosen to generate.  These are used
        to estimate the exposed surface area to any point in space'''
        
        self.PMTFacePosns = []
        self.PMTFaceDirs = []
            
        PMTSphere = pg.SpherePoints(numpoints=pointresolution,radius=self.PMTRadius)
        PMTSphere.PopulatePositions()
        SpherePosns = np.array(PMTSphere.positions)
        PMTDirections = np.array(self.directions)
        smag = self.PMTRadius
        thedot_test =np.sum(SpherePosns*PMTDirections[:,np.newaxis],axis=2)
        exp_check = (thedot_test/smag)
        ExposedPMTPoints = []
        PMTFaceDirs = []
        PMTPointRadiuses = np.sqrt(4.0 * self.PMTRadius**2 / pointresolution)
        for j,PMTSphereDotProds in enumerate(exp_check):
            indices = np.where(PMTSphereDotProds>=(1.0-2*self.frac_exposed))[0]
            ExposedPMT = np.array(SpherePosns[indices])
            #if each point is a disk on the PMT's surface, calculate the radius of each point assuming
            #it represents a disk that's some fraction of the surface area
            ExposedDir= ExposedPMT/(smag)
            #Normalized and points point outward from sphere in ExposedPMT,
            #So they're ready to represent a point's exposure direction
            PMTFaceDirs.append(ExposedDir)
            #Now, shift the points on the sphere all to be located around the PMT position
            ExposedPMT = ExposedPMT + self.positions[j]
            ExposedPMTPoints.append(ExposedPMT)
        self.PMTFacePosns=ExposedPMTPoints
        self.PMTFaceDirs=PMTFaceDirs
        self.PMTPointRadiuses = PMTPointRadiuses
        self.ShowPMTFacePosns()
#        self.ShowPMTFaceDirs()
#        for j,direc in enumerate(self.directions):
#            #Generate a sphere of points
#            PMTSphere = pg.SpherePoints(numpoints=1000,radius=self.PMTRadius)
#            PMTSphere.PopulatePositions()
#            SpherePosns = np.array(PMTSphere.positions)
#            smags = np.sqrt(np.sum(SpherePosns*SpherePosns,axis=1))
#            indices = np.where((np.sum(direc*SpherePosns,axis=1)/smags)<=
#                -(self.frac_exposed/2))[0]
#            PMTFacePosns = SpherePosns[indices]
#            PMTFaceDirs = PMTFacePosns*(1.0/np.sqrt(np.sum(PMTFacePosns*PMTFacePosns,axis=1)))[:,np.newaxis]
#            #Now, shift the points on the sphere all to be located around the PMT position
#            PMTFacePosns = PMTFacePosns + self.positions[j]
#            self.PMTFacePosns.append(PMTFacePosns)
#            self.PMTFaceDirs.append(PMTFaceDirs)
