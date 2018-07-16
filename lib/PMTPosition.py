from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

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
