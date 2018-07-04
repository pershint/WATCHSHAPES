import numpy as np
import playDarts as pd
import sys

class DetectorWorkshop(object):
    def __init__(self, PMTInformation=None, geometry=None, fiducial_mass=None):
        self.shape = {'geometry': geometry,'radius': None, 'height': None}
        self.fiducial_shape = {'geometry': geometry, 'mass':fiducial_mass, \
                'radius': None, 'height': None}
        if self.fiducial_shape['mass'] is not None:
            self._SetFiducialVolumeDimensions()
        if PMTInformation is not None: 
            self.PMTs = PMTInformation
            self.a = self.PMTs.PMTRadius
        #TODO: Associated this with a DB and a selected medium 
        self.attncf = (1/80000.0)

    def SetShape(self, shape):
        self.shape = shape

    def LoadPMTInfo(self, PMTInfo):
        self.PMTs = PMTInfo
        self.a = self.PMTs.PMTRadius
    
    def SetFiducialMassAndDimensions(self, FV):
        '''Define the detector's fiducial mass in tons'''
        self.fiducial_shape['mass'] = FV
        self._SetFiducialVolumeDimensions()

    def SetRadius(self,radius):
        self.shape['radius'] = radius

    def SetHeight(self, height):
        if self.shape['geometry'] is not "CYLINDER":
            print("Why are you adding a height if not a cylinder..")
        self.shape['height'] = height


    def _SetFiducialVolumeDimensions(self):
        '''Assuming water volume, define a fiducial volume in the
        center of the detector that's a scaled down version of the
        defined detector shape'''
        if self.fiducial_shape['mass'] is None:
            print("YOU NEED TO DEFINE YOUR FIDUCIAL MASS")
        #TODO: Add option for change in density here
        FV = self.fiducial_shape['mass']*1.0 #Units are m**3
        print("HEIGHT,RADIUS: %s,%s"%(self.shape['height'],self.shape['radius']))
        if self.shape['geometry'] == 'CYLINDER':
            height_radius_ratio = self.shape['height']/self.shape['radius']
            self.fiducial_shape['radius'] = (FV*1.0E9/(np.pi*height_radius_ratio))**\
                    (1.0/3.0)    
            self.fiducial_shape['height'] = self.fiducial_shape['radius'] * \
                    height_radius_ratio
        if self.shape['geometry'] == 'SPHERE':
            self.fiducial_shape['radius'] = pow(3*FV*1.0E9/(4*np.pi),(1.0/3.0))
        print("FIDUCIAL INFO: %s"%(str(self.fiducial_shape)))
        if self.fiducial_shape['radius'] > self.shape['radius']:
            warning = ("WARNING: YOU'VE DEFINED A FIDUCIAL VOLUME THATS"+
            "BIGGER THAN THE VOLUME ENCLOSED BY THE PMTS.  THE LIGHT"+
            "COLLECTION ALGORITHM WILL NOT WORK PROPERLY.")
            print(warning)

    def ShootFVPosition(self):
        '''Randomly shoot a position in the fiducial volume defined'''
        self._CheckFiducialDefinition()
        u = np.random.random()
        if self.shape['geometry'] == 'SPHERE':
            #Shooting from three gaussians of mean 0 and variance 1
            #Gives points evenly distributed on the unit sphere
            #after renormalization
            xyz = pd.RandShoot(0.0, 1.0, 3)
            xyz = xyz/np.sqrt(np.dot(xyz,xyz))
            #Now, You can't just shoot from 0 to 1 and multiply by
            #radius; you'll get too many events out at far radii.
            shot_radius = self.fiducial_shape['radius']*(u**(1.0/3.0))
            return shot_radius * xyz
        if self.shape['geometry'] == 'CYLINDER':
            #Shoot for radius, phi, and z.  Convert r,phi to x,y.
            r = self.fiducial_shape['radius']*(u**(1.0/2.0))
            phi = np.random.random() * 2.0 * np.pi
            z = np.random.random()*self.fiducial_shape['height']
            z = z - (self.fiducial_shape['height']/2.0)
            x = r * np.cos(phi)
            y = r * np.sin(phi)
            return np.array([x,y,z])

    def EvaluateLightCollection(self, position):
        '''For a single position, evaluate the light collection parameter'''
        print("EVALUATING LIGHT COLLECTION")
        light_collection_metric = 0.0 
        for i in xrange(len(self.PMTs.positions)):
            PMTPos = np.array(self.PMTs.positions[i])
            PMTDir = np.array(self.PMTs.directions[i])
            ObsPos = np.array(position)
            atn_factor = self._GetAttenuationFactor(PMTPos,ObsPos)
            solid_angle = self._GetExposedSolidAngle(PMTPos, PMTDir, ObsPos)
            light_collection_metric += (atn_factor*solid_angle)
        return light_collection_metric

    def _GetAttenuationFactor(self, PMTPos, ObsPos):
        '''Return the light reduction factor associated with
        the distance between a PMT and the Observation position'''
        pp = PMTPos 
        op = ObsPos 
        r = -pp + op
        r_mag = np.sqrt(np.dot(r,r))
        return np.exp(-r_mag*self.attncf) 
    
    def _GetExposedSolidAngle(self, PMTPos, PMTDir, ObsPos):
        '''Return the exposed solid angle of PMT surface from the
        PMT at position PMTPos [x,y,z] relative to ObsPos [x,y,z].'''
        #Turn given vectors into numpy arrays
        pd = PMTDir
        pp = PMTPos
        op = ObsPos
        r = -pp + op
        r_mag = np.sqrt(np.dot(r,r))
        pd_mag = np.sqrt(np.dot(pd,pd))
        #Define the angle and critical angles determining if theres
        #Surface area exposure loss for hemisphere glass PMTs
        theta = np.arccos(np.dot(r,pd)/(r_mag*pd_mag))
        theta_crit_low = np.arccos(np.sqrt(1 - (self.a**2/r_mag**2)))
        theta_crit_high = np.pi - theta_crit_low
        #We've gone through the calculations for the exposed
        #Solid angle based on theta and the theta_crits. Implement
        if theta > theta_crit_high:
            exp_solid_angle = 0
        elif theta < theta_crit_low:
            exp_solid_angle = (self.a**2/r_mag**2)*(1-(self.a**2/r_mag**2))
        else:
            exp_solid_angle = (self.a**2/2.0*r_mag**2)*\
                    (1-(self.a**2/r_mag**2))*(np.cos(theta)*\
                    np.sqrt((1-(self.a**2/r_mag**2)))+1)
        print("EXP_SA: %s"%(str(exp_solid_angle))) 
        return exp_solid_angle

    def _CheckFiducialDefinition(self):
        #Check necessary dimensions are defined
        VolumeNotDefined = False
        if self.fiducial_shape['mass'] is None:
            print("You must first define your detectors fiducial mass")
            VolumeNotDefined = True
        if self.fiducial_shape['radius'] is None:
            print("You must first define your detector's radius.")
            VolumeNotDefined = True
        if self.shape['geometry'] == 'CYLINDER':
            if self.fiducial_shape['height'] is None:
                print("You must first define your detector's cylinder height.")
                VolumeNotDefined = True
        if VolumeNotDefined:
            sys.exit(0)

