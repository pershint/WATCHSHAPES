import numpy as np
import playDarts as pd
import sys

class DetectorWorkshop(object):
    def __init__(self, PMTInformation=None, geometry=None, fiducial_mass=None,density=None,
            attenuation_length=None):
        self.shape = {'geometry': geometry,'radius': None, 'height': None,
                'PMT_center_facing':None}
        self.fiducial_shape = {'geometry': geometry, 'mass':fiducial_mass, \
                'radius': None, 'height': None}
        self.medium = {'density': density, 'attenuation_length':attenuation_length}
        if self.fiducial_shape['mass'] is not None:
            self._SetFiducialVolumeDimensions()
        if PMTInformation is not None: 
            self.PMTs = PMTInformation
            self.a = self.PMTs.PMTRadius
            self.shape['PMT_center_facing']= self.PMTs.PMT_center_facing
        #TODO: Associated this with a DB and a selected medium 
        if density is not None:
            self.density = density
            self.medium['density'] = density
        if attenuation_length is not None:
            self.medium['attenuation_length'] = attenuation_length
            self.attncf = (1/attenuation_length)

    def LoadPMTInfo(self, PMTInfo):
        self.PMTs = PMTInfo
        self.a = self.PMTs.PMTRadius
        self.shape['PMT_center_facing'] = self.PMTs.PMT_center_facing

    def SetAttenuationLength(self, attnlength):
        '''set the mediums attenuation length in mm''' 
        self.attncf = (1.0/attnlength)
        self.medium['attenuation_length'] = attnlength

    def SetMediumDensity(self, dens):
        '''set the detector's medium density in kg/m**3'''
        self.density = dens
        self.medium['density'] = dens

    def SetGeometry(self, geo):
        self.shape['geometry'] = geo

    def SetFiducialMassAndDimensions(self, FV):
        '''Define the detector's fiducial mass in tons'''
        self.fiducial_shape['mass'] = FV
        self._SetFiducialVolumeDimensions()

    def SetRadius(self,radius):
        self.shape['radius'] = radius

    def GetDetectorRadius(self, buff):
        '''Given the input buffer volume in mm, returns
        the detector radius in mm'''
        self._CheckFiducialDefinition()
        return self.fiducial_shape['radius'] + buff
    
    def GetDetectorHeight(self, buff):
        '''Given the calculated fiducial volume dimensions, returns
        the detector radius for a cylinder configuration'''
        self._CheckFiducialDefinition()
        return self.fiducial_shape['height'] + 2.0*buff

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
            return
        #TODO: Add option for change in density here
        if self.density is None:
            print("YOU MUST FIRST DEFINE YOUR DETECTOR DENSITY")
            sys.exit(0)
        FV = self.fiducial_shape['mass']*self.density #Units are m**3
        print("HEIGHT,RADIUS: %s,%s"%(self.shape['height'],self.shape['radius']))
        if self.shape['geometry'] == 'CYLINDER':
            if self.shape['radius'] is None or self.shape['height'] is None:
                print("NO RADIUS AND/OR HEIGHT DEFINED IN WORKSHOP. ASSUMING RIGHT "+\
                        "CYLINDER WITH HEIGHT=2*RADIUS")
                height_radius_ratio = 2.0
            else:
                height_radius_ratio = self.shape['height']/self.shape['radius']
            self.fiducial_shape['radius'] = (FV*1.0E9/(np.pi*height_radius_ratio))**\
                    (1.0/3.0)    
            self.fiducial_shape['height'] = self.fiducial_shape['radius'] * \
                    height_radius_ratio
        if self.shape['geometry'] == 'SPHERE':
            self.fiducial_shape['radius'] = pow(3*FV*1.0E9/(4*np.pi),(1.0/3.0))
        print("FIDUCIAL INFO: %s"%(str(self.fiducial_shape)))
        if self.shape['radius'] is not None: 
            if self.fiducial_shape['radius'] > self.shape['radius']:
                warning = ("WARNING: YOU'VE DEFINED A FIDUCIAL VOLUME THATS"+
                "BIGGER THAN THE VOLUME ENCLOSED BY THE PMTS.  THE LIGHT"+
                "COLLECTION ALGORITHM WILL NOT WORK PROPERLY.")
                print(warning)

    def ShootPosition(self, plane=None, axis=None,region='FV'):
        '''Randomly shoot a position in the fiducial volume defined'''
        self._CheckFiducialDefinition()
        u = np.random.random()
        axis_dict = {'x':0, 'y':1, 'z':2} 
        if region=='FV':
            radius = self.fiducial_shape['radius']
            height = self.fiducial_shape['height']
        if region=='PMT':
            radius = self.shape['radius']
            height = self.shape['height']
        if self.shape['geometry'] == 'SPHERE':
            if plane=='zeq0':
                #Sample from a disc of radius R 
                z = 0.0
                r = radius*(u**(1.0/2.0))
                phi = np.random.random() * 2.0 * np.pi
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                xyz = np.array([x,y,z])
            elif axis is not None:
                for a in axis_dict: 
                    if axis == a:
                        xyz = np.zeros(3) 
                        thevar = ((np.random.random() * radius*2.0) - radius)
                        xyz[axis_dict[a]] = thevar 
            else:
                #Shooting from three gaussians of mean 0 and variance 1
                #Gives points evenly distributed on the unit sphere
                #after renormalization
                xyz = pd.RandShoot(0.0, 1.0, 3)
                xyz = xyz/np.sqrt(np.dot(xyz,xyz))
                #Now, You can't just shoot from 0 to 1 and multiply by
                #radius; you'll get too many events out at far radii.
                shot_radius = radius*(u**(1.0/3.0))
                xyz = shot_radius*xyz 
            return xyz
        if self.shape['geometry'] == 'CYLINDER':
            #Shoot for radius, phi, and z.  Convert r,phi to x,y.
            r = radius*(u**(1.0/2.0))
            phi = np.random.random() * 2.0 * np.pi
            if plane is not None:
                if plane=='side':
                    z = np.random.random()*height
                    z = z - (height/2.0)
                    r = radius
                if plane=='top':
                    z = height/2.0
                if plane=='bottom':
                    z = -height/2.0
                if plane=='zeq0':
                    z = 0.0 
                x = r * np.cos(phi)
                y = r * np.sin(phi) 
                return np.array([x,y,z])
            elif axis is not None:
                for a in axis_dict: 
                    if axis == a:
                        xyz = np.zeros(3) 
                        if a == 'x' or a == 'y': 
                            thevar = ((np.random.random() * radius*2.0) - radius)
                        elif a == 'z':
                            thevar = ((np.random.random() * height) - height/2.0)
                        xyz[axis_dict[a]] = thevar 
                return xyz
            elif axis is None and plane is None:
                z = np.random.random()*height
                z = z - (height/2.0)
                x = r * np.cos(phi)
                y = r * np.sin(phi)
                return np.array([x,y,z])
            else:
                print("ONLY DEFINE AN AXIS OR PLANE TO SHOOT POINTS ON")

    def EvaluateLightCollection(self, position,UseSegSurface=False):
        '''For a single position, evaluate the light collection parameter'''
        PMTPos = np.array(self.PMTs.positions)
        PMTDir = np.array(self.PMTs.directions)
        ObsPos = np.array(position)
        atn_factors = self._GetAttenuationFactors(PMTPos,ObsPos)
        if UseSegSurface is False:
            solid_angles = self._GetExposedSolidAngles(PMTPos, PMTDir, ObsPos)
        else:
            PMTFacePosns= self.PMTs.PMTFacePosns
            PMTFaceDirs = self.PMTs.PMTFaceDirs
            PMTFaceExposed = self.PMTs.frac_exposed
            PMTPointRadiuses = self.PMTs.PMTPointRadiuses
            solid_angles = self._GetExposedSolidAngles_SS(PMTFacePosns,PMTFaceDirs,
                    PMTPointRadiuses,ObsPos,PMTFaceExposed)
        light_collection_metric = np.sum(atn_factors*solid_angles)
        return light_collection_metric

    def _GetAttenuationFactors(self, PMTPos, ObsPos):
        '''Return the light reduction factors associated with
        the distance between each PMT and the Observation position'''
        pp = PMTPos 
        op = ObsPos 
        r = -pp + op
        r_mags = np.sqrt(np.sum(r*r,axis=1)) 
        attn_factors = np.ones(len(r_mags))*self.attncf
        return np.exp(-r_mags*attn_factors) 
   
    def _GetExposedSolidAngles_SS(self,PMTFaces,PMTFaceDirs,PMTPointRadiuses,ObsPos,frac_exposed):
        '''Return the exposed solid angle of PMT surface from each PMT
        face using the populated PMT dot positions relative to ObsPos [x,y,z]'''
        solid_angles = []
        #FIXME: Have to speed up shit using the [:,np.newaxis] so there's
        #no for loop
        for j,face in enumerate(PMTFaces):
            vecs_to =  (ObsPos - face)
            vecs_to_dist = np.array(np.sqrt(np.sum(vecs_to*vecs_to,axis=1)))
            vecs_to_norm = vecs_to*(1/np.sqrt(np.sum(vecs_to*vecs_to,axis=1)))[:,np.newaxis]
            vecdotdir = np.sum(vecs_to_norm*PMTFaceDirs[j],axis=1)
            isvis_indices =np.where(vecdotdir>=0)[0]
            
            visible_dists = vecs_to_dist[isvis_indices]
            visible_costhetas = np.array(vecdotdir[isvis_indices])
            #Now, assume each point on the sphere is a disk with it's respective
            #fraction of the total PMT's surface area.  Calculate the exposed
            #solid angle for that disk
            solid_angle_pmtsurf = 2*np.pi*(1 - (visible_dists/(np.sqrt(PMTPointRadiuses**2 + \
                    visible_dists**2))))*visible_costhetas
            #Total solid angle visible on this PMT is the sum of all the above points
            solid_angles.append(np.sum(solid_angle_pmtsurf))
        return solid_angles


    def _GetExposedSolidAngles(self, PMTPos, PMTDir, ObsPos):
        '''Return the exposed solid angle of PMT surface from each PMT
        position PMTPos (array of [x,y,z]) relative to the ObsPos [x,y,z].'''
        #Turn given vectors into numpy arrays   
        pd = PMTDir
        pp = PMTPos
        op = ObsPos
        r = -pp + op
        r_mags = np.sqrt(np.sum(r*r,axis=1)) 
        pd_mags = np.sqrt(np.sum(pd*pd,axis=1)) 
        a = np.ones(len(PMTPos))*self.a
        #Define the angle and critical angles determining if theres
        #Surface area exposure loss for hemisphere glass PMTs
        rdotPMTs = np.array(np.sum(r*pd,axis=1))/(r_mags*pd_mags)
        #There's rounding errors that can make some values greater than 1. 
        GT1 = np.where(rdotPMTs>1)[0]
        rdotPMTs[GT1] = 1.0 
        thetas = np.arccos(rdotPMTs)
        theta_crit_lows = np.arccos(np.sqrt(1 - (a**2/r_mags**2)))
        theta_crit_highs = np.pi - theta_crit_lows
        #We've gone through the calculations for the exposed
        #Solid angle based on theta and the theta_crits. Implement
        exp_solid_angles = np.zeros(len(PMTPos))
        #Now, we fill in the solid angles at each entry based on thetas
        crit_low_indices = np.where(thetas < theta_crit_lows)[0]
        cl_thetas = thetas[crit_low_indices] 
        cl_solid_angles = 2.0*np.pi*(1-np.sqrt(1-(a[crit_low_indices]**2/
            r_mags[crit_low_indices]**2)))
        exp_solid_angles[crit_low_indices] = cl_solid_angles
        
        crit_mid_a = np.where(thetas >= theta_crit_lows)[0]
        crit_mid_b = np.where(thetas <= theta_crit_highs)[0]
        crit_mid_indices = np.intersect1d(crit_mid_a,crit_mid_b)
        cm_thetas = thetas[crit_mid_indices]
        cm_alphas = theta_crit_lows[crit_mid_indices]
        ang_factors = ((np.cos(cm_thetas)/np.cos(cm_alphas))+1)
        cm_solid_angles = np.pi*(1-np.sqrt(1-(a[crit_mid_indices]**2/
                r_mags[crit_mid_indices]**2)))*ang_factors

        exp_solid_angles[crit_mid_indices] = cm_solid_angles 
        #NOTE: Exposed solid angles above theta_crit_high are zero 
        return exp_solid_angles

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

