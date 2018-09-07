#So we'll have the main code here.

#ideas for structure:
#We have the point geometry class and PMT info class
#We now need two classes I think:
#   1. A class that takes in PMTInformation, and will also have the
#      Fiducial volume information and medium information.
#      Let's call it the Detector class.
#   2. The 
import lib.PointGeometries as pp
import lib.PMTPosition as pmt
import lib.Detector as d
import numpy as np
import json

if __name__=='__main__':
    print("########### WELCOME TO WATCH SHAPES ############")
    print("###### LET'S HAVE SOME FUN WITH GEOMETRY #######")
    
    
    ####### BEGIN MODIFIERS #######
    #Use point segment monte carlo?  If False, hemispheres used
    USESEG = True
    frac_exposed = 0.5  #Hemisphere cathode; should match analytical
    #PMT Specs 
    NUMPMTS = 4800
    PR = 127.0   #PMT RADIUS in millimeters
    CENTRAL_POINTING = False 
    #Detector Specs 
    geometry = "SPHERE"
    DEFINE_DETECTOR_VIA_BUFFER = True #if true, RADUIS & HEIGHT defns overwritten with
                                      #the FV + BUFFER dimension
    BUFFER_STANDOFF = 2000.0 #Desired buffer volume 
    TONS_FV = 1000.0 #Tons of medium in fiducial volume
    DENSITY = 1.0 #Density of detector medium in kg/m**3
    ATTENUATION_LENGTH = 80000.0  #Attenuation length of medium in millimeters 
    RADIUS = 8000.0 #Radius of detector (sphere or cylinder) in millimeters
    HEIGHT = 16000.0 #Height of detector (cylinder only) in millimeters
    #Results Saved specs 
    OUTFILE = 'results_%s.json'%(geometry) 
    if USESEG is True:
        OUTFILE = 'results_%s_useseg.json'%(geometry)
    NUM_SAMPLES = 10000
    ######## END MODIFIERS #########

    
    WATCHMAN = d.DetectorWorkshop(geometry=geometry)
    WATCHMAN.SetMediumDensity(DENSITY) 
    WATCHMAN.SetAttenuationLength(ATTENUATION_LENGTH)
    WATCHMAN.SetFiducialMassAndDimensions(TONS_FV)
    if DEFINE_DETECTOR_VIA_BUFFER is True:
        RADIUS = WATCHMAN.GetDetectorRadius(BUFFER_STANDOFF)
    if geometry == "SPHERE":
        pointGen = pp.SpherePoints(numpoints = NUMPMTS, radius=RADIUS)
        pointGen.PopulatePositions()
        pointGen.ShowPositions()
    if geometry == "CYLINDER":
        if DEFINE_DETECTOR_VIA_BUFFER is True: 
            HEIGHT = WATCHMAN.GetDetectorHeight(BUFFER_STANDOFF)
        pointGen = pp.CylinderPoints(numpoints = NUMPMTS, radius=RADIUS,height=HEIGHT)
        pointGen.PopulatePositions()
        pointGen.ShowPositions()
    
    PMTInformation = None
    if USESEG is False:
        PMTInformation = pmt.PMTGen(PointGeometry=pointGen,PMTRadius=PR)
        PMTInformation.BuildPMTDirections(face_center=CENTRAL_POINTING)

    elif USESEG is True:
        PMTInformation = pmt.SegSurfaceGen(PointGeometry=pointGen,PMTRadius=PR,
                frac_exposed = frac_exposed)
        PMTInformation.BuildPMTDirections(face_center=CENTRAL_POINTING)
        PMTInformation.BuildPMTSurfacePoints() 
    
    WATCHMAN.LoadPMTInfo(PMTInformation) 
    WATCHMAN.SetRadius(RADIUS)
    if geometry == "CYLINDER":
        WATCHMAN.SetHeight(HEIGHT)
    fv_positions = []
    fv_positions_outsave = []
    light_factors = []
    for j in xrange(NUM_SAMPLES):
        FVposition = WATCHMAN.ShootPosition(axis='x', region='FV')
        fv_positions.append(FVposition)
        fv_positions_outsave.append(list(FVposition))
        light_factors.append(WATCHMAN.EvaluateLightCollection(fv_positions[j],UseSegSurface=USESEG))
    results_dict = {'use_seg_mc':USESEG, 'detector_specs':WATCHMAN.shape, 'FV_specs':WATCHMAN.fiducial_shape,
            'FV_positions': fv_positions_outsave, 'light_factors': light_factors,
            'PMT_radius': WATCHMAN.PMTs.PMTRadius, 'Attenuation_coeff':WATCHMAN.attncf,
            'Medium_properties':WATCHMAN.medium}
    with open(OUTFILE,"w") as f:
        json.dump(results_dict, f, sort_keys=True, indent=4)
