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
import json

if __name__=='__main__':
    print("########### WELCOME TO WATCH SHAPES ############")
    print("###### LET'S HAVE SOME FUN WITH GEOMETRY #######")
    ####### BEGIN MODIFIERS #######
    #PMT Specs 
    NUMPMTS = 4400
    PR = 127.0   #PMT RADIUS in millimeters
    CENTRAL_POINTING = False 
    #Detector Specs 
    geometry = "SPHERE"
    DEFINE_VIA_BUFFER = True #if true, RADUIS & HEIGHT defns overwritten with
                             #the FV + BUFFER dimension
    BUFFER_VOLUME = 2000.0 #Desired buffer volume 
    TONS_WATER_FV = 1000.0 #Tons of water in fiducial volume
    ATTENUATION_LENGTH = 1.0E9  #Attenuation length of medium in millimeters 
    RADIUS = 8000.0 #Radius of detector (sphere or cylinder) in millimeters
    HEIGHT = 16000.0 #Height of detector (cylinder only) in millimeters
    #Results Saved specs 
    OUTFILE = 'results_%s.json'%(geometry) 
    NUM_FVSAMPLES = 1000
    
    ######## END MODIFIERS #########
    WATCHMAN = d.DetectorWorkshop(geometry=geometry)
    WATCHMAN.SetFiducialMassAndDimensions(TONS_WATER_FV)
    if DEFINE_VIA_BUFFER is True:
        RADIUS = WATCHMAN.GetDetectorRadius(BUFFER_VOLUME)
    if geometry == "SPHERE":
        pointGen = pp.SpherePoints(numpoints = NUMPMTS, radius=RADIUS)
        pointGen.PopulatePositions()
        pointGen.ShowPositions()
    if geometry == "CYLINDER":
        if DEFINE_VIA_BUFFER is True: 
            HEIGHT = WATCHMAN.GetDetectorHeight(BUFFER_VOLUME)
        pointGen = pp.CylinderPoints(numpoints = NUMPMTS, radius=RADIUS,height=HEIGHT)
        pointGen.PopulatePositions()
        pointGen.ShowPositions()
    PMTInformation = pmt.PMTGen(PointGeometry=pointGen,PMTRadius=PR)
    PMTInformation.BuildPMTDirections(face_center=CENTRAL_POINTING)
    
    WATCHMAN.LoadPMTInfo(PMTInformation) 
    WATCHMAN.SetRadius(RADIUS)
    WATCHMAN.SetAttenuationLength(ATTENUATION_LENGTH)
    if geometry == "CYLINDER":
        WATCHMAN.SetHeight(HEIGHT)
    fv_positions = []
    fv_positions_outsave = []
    light_factors = []
    for j in xrange(NUM_FVSAMPLES):
        FVposition = WATCHMAN.ShootFVPosition(axis='z')
        fv_positions.append(FVposition)
        fv_positions_outsave.append(list(FVposition))
        light_factors.append(WATCHMAN.EvaluateLightCollection(fv_positions[j]))
    results_dict = {'detector_specs':WATCHMAN.shape, 'FV_specs':WATCHMAN.fiducial_shape,
            'FV_positions': fv_positions_outsave, 'light_factors': light_factors}
    with open(OUTFILE,"w") as f:
        json.dump(results_dict, f, sort_keys=True, indent=4)
