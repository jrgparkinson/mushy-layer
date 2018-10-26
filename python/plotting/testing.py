import sys
#from visit_utils import *

# To run:
# visit -nowin -cli -s <script.py>

# Make sure noise.silo is the active source

database = "/home/parkinsonjl/convection-in-sea-ice/MushyLayer/execSubcycle/output/fixedChill-periodic-CR20.000RaC10Le100KozenyPermeabilityDa1.0e-02R1.0e-03pts256-St5-domWidth32.0-001024.2d.hdf5"
OpenDatabase(database)


s=SaveWindowAttributes()
s.format=s.JPEG
s.screenCapture = 0
s.resConstraint = s.NoConstraint
s.width = 2048
s.height = 2048
s.outputToCurrentDirectory=0
s.outputDirectory="/home/parkinsonjl/convection-in-sea-ice/MushyLayer/python/plotting/"
s.fileName="test-"
SetSaveWindowAttributes(s)


#ga = GetGlobalAttributes()
#for src in ga.sources:
#    print(str(src))
#    if src.count("noise.silo"):
#        ActivateDatabase(src)
#    else:
#        CloseDatabase(src)


# Clear any previous plots
DeleteAllPlots()
# Create a plot of the scalar field 'temp'
AddPlot("Pseudocolor","Temperature")
AddPlot("Subset", "levels") 
# Slice the volume to show only three
# external faces.
#AddOperator("ThreeSlice")
#tatts = ThreeSliceAttributes()
#tatts.x = -10
#tatts.y = -10
#tatts.z = -10
#SetOperatorOptions(tatts)

AnnotationAtts = AnnotationAttributes()
AnnotationAtts.databaseInfoFlag = 0
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.axes3D.visible = 0
AnnotationAtts.axes3D.triadFlag = 0
AnnotationAtts.axes3D.bboxFlag = 0
AnnotationAtts.legendInfoFlag = 1
AnnotationAtts.backgroundColor = (255, 255, 255, 255)
AnnotationAtts.foregroundColor = (0,0,0,0)
SetAnnotationAttributes(AnnotationAtts)



DrawPlots()



# Find the maximum value of the field 'temp'
Query("Max")
val = GetQueryOutputValue()
print "Max value of 'temp' = ", val






SaveWindow()

DeleteAllPlots()
CloseDatabase(database)

sys.exit()
