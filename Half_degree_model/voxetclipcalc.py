# This script will read in the FEILDS Australia GoCAD voxet files and output CSV text files for a selected region, defined by latitude, longitude and h. 
# Options exist also to compute the "residual" gravity field from density variations within the selected region and also the "regional" gravity field 
# from the model (excluding the selected region)

# It requires a functional installation of escript, which can be downloaded here(https://launchpad.net/escript-finley). 
#Some functions (e.g. Silo output and MPI support) will require the appropriate functionality to be enabled in escript  

# Please note this script is intended ONLY for the FEILDS voxet files, and are resolution specific.
# Please note we also assume the directory structure and directory names are not modified  
# Any use apart from these specific voxet files will require significant modifications, which are up to you to implement correctly.  

# When using the FEILDS Australia models or this script please reference the following article:

# Aitken, A.R.A., Altinay, C. and Gross, L. 2016. Australia's lithospheric density field, and its isostatic equilibration. Geophysical Journal International

#Begin script
#import the required escript and python modules 
import os
import sys
import numpy as np
from esys.ripley import *
from time import time
from esys.escript import ReducedFunction
from esys.downunder import *
from esys.escript import unitsSI as U
from esys.escript import *
from esys.escript.linearPDEs import SolverOptions
from esys.ripley import *
from esys.weipa import *

#Define variables to be used in the clipping process, these should be changed

#Clipping Parameters

# Selection in the escript domain is inclusive of the limit value and uses nodal values to define the region
# i.e. only nodes within the selected region, or equal to the limit value, will be selected
# The CSV output will output cell-centroid registered values for cells that have all 8 nodes selected.
# To include the entire FEILDS model use values outside of the domain extents.

# Western, eastern, northern and southern limits of the desired region in decimal degrees.
Longitude_W= 100. # longitude E
Longitude_E= 165. # longitude E
Latitude_N= 5.  # latitude, negative for southern hemisphere
Latitude_S= -55.  # latitude, negative for southern hemisphere

#The top and base of the desired region in kilometres below top of crust (i.e. depth). It is permissible to include some of the air-layer.	
h_top= 0.
h_base= 305.

#Choose the model property (one only) you wish to extract.  If not selecting from the list below, the string must equal 
#the property listed in the voxet header file in every way.

MODEL_PROPERTY = 'True_density'
#MODEL_PROPERTY = 'Density_change_m'
#MODEL_PROPERTY = 'Pressure'
#MODEL_PROPERTY = 'Pressure_Anomaly'

#Gravity calculation parameters

#Set density options for masked region. If UseMean is True the mean density of the selected region will be used. Or if False you can specify a preferred value
UseMean = True
BackgroundDensity = 2670. # in kg/m^3

# Set as True to calculate the gravity effect of the selected region
CalcResidual = False 

# Set as True to calculate the gravity effect of the model excluding the 
# selected region
CalcRegional = False 

# Set the altitude of gravity calculation in meters above the top of the crust. 
# Negative values are acceptable if that's really what you want, but the upper limit is +55000m (the top of the escript domain) 
# Data will be extracted in the centre of the cell in which the data exists (i.e 4000m will extract data at 2500m with 5 km cells)

Altitude = 4000.

#Advanced Options

#MPI support enabled for faster running on computers with MPI support. MPI must be enabled in escript installation
MPIsupport = False
#Allows output of certain properties for the ENTIRE escript domain in Silo database format. This is good for visualisation in, for example, VisIT. 
#Silo support must be enabled in escript installation
SiloOutput = True

#reference to Data, Model and Output directories, only change if you have modified the directory structure
DATADIR='Data'
MODELDIR='Model'
OUTPUTDIR='Output'

# 3d Voxet dataset, only change if you have renamed the voxet header file
MODEL_DATASET = os.path.join(MODELDIR,'Half_degree_model.vo')

#Specify gravity units for output, default is milligals
GRAV_UNITS = 1. * U.mgal

#Reverse selection. This will allow one to output everything APART from the defined region as CSV. Calculations will be valid still, but reversed.
ReverseSelection = False

#Best not to change anything below this line unless you know exactly what it does

#-------------------------------------------------------------------------------------------------------------------------------------

#make the output directory if not existing already
try:
    os.mkdir(OUTPUTDIR)
except:
    pass

#Error reporting for people who make logical mistakes with input coordinates

if Longitude_E <= Longitude_W:
    raise ValueError("Longitude West was >= than longitude East")
if Latitude_N <= Latitude_S:
    raise ValueError("Latitude N was <= Latitude S. Remember negative values for southern hemisphere")
if h_top >= h_base:
    raise ValueError("h_top was >= h_base, Remember h is kilometres below the top of the crust")	

# This subroutine reads in the Voxet
def readVoxet(domain, filename, voproperty=0, origin=None, fillValue=0., reference_system=None):
    header=open(filename).readlines()
    if not header[0].startswith('GOCAD Voxet'):
        raise ValueError("Invalid Voxet file")
    NE=None
    axis_uvw=[None,None,None]
    axis_min=[0.,0.,0.]
    axis_max=[1.,1.,1.]
    # props[id]=[name,file,datatype]
    props={}
    for line in header:
        if line.startswith('AXIS_O '):
            if origin is None:
                origin=[float(i) for i in line.split()[1:4]]
        elif line.startswith('AXIS_U '):
            u=[float(i) for i in line.split()[1:4]]
            if (u[1] != 0) or (u[2] != 0):
                raise ValueError('This coordinate system is not supported')
            axis_uvw[0]=u[0]
        elif line.startswith('AXIS_V '):
            v=[float(i) for i in line.split()[1:4]]
            if (v[0] != 0) or (v[2] != 0):
                raise ValueError('This coordinate system is not supported')
            axis_uvw[1]=v[1]
        elif line.startswith('AXIS_W '):
            w=[float(i) for i in line.split()[1:4]]
            if (w[0] != 0) or (w[1] != 0):
                raise ValueError('This coordinate system is not supported')
            axis_uvw[2]=w[2]
        elif line.startswith('AXIS_MIN '):
            axis_min=[float(i) for i in line.split()[1:4]]
        elif line.startswith('AXIS_MAX '):
            axis_max=[float(i) for i in line.split()[1:4]]
        elif line.startswith('AXIS_N '):
            NE=[int(i) for i in line.split()[1:4]]
        elif line.startswith('PROPERTY '):
            propid=int(line.split()[1])
            if not props.has_key(propid):
                props[propid]=[None,None,None]
            props[propid][0]=line.split()[2].strip()
        elif line.startswith('PROP_ESIZE '):
            propid=int(line.split()[1])
            t=int(line.split()[2])
            if t==4:
                props[propid][2]=DATATYPE_FLOAT32
            elif t==8:
                props[propid][2]=DATATYPE_FLOAT64
            else:
                raise ValueError('Unsupported data size '+t)
        elif line.startswith('PROP_ETYPE '):
            t=line.split()[2].strip()
            if t != 'IEEE':
                raise ValueError('Unsupported data type '+t)
        elif line.startswith('PROP_FORMAT '):
            t=line.split()[2].strip()
            if t != 'RAW':
                raise ValueError('Unsupported data format '+t)
        elif line.startswith('PROP_OFFSET '):
            dataoffset=int(line.split()[2])
            if dataoffset != 0:
                raise ValueError('data offset != 0 not supported yet')
        elif line.startswith('PROP_FILE '):
            propid=int(line.split()[1])
            props[propid][1]=line.split()[2].strip()

    if (axis_uvw[0] is None) or (axis_uvw[1] is None) or (axis_uvw[2] is None)\
            or (NE is None) or (origin is None):
        raise ValueError('Could not determine data configuration. Invalid file?!')
    if len(props)==0:
        raise ValueError('No properties found.')

    #l=[(origin[i], origin[i]+axis_uvw[i]*axis_max[i]) for i in range(3)]
    gridorigin, gridspacing, gridNE = domain.getGridParameters()
    if not reference_system:
         reference_system = CartesianReferenceSystem()

    if reference_system.isCartesian():
        v_scale=1.
    else:
        v_scale=1./reference_system.getHeightUnit()

    # hack to add 2x half cells (top and bottom)
    dz = axis_uvw[-1] / (NE[-1]-1)
    axis_uvw = [x + dz for x in axis_uvw]
    #print axis_uvw
    origin[-1] = origin[-1]*v_scale

    # determine base location of this dataset within the domain
    first=[int((origin[i]-gridorigin[i])/gridspacing[i]) for i in range(domain.getDim())]

    # determine the resolution difference between domain and data. If domain has twice the resolution we can double up the data etc.
	#Not used in this example
    #multiplier=[int(round((abs(axis_uvw[i])/NE[i])/gridspacing[i])) for i in range(domain.getDim())]
    multiplier = [1,1,1]

    datatype=None
    propfile=None
    for p in props.values():
        if (isinstance(voproperty, int) and p[1] == voproperty) or \
           (isinstance(voproperty, str) and p[0] == voproperty):
            datatype=p[2]
            name=p[1]
            #remove quotes which GoCAD introduces for filenames with spaces
            if name.startswith('"') and name.endswith('"'):
                name=name[1:-1]
            propfile=os.path.join(os.path.dirname(filename), name)
            print("Voxet property file: %s"%propfile)
            break

    if propfile is None or datatype is None:
        raise ValueError("Invalid property "+str(voproperty))

    reverse=[0]*domain.getDim()
    if axis_uvw[-1] < 0:
        reverse[-1]=1

    print("calling readBinaryGrid with first=%s, nValues=%s, multiplier=%s, reverse=%s"%(str(first),str(NE),str(multiplier),str(reverse)))
    data=readBinaryGrid(propfile, ReducedFunction(domain), shape=(),
            fill=fillValue, byteOrder=BYTEORDER_BIG_ENDIAN,
            dataType=p[2], first=first, numValues=NE, multiplier=multiplier,
            reverse=reverse)

    return data

# 3d domain properties

resolution = 0.5       # grid cell size (we are assuming equal in x and y)
originX = 100.0         # longitude origin (cell centre!)
originY = -55.0         # latitude origin (cell centre!)
thickness = 305. * U.km # below surface
l_air = 55. * U.km      # above surface
NX = 131               # number of grid cells in longitude
NY = 121                # number of grid cells in latitude
n_cells_v = 72         # total number of cells in the vertical direction (inc air)
DZ=(thickness+l_air)/n_cells_v #cell height
cell_at_altitude = int(np.floor(thickness+Altitude/1000)/DZ) #Figure out which cell to calculate the data within
CaCentre = DZ/2-thickness+cell_at_altitude*DZ # Only used for reporting purposes
print ('cell at data altitude is #%s which has a central height of %s'%(cell_at_altitude, CaCentre))

# shift origin because data is cell centred and we need "real" grid origin
originX = originX - resolution/2.
originY = originY - resolution/2.

# reference model origin (west, south, bottom)
MODEL_ORIGIN = [originX, originY, -thickness]

# we are using geodetic coordinates
COORDINATES=WGS84ReferenceSystem()

h_unit = COORDINATES.getHeightUnit()
dom_origin=[originX, originY, -thickness/h_unit]
dom_origin=[1e-5*np.floor(oi*1e5) for oi in dom_origin]
spacing=[resolution, resolution, np.floor((thickness+l_air)/n_cells_v)/h_unit]
l0=(dom_origin[0], dom_origin[0]+spacing[0]*NX)
l1=(dom_origin[1], dom_origin[1]+spacing[1]*NY)
l2=(dom_origin[2], dom_origin[2]+spacing[2]*n_cells_v)
d0,d1,d2=1,1,1

if MPIsupport is True:
     # ensure ripley does not add additional elements in horizontal direction
    
     mpifactors=list(factor(d2))
     while len(mpifactors)>0:
         f=mpifactors.pop()
         if d1>=d0:
            if (NX+1)%(d0*f) == 0:
                 d0 = d0 * f
                 d2 = d2 / f
            elif (NY+1)%(d1*f) == 0:
                 d1 = d1 * f
                 d2 = d2 / f
         else:
            if (NY+1)%(d1*f) == 0:
                 d1 = d1 * f
                 d2 = d2 / f
            elif (NX+1)%(d0*f) == 0:
                 d0 = d0 * f
                 d2 = d2 / f

# create domain
print("Domain subdivisions: %d x %d x %d"%(d0,d1,d2))
dom=Brick(NX, NY, n_cells_v, l0, l1, l2, d0, d1, d2)

dom_len=[sup(dom.getX()[i])-inf(dom.getX()[i]) for i in range(dom.getDim())]

# report domain setup
print("Domain size: "+str([NX, NY, n_cells_v]))
print("     length: "+str(dom_len))
print("     origin: "+str(dom_origin))

DIM=dom.getDim()

datacoords = ReducedFunction(dom).getX()

#read in the FEILDS Voxet
initial_model = readVoxet(dom, MODEL_DATASET, MODEL_PROPERTY, MODEL_ORIGIN[:],
                        0., reference_system=COORDINATES)


# set the extents of the desired "region" for clipping and computation, mask is = 1 OUTSIDE area of interest

mask_air = whereNonNegative(dom.getX()[2]+spacing[2]/2) #reg_mask for air layer

mask_LONG = wherePositive(Longitude_W-dom.getX()[0]) + whereNegative(Longitude_E-dom.getX()[0]) #reg_mask for longitude
mask_LAT = whereNegative(Latitude_N-dom.getX()[1]) + wherePositive(Latitude_S-dom.getX()[1])    #reg_mask for latitude
mask_h = wherePositive(dom.getX()[2]+(h_top/100)) + whereNegative(dom.getX()[2]+(h_base/100))   #reg_mask for depth

reg_mask = whereZero(mask_LONG+mask_LAT+mask_h)

if ReverseSelection is True:
     reg_mask = whereNonZero(mask_LONG+mask_LAT+mask_h)

	 
#prior to any computation, write out the selected region model as CSV format
fn=os.path.join(OUTPUTDIR,"region_%s")%(MODEL_PROPERTY)
#OK this doesn't seem to be working with the mask....perhaps I need a reduced function somewhere (mask? or initial_model?).
saveDataCSV(fn+".csv", Long=datacoords[0], Lat=datacoords[1], h=datacoords[2], PROPERTY=initial_model, masking=reg_mask, mask=reg_mask)
print('CSV file written with the following fields: Longitude (degrees) Latitude (degrees), h (100km), Property (kg/m^3 or Pa)')

if SiloOutput is True:
     saveSilo(fn, PROPERTY=initial_model, mask=reg_mask) #mask is fine in here - do we need to define prop
     print('SILO file written with the following fields: Property (kg/m^3 or Pa), mask')

#calculates the "residual" from the selected region for the entire FEILDS model region and outputs gravity at the specified altitude...see below for the "regional"
def ResidualCalculation(): 
        
    #read in a gravity data grid to define data computation - not really active - perhaps can be omitted?
    G_DATA = os.path.join(DATADIR,'Final_BouguerTC_UC25K_halfdeg')
    FS=ReducedFunction(dom)
    nValues=[NX, NY, 1]
    first = [0, 0, cell_at_altitude]
    multiplier = [1, 1, 1]
    reverse = [0, 0, 0]
    byteorder = BYTEORDER_NATIVE
    print("calling readBinaryGrid with first=%s, nValues=%s, multiplier=%s, reverse=%s"%(str(first),str(nValues),str(multiplier),str(reverse)))

    gdata = readBinaryGrid(G_DATA, FS, shape=(),
                fill=-999999, byteOrder=byteorder,
                dataType=DATATYPE_FLOAT32, first=first, numValues=nValues,
                multiplier=multiplier, reverse=reverse)
    print("grid successfully read")
        
    #get the masking and units sorted out for the data-space
    g_mask = whereNonZero(gdata+999999)

    gdata=gdata*g_mask * GRAV_UNITS 	
        
    global reg_mask  
        
    # if for some reason people choose to have air in their region we exclude it from the gravity calculation region
    if h_top < 0.:
	     reg_mask = reg_mask+mask_air
         
    live_model = initial_model* wherePositive(reg_mask)
		        
    global BackgroundDensity
		
	if UseMean is True:
         #calculate the mean density within the selected region
         BackgroundDensity = integrate(live_model)/integrate(wherePositive(reg_mask)) 
	print("Density mean for selected region equals = %s"%BackgroundDensity)
        
	live_model = live_model + BackgroundDensity * whereNonPositive(reg_mask)   

    # create mapping
    rho_mapping=DensityMapping(dom, rho0=live_model) 

    # create regularization with two level set functions....necessary to run forward model? but not "active" so the weightings don't matter.
	regularization=Regularization(dom, numLevelSets=1,
                           w0=1,              # weighting factor for m**2 term
                           w1=np.ones((DIM,)), # consider gradient terms
                           coordinates=COORDINATES,
                           location_of_set_m=reg_mask)

	# invert sign of gravity field to account for escript's coordinate system
	gdata = -GRAV_UNITS * gdata

	# turn the scalars into vectors (vertical direction)
	d=kronecker(DIM)[DIM-1]
	w=safeDiv(1., g_mask)
	gravity_model=GravityModel(dom, w*d, gdata*d, fixPotentialAtBottom=False, coordinates=COORDINATES)
    gravity_model.rescaleWeights(rho_scale=rho_mapping.getTypicalDerivative())
    cf=InversionCostFunction(regularization, rho_mapping, gravity_model)
    phi,_ = cf.getForwardModel().getArguments(live_model)
    g_init = -regularization.getCoordinateTransformation().getGradient(phi)
    g_init = interpolate(g_init, gdata.getFunctionSpace())
    print("Computed gravity: %s"%(g_init[2]))
        
    fn=os.path.join(OUTPUTDIR,'residual-gravity')
    if SiloOutput is True:
	     saveSilo(fn, density=live_model, gravity_init=g_init, g_initz=-g_init[2], gravitymask=g_mask, modelmask=reg_mask)
         print('SILO file written with the following fields: density (kg/m^3), gravity (m/s^2), gz (m/s^2), gravitymask, modelmask')
        
	# to compare calculated data against input dataset. Not used by default
	#gslice = g_init[2]*wherePositive(g_mask)
    #g_dash = integrate(gslice)/integrate(wherePositive(g_mask))
    #gdataslice = gdata*wherePositive(g_mask)
    #gdata_dash = integrate(gdataslice)/integrate(wherePositive(g_mask))
    #misfit=(gdataslice-gdata_dash)-(gslice-g_dash)
    saveDataCSV(fn+".csv", mask=g_mask, gz=-g_init[2], Long=datacoords[0], Lat=datacoords[1], h=datacoords[2]) 
	print('CSV file written with the following fields: Longitude (degrees) Latitude (degrees), h (100km), gz (m/s^2)')
	# As above but including g-data and misfit
	#saveDataCSV(fn+".csv", mask=g_mask, ginit=-g_init[2], X=datacoords[0], Y=datacoords[1], Z=datacoords[2], gdata=gdata, misfit=misfit)
		
#calculates the "regional" from the entire FEILDS model excluding the selected region and outputs gravity at the specified altitude...see above for the "residual"		
def RegionalCalculation(): 
                
    #read in a gravity data grid to define data computation - not really active - perhaps can be omitted?
    G_DATA = os.path.join(DATADIR,'Final_BouguerTC_UC25K_halfdeg')
    FS=ReducedFunction(dom)
    nValues=[NX, NY, 1]
    first = [0, 0, cell_at_altitude]
    multiplier = [1, 1, 1]
    reverse = [0, 0, 0]
    byteorder = BYTEORDER_NATIVE
    print("calling readBinaryGrid with first=%s, nValues=%s, multiplier=%s, reverse=%s"%(str(first),str(nValues),str(multiplier),str(reverse)))

    gdata = readBinaryGrid(G_DATA, FS, shape=(),
                fill=-999999, byteOrder=byteorder,
                dataType=DATATYPE_FLOAT32, first=first, numValues=nValues,
                multiplier=multiplier, reverse=reverse)
    print("grid successfully read")
        
    #get the masking and units sorted out for the data-space
    g_mask = whereNonZero(gdata+999999)

    gdata=gdata*g_mask * GRAV_UNITS 	
        
    global reg_mask  
        
    # if people choose to have air in their region we exclude it from the specified gravity calculation region
    if h_top < 0.:
	    reg_mask = reg_mask+mask_air
         
    live_model = initial_model* whereNonPositive(reg_mask)
	dead_model = initial_model* wherePositive(reg_mask)        
	
	global BackgroundDensity
		
	if UseMean is True:
        #calculate the mean density within the selected region
        BackgroundDensity = integrate(dead_model)/integrate(wherePositive(reg_mask)) 
	    print("Density mean for selected region equals = %s"%BackgroundDensity)
        
	live_model = live_model + BackgroundDensity * wherePositive(reg_mask)   

    # create mapping
    rho_mapping=DensityMapping(dom, rho0=live_model) 

    # create regularization with two level set functions....necessary to run forward model? but not "active" so the weightings don't matter.
	regularization=Regularization(dom, numLevelSets=1,
                           w0=1,              # weighting factor for m**2 term
                           w1=np.ones((DIM,)), # consider gradient terms
                           coordinates=COORDINATES,
                           location_of_set_m=reg_mask)

	# invert sign of gravity field to account for escript's coordinate system
	gdata = -GRAV_UNITS * gdata

	# turn the scalars into vectors (vertical direction)
	d=kronecker(DIM)[DIM-1]
	w=safeDiv(1., g_mask)
	gravity_model=GravityModel(dom, w*d, gdata*d, fixPotentialAtBottom=False, coordinates=COORDINATES)
    gravity_model.rescaleWeights(rho_scale=rho_mapping.getTypicalDerivative())
    cf=InversionCostFunction(regularization, rho_mapping, gravity_model)
    phi,_ = cf.getForwardModel().getArguments(live_model)
    g_init = -regularization.getCoordinateTransformation().getGradient(phi)
    g_init = interpolate(g_init, gdata.getFunctionSpace())
    print("Computed gravity: %s"%(g_init[2]))
        
    fn=os.path.join(OUTPUTDIR,'regional-gravity')
    if SiloOutput is True:
	    saveSilo(fn, density=live_model, gravity_init=g_init, g_initz=-g_init[2], gravitymask=g_mask, modelmask=reg_mask)
        print('SILO file written with the following fields: density (kg/m^3), gravity vector (m/s^2), gz (m/s^2), gravitymask, modelmask')
        
	# to compare calculated data against input dataset. Not used by default but should work if the input dataset is correct
	#gslice = g_init[2]*wherePositive(g_mask)
    #g_dash = integrate(gslice)/integrate(wherePositive(g_mask))
    #gdataslice = gdata*wherePositive(g_mask)
    #gdata_dash = integrate(gdataslice)/integrate(wherePositive(g_mask))
    #misfit=(gdataslice-gdata_dash)-(gslice-g_dash)
    saveDataCSV(fn+".csv", mask=g_mask, gz=-g_init[2], Long=datacoords[0], Lat=datacoords[1], h=datacoords[2]) 
	print('CSV file written with the following fields: Longitude (degrees) Latitude (degrees), h (100km), gz (m/s^2)')
	# As above but including g-data and misfit
	#saveDataCSV(fn+".csv", mask=g_mask, ginit=-g_init[2], X=datacoords[0], Y=datacoords[1], Z=datacoords[2], gdata=gdata, misfit=misfit)
		
if CalcResidual is True:
    ResidualCalculation()
    print('Residual gravity calculation is complete')
if CalcRegional is True:
    RegionalCalculation()
    print('Regional gravity calculation is complete')
print "Finished, Thanks for using the FEILDS Australia model"

