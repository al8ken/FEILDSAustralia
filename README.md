# FEILDS_Australia
Results of the FEILDS Australia project at 1/8 1/4 and 1/2 degree resolution.

This project contains the data and models resulting from the FEILDS Australia gravity inversion program.

Use of these data and models should be acknowledged by citation of the following manuscript:

Aitken, A.R.A., Altinay, C. and Gross, L. 2016. Australia's lithospheric density field, and its isostatic equilibration. Geophysical Journal International

Please read this manuscript carefully to understand the uncertainties and biases inherent in this dataset prior to using these models.

# Model Directory Structure and contents

For these models we include both Data and Model directories, as well as a model extraction and calculation tool

#DATA:

This folder has three files:

1) Final_BouguerTC_UC25K_halfdeg.nc 
	NetCDF grid of Bouguer and Terrain corrected gravity data used in the inversion in mGal. UCxK denotes upwards continuation distance
2) Data_Elevation_halfdeg.nc 
	NetCDF grid of Data elevations for the Bouguer Gravity data in metres...NB these are subsampled onto the mesh, so are less vertically precise in the model than in this file. 
3) final-mu0.005-wo1e-07.csv 
	Data in CSV format.  Long - Longitude (deg), Lat - Latitude (deg), h - Height above ellipsoid (100skm), g_data - Input data in m/s^2, g_calc - Calculated Data in m/s^2
	
#MODEL:

This folder has three files (one for very high resolution), and one folder (VoxFiles):

1) Quarter_degree_model.vo
	Gocad Voxet file containing 14 properties of the model. Coordinates are u=Latitude and v=longitude (in degrees) and w=height above ellipsoid (in 100s km)
	Individual property binaries are in the VoxFiles directory. These are:
	Density Parameters (all in kg/m^3)
		Density --- Density returned from the inversion
		True_density --- is Density plus the mean density removed prior to inversion (Half degree: 3421.461; Quarter degree: 3422.180)
		Density_initial --- the initial density model
		Density_change_m --- m, the density correction resulting from inversion
		Mu0Range_m --- Variability in m as a result of Mu0
		Mu1Range_m --- Variability in m as a result of Mu1
	Gravity Parameters (all in m/s^2)
		gx --- calculated anomalous gravitational acceleration in the meridional direction (+ve indicates eastwards acceleration)
		gy --- calculated anomalous gravitational acceleration in the latitudinal direction (+ve indicates northwards acceleration)
		gz --- calculated anomalous gravitational acceleration in the vertical direction (+ve indicates upwards acceleration)
	Pressure Parameters (all in Pa)	
		p0 --- Topographic pressure boundary condition
		Pressure --- Total Pressure returned from the model
		Pb_at_PseudoElevation --- Background pressure calculated for PseudoElevations
		Pressure_Anomaly --- is Pressure - Pb_at_PseudoElevation
	    PseudoElevation --- PseudoElevation in 100skm, see text for details of derivation
2) Quarter_degree_model_density.xyz
	The entire density model realised in space-delimited ASCII format. This is the density ANOMALY, i.e. it is necessary to add the mean density for real-valued data. Format: 1 header line, columns are Longitude(degrees), Latitude(degrees), h(100s km), density(kg/m^3). Due to size restrictions, this file is NOT provided for the very-high resolution model.
3) Quarter_degree_model_pressure_anomaly.xyz	
	The entire pressure model realised in space-delimited ASCII format. This is the pressure ANOMALY, i.e. it is necessary to add the background pressure for real-valued data. Format: 1 header line, columns are Longitude(degrees), Latitude(degrees), h(100s km), pressure anomaly(Pa). Due to size restrictions, this file is NOT provided for the very-high resolution model.

In the Model_Grids directory, the very high resolution model also includes raster grids (in ERMapper format) of True_density at the same depths as the initial AuSREM model.  

# Voxet extraction and calculation tool

For ease of use we also supply for each model resolution a voxet extraction and calculation tool in the script voxetclipcalc.py. This script requires escript to be installed. This has three main functions:

1 - Extraction of model sub-volumes as CSV format (e.g. to extract only a desired Lat/Long and/or depth-range or to extract all but a selected volume )
2 - Calculation of the regional field excluding a selected model volume (e.g. as a precursor to more detailed modelling of that volume) 
3 - Calculation of the residual field from a selected model volume (e.g. as a precursor to more detailed modelling of that volume)

Please see comments in the script for further details of this tool. Note that, as we implement it here, it is model specific and relies on unchanged filenames and directory structure.
=======
# FEILDSAustralia
Data and results of the Finite Element Inversion of Lithospheric Density Structure (FEILDS) project for Australia at 1/8 1/4 and 1/2 degree resolution.
