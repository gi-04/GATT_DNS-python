# This script is the main preprocessing routine, which will generate the mesh, prepare the boundaries, the matrices and write the Fortran files

# Required python packages
import os
import scipy.io as sio
# Additional functions and files
from __main__ import *
from meshAddFixedPoints import *

# Generate mesh
# Add fixed points to mesh structure
# The mesh will be slightly transformed so that key points are present at
# the correct positions if mesh.x.matchFixed = true
# X = 0 and Y = 0 are always added
meshAddFixedPoints()

# Run mesh generator
mesh['X'], mesh['x'], mesh['nx'] = generateMesh(domain['xi'], domain['xf'], mesh['x'], 'X')
mesh['Y'], mesh['y'], mesh['ny'] = generateMesh(domain['yi'], domain['yf'], mesh['y'], 'Y')
mesh['Z'], mesh['z'], mesh['nz'] = generateMesh(domain['zi'], domain['zf'], mesh['z'], 'Z')

# Select boundary conditions
boundary, mesh = getBoundaryConditions(flowType, mesh, flowParameters, [numMethods['neumannOrder'], numMethods['neumann2Order']])

domainSlicesY = getDomainSlices(mesh['ny'], p_row)
domainSlicesZ = getDomainSlices(mesh['nz'], p_col)

boundaryInfo = initBoundaries(boundary, mesh, domainSlicesY, domainSlicesZ, p_row, p_col)

# Make matrices for spatial derivatives and for spatial filters
matrices = makeMatrices(mesh, domain, boundary, numMethods)

matrices['x'] = prepareThomas(matrices['x'])
matrices['x']['blocks'] = getMatrixTypeBlocks(matrices['x']['types'], p_row, p_col)
matrices['y'] = prepareThomas(matrices['y'])
matrices['y']['blocks'] = getMatrixTypeBlocks(matrices['y']['types'], p_row, p_col)
if mesh['nz'] > 1:
    matrices['z'] = prepareThomas(matrices['z'])
    matrices['z']['blocks'] = getMatrixTypeBlocks(matrices['z']['types'], p_row, p_col)

matrices['neumannCoeffs'] = boundary['neumannCoeffs']
matrices['neumann2Coeffs'] = boundary['neumann2Coeffs']

# Write files
if not os.path.exists(caseName):
    os.makedirs(caseName)

# Prepare SFD
if 'SFD' in numMethods:
    if numMethods['SFD']['type'] == 2:
        calcSFDregion()
    if numMethods['SFD']['Delta'] == float('inf'):
        numMethods['SFD']['Delta'] = -1
    if numMethods['SFD']['type'] > 0 and (os.path.exists(f"{caseName}/meanflowSFD.mat") or 'meanFile' in flowType['initial']):
        if 'resume' not in numMethods['SFD']:
            numMethods['SFD']['resume'] = 1
    else:
        numMethods['SFD']['resume'] = 0

# Mesh file
X = mesh['X']
Y = mesh['Y']
Z = mesh['Z']
wall = boundary['insideWall']
sio.savemat(f"{caseName}/mesh.mat", {'X': X, 'Y': Y, 'Z': Z, 'wall': wall, 'flowParameters': flowParameters, 'flowType': flowType})

# Check for previous save files
if 'runningLST' not in globals():
    nStep, nx, ny, nz = checkPreviousRun(caseName)
    if nStep:
        genInitialFlow = False
        time['nStep'] = nStep
        if (nx != mesh['nx'] or ny != mesh['ny'] or nz != mesh['nz']):
            raise ValueError(f"Mesh size has changed since last run: Parameters file indicates {mesh['nx']}x{mesh['ny']}x{mesh['nz']} but mat file contains data for {nx}x{ny}x{nz}")
    else:
        genInitialFlow = True
        time['nStep'] = 0
else:
    genInitialFlow = False
    time['nStep'] = 0
    numMethods['SFD']['resume'] = 0

# Fortran files
if not os.path.exists(f"{caseName}/bin"):
    os.makedirs(f"{caseName}/bin")

disturbTypes = writeFortranDisturbances(caseName, boundaryInfo, tridimensional)

writeFortranParameters(caseName, mesh, flowParameters, time, numMethods, logAll, p_row, p_col)

writeFortranMatrices(caseName, matrices, numMethods, mesh)

writeFortranBoundaries(caseName, boundaryInfo)

# SFD
if numMethods['SFD']['type'] == 2:
    sio.savemat(f"{caseName}/bin/SFD.mat", {'SFD_X': SFD_X})
