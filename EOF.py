# Usage: python EOF.py
# Input: Reads in a text file with extension ".VSGR"
# Input: Format of input file as follows
# Input: STL start
# Input: <STL Filenames> <levels of refinement of each | optional>
# Input: STL end
# Input: Refinement start
# Input: <STL Filenames around which refinement needed> <mode | optional> <levels | optional>
# Input: Refinement end
# Input: If <levels of refinement of each> not specified, then defaults to (2 3)
# Input: If <mode> not specified, defaults to "inside"
# Input: If <levels> not specified, defaults to (1.00 3)
# Output: OpenFOAM files created in directory from which EOF.py was called
# Notes: blockMeshDict created such that there are 10 cells in Z direction. 
# Notes: Number of cells in X & Y computed to ensure aspect ratio as close to 1 as possible
# Created by: Venugopalan Raghavan

## Setup basic stuff

header = list()
header.append("/*------------------------------------------------------------------------*\\\n")
header.append("|=========                 |                                               |\n")
header.append("|\\\\      /   F ield        | OpenFOAM: The Open Source CFD Toolbox         |\n")
header.append("| \\\\    /    O peration    | Version:  2.1.0                               |\n")
header.append("|  \\\\  /     A nd          | Web:      www.OpenFOAM.org                    |\n")
header.append("|   \\\\/      M anipulation |                                               |\n")
header.append("\*------------------------------------------------------------------------*/\n")
header.append("\n")
header.append("/* Exported using EOFv */\n")
header.append("\n")
header.append("FoamFile\n")
header.append("{\n")
header.append("\tversion\t2.1;\n")
header.append("\tformat\tascii;\n")

blockMesh = ["\tclass\tdictionary;\n","\tlocation\tconstant;\n","\tobject\tblockMeshDict;\n","}\n\n","convertToMeters\t1;\n"]

snappyMesh = ["\tclass\tdictionary;\n","\tlocation\tsystem;\n","\tobject\tsnappyHexMeshDict;\n","}\n\n","\ncastellatedMesh\ttrue;","\nsnap\ttrue;","\naddLayers\tfalse;"]

snappyAdd = list()
snappyAdd.append("\nsnapControls")
snappyAdd.append("\n{")
snappyAdd.append("\n\tnSmoothPatch\t2;")
snappyAdd.append("\n\ttolerance\t4;")
snappyAdd.append("\n\tnSolveIter\t20;")
snappyAdd.append("\n\tnRelaxIter\t4;")
snappyAdd.append("\n\tnFeatureSnapIter\t10;")
snappyAdd.append("\n}")
snappyAdd.append("\naddLayersControls")
snappyAdd.append("\n{")
snappyAdd.append("\n\trelativeSizes\ttrue;")
snappyAdd.append("\n\texpansionRatio\t1.3;")
snappyAdd.append("\n\tfinalLayerThickness\t0.4;")
snappyAdd.append("\n\tminThickness\t0.3;")
snappyAdd.append("\n\tnGrow\t0;")
snappyAdd.append("\n\tfeatureAngle\t45;")
snappyAdd.append("\n\tnRelaxIter\t4;")
snappyAdd.append("\n\tnSmoothSurfaceNormals\t1;")
snappyAdd.append("\n\tnSmoothNormals\t1;")
snappyAdd.append("\n\tnSmoothThickness\t10;")
snappyAdd.append("\n\tmaxFaceThicknessRatio\t0.4;")
snappyAdd.append("\n\tmaxThicknesstoMedialRatio\t4;")
snappyAdd.append("\n\tminMedianAxisAngle\t130;")
snappyAdd.append("\n\tnBufferCellsNoExtrude\t0;")
snappyAdd.append("\n\tnLayerIter\t30;")
snappyAdd.append("\n\tlayers")
snappyAdd.append("\n\t{")

snappyQual = list()
snappyQual.append("\n\nmeshQualityControls")
snappyQual.append("\n{")
snappyQual.append("\n\tmaxNonOrtho\t65;")
snappyQual.append("\n\tmaxBoundarySkewness\t20;")
snappyQual.append("\n\tmaxInternalSkewness\t4;")
snappyQual.append("\n\tmaxConcave\t80;")
snappyQual.append("\n\tminFlatness\t0.5;")
snappyQual.append("\n\tminVol\t1E-13;")
snappyQual.append("\n\tminArea\t1E-13;")
snappyQual.append("\n\tminTwist\t0.05;")
snappyQual.append("\n\tminDeterminant\t0.001;")
snappyQual.append("\n\tminFaceWeight\t0.06;")
snappyQual.append("\n\tminVolRatio\t0.025;")
snappyQual.append("\n\tminTriangleTwist\t-0.99;")
snappyQual.append("\n\tnSmoothScale\t4;")
snappyQual.append("\n\terrorReduction\t0.75;")
snappyQual.append("\n\tminTetQuality\t1E-30;")
snappyQual.append("\n}")

import os as os
cwd = os.getcwd().replace("\\","/")

os.mkdir("0")
os.mkdir("constant")
os.mkdir("constant/polyMesh")
os.mkdir("constant/triSurface")
os.mkdir("system")

iF = cwd+"/Input.VSGR"
iFR = open(iF,"r")
currline = iFR.readline()

STLFileList = list()
RefinementList = list()

## Find out the names of the STL files and refinement regions and the levels of refinement

while currline!='':
	if "STL start" in currline:
		currline = iFR.readline()
		while "STL end" not in currline:
			dummy = currline.strip().split(" ",1)
			print dummy,len(dummy)
			name = dummy[0]
			if len(dummy) == 2:
				refSurf = dummy[1]
			else:	
				refSurf = "(2 3)"
			STLFileList.append(name + "\t" + refSurf)
			currline = iFR.readline()
	elif "Refinement start" in currline:
		currline = iFR.readline()
		while "Refinement end" not in currline:
			dummy = currline.strip().split()
			name = dummy[0]
			if len(dummy) == 6:
				mode = dummy[1]
				refVol = ' '.join(dummy[2:])
			elif len(dummy) == 3:
				mode = "inside"
				refVol = ' '.join(dummy[1:])
			else:
				mode = "inside"
				refVol = "((1.00 3))"
			RefinementList.append(name + "\t" + mode + "\t" + refVol)
			currline = iFR.readline()
	currline = iFR.readline()

iFR.close()

## Declare variables for the global maximum and minimum coordinates

gXMin = gYMin = gZMin = 1e10
gXMax = gYMax = gZMax = -1e10

## Find out the names of the solids in each STL file and the extents and write to MasterList

oF = cwd + "/MasterSTLList"
oFW = open(oF,"w")

for ii in xrange(0,len(STLFileList)):
	solidNames = list()
	localMinBounds = list()
	localMaxBounds = list()
	parts = STLFileList[ii].split("\t")
	name = parts[0]
	refSurf = parts[1]
	iF = cwd+"/"+name
	iFR = open(iF,"r")
	currline = iFR.readline()
	while currline!='':
		if ("solid" in currline) and ("endsolid" not in currline):
			xMin = yMin = zMin = 1e10
			xMax = yMax = zMax = -1e10
			solidNames.append(currline.strip().split()[1])
			currline = iFR.readline()
			while "endsolid" not in currline:		
				if "vertex" in currline:
					dummy = currline.strip().split()
					x = float(dummy[1])
					y = float(dummy[2])
					z = float(dummy[3])
					xMin = min(xMin,x)
					xMax = max(xMax,x)
					yMin = min(yMin,y)
					yMax = max(yMax,y)
					zMin = min(zMin,z)
					zMax = max(zMax,z)
				currline = iFR.readline()
			localMinBounds.append([xMin,yMin,zMin])
			localMaxBounds.append([xMax,yMax,zMax])
			gXMin = min(gXMin,xMin)
			gXMax = max(gXMax,xMax)
			gYMin = min(gYMin,yMin)
			gYMax = max(gYMax,yMax)
			gZMin = min(gZMin,zMin)
			gZMax = max(gZMax,zMax)
		currline = iFR.readline()
	iFR.close()
	for qq in xrange(0,len(solidNames)):
		oFW.write(name+"\t"+refSurf+"\t"+solidNames[qq]+"\t"+str(localMinBounds)+"\t"+str(localMaxBounds)+"\n")

oFW.close()
iFR.close()

## Set up the numbers for blockMeshDict

eps = 1e-3

gXMin = 0 if abs(gXMin) < eps else gXMin
gYMin = 0 if abs(gYMin) < eps else gYMin
gZMin = 0 if abs(gZMin) < eps else gZMin

gXMin = gXMin - 100
gYMin = gYMin - 100

gXMax = gXMax + 100
gYMax = gYMax + 100

gZMin = gZMin - 10 if gZMin!=0 else gZMin
gZMax = gZMax * 4

rX = gXMax - gXMin
rY = gYMax - gYMin
rZ = gZMax - gZMin

avgX = (gXMax + gXMin)/2
avgY = (gYMax + gYMin)/2
avgZ = (gZMax + gZMin)/2

dZ = 10
dX = int(rX/(rZ/dZ))
dY = int(rY/(rZ/dZ))

oF = cwd + "/constant/polyMesh/blockMeshDict"
oFW = open(oF,"w")

for ii in xrange(0,len(header)):
	oFW.write(header[ii])

for ii in xrange(0,len(blockMesh)):
	oFW.write(blockMesh[ii])

oFW.write("\n\nvertices\n(")
oFW.write("\n\t(" + str(gXMin) + "\t" + str(gYMax) + "\t" + str(gZMax) + ")")
oFW.write("\n\t(" + str(gXMax) + "\t" + str(gYMax) + "\t" + str(gZMax) + ")")
oFW.write("\n\t(" + str(gXMax) + "\t" + str(gYMin) + "\t" + str(gZMax) + ")")
oFW.write("\n\t(" + str(gXMin) + "\t" + str(gYMin) + "\t" + str(gZMax) + ")")
oFW.write("\n\t(" + str(gXMin) + "\t" + str(gYMax) + "\t" + str(gZMin) + ")")
oFW.write("\n\t(" + str(gXMax) + "\t" + str(gYMax) + "\t" + str(gZMin) + ")")
oFW.write("\n\t(" + str(gXMax) + "\t" + str(gYMin) + "\t" + str(gZMin) + ")")
oFW.write("\n\t(" + str(gXMin) + "\t" + str(gYMin) + "\t" + str(gZMin) + ")")
oFW.write("\n);")
oFW.write("\nblocks\n(\n\thex\t(0 1 2 3 4 5 6 7)\t");
oFW.write("("+str(dX)+" "+str(dY)+" "+str(dZ)+")"+"\tsimpleGrading\t(1 1 1)\n);")
oFW.write("\nedges\n(\n);")
oFW.write("\nboundary\n(")
oFW.write("\n\tLeft\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(0 4 7 3)\n\t\t);\n\t}")
oFW.write("\n\tRight\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(1 2 6 5)\n\t\t);\n\t}")
oFW.write("\n\tTop\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(0 3 2 1)\n\t\t);\n\t}")
oFW.write("\n\tBottom\n\t{\n\t\ttype\twall;\n\t\tfaces\n\t\t(\n\t\t\t(4 5 6 7)\n\t\t);\n\t}")
oFW.write("\n\tFront\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(0 1 5 4)\n\t\t);\n\t}")
oFW.write("\n\tBack\n\t{\n\t\ttype\tpatch;\n\t\tfaces\n\t\t(\n\t\t\t(3 7 6 2)\n\t\t);\n\t}")
oFW.write("\n);")
oFW.write("\nmergePatchPairs\n(\n);")
oFW.close()

## Setup and write snappyHexMeshDict

oF = cwd + "/system/snappyHexMeshDict" #"/constant/polyMesh/blockMeshDict"
oFW = open(oF,"w")

for ii in xrange(0,len(header)):
	oFW.write(header[ii])

for ii in xrange(0,len(snappyMesh)):
	oFW.write(snappyMesh[ii])

oFW.write("\n\ngeometry\n{\n")

iF = cwd + "/MasterSTLList"
iFR = open(iF,"r")
currline = iFR.readline()

namesSoFar = list()

ii = 0

wSTL = list()
wSolid = list()
wRefSurf = list()
wMinBd = list()
wMaxBd = list()

while currline!='':
	part = currline.strip().split("\t")
	name = (parts[0].split(".stl"))[0]
	refSurf = parts[1]
	solid = part[2]
	minBd = part[3]
	maxBd = part[4]
	wSTL.append(name)
	wSolid.append(solid)
	wRefSurf.append(refSurf)
	wMinBd.append(minBd)
	wMaxBd.append(maxBd)
	currline = iFR.readline()

c1 = 0

for ii in xrange(0,len(wSTL)):
	name = wSTL[ii]
	refSurf = wRefSurf[ii]
	solid = wSolid[ii]
	minBd = wMinBd[ii]
	maxBd = wMaxBd[ii]
	if name not in namesSoFar:
		if ii!=0:
			oFW.write("\n\t\t}\n\t}\n");
		namesSoFar.append(name)
		oFW.write("\n\t"+ name + ".stl\n\t{")
		oFW.write("\n\t\ttype\ttriSurfaceMesh;")
		oFW.write("\n\t\tname\t"+name+";")
		oFW.write("\n\t\tregions\n\t\t{")
		oFW.write("\n\t\t\t"+solid)
		oFW.write("\n\t\t\t{")
		oFW.write("\n\t\t\t\tname\t"+name+"_"+solid+";"+"\n\t\t\t}")
	else:
		oFW.write("\n\t\t\t{")
		oFW.write("\n\t\t\t\tname\t"+name+"_"+solid+";"+"\n\t\t\t}")
	if ii==len(wSTL)-1:
		oFW.write("\n\t\t}\n\t}\n");

print wSTL

for ii in xrange(0,len(RefinementList)):
	parts = RefinementList[ii]
	name = (parts.split("\t")[0]).split(".stl")[0]
	pos = wSTL.index(name)
	minBd = wMinBd[pos]
	maxBd = wMaxBd[pos]
	verticesMin = minBd.split()
	verticesMax = maxBd.split()
	minX = (verticesMin[0].split(",")[0]).split("[[")[1]
	minY = (verticesMin[1].split(","))[0]
	minZ = (verticesMin[2].split("]]"))[0]
	maxX = (verticesMax[0].split(",")[0]).split("[[")[1]
	maxY = (verticesMax[1].split(","))[0]
	maxZ = (verticesMax[2].split("]]"))[0]
	oFW.write("\n\t"+ name + "_Ref\n\t{")
	oFW.write("\n\t\ttype\tsearchableBox;")
	oFW.write("\n\t\tmin\t("+minX+" "+minY+" "+minZ+");")
	oFW.write("\n\t\tmax\t("+maxX+" "+maxY+" "+maxZ+");")
	oFW.write("\n\t}")

oFW.write("\n}")

oFW.write("\ncastellatedMeshControls")
oFW.write("\n{")
oFW.write("\n\tlocationInMesh\t("+str(avgX)+" "+str(avgY)+" "+str(avgZ)+");")
oFW.write("\n\tmaxLocalCells\t6000000;")
oFW.write("\n\tmaxGlobalCells\t20000000;")
oFW.write("\n\tminRefinementCells\t50;")
oFW.write("\n\tnCellsBetweenLevels\t3;")
oFW.write("\n\tresolveFeatureAngle\t60;")
oFW.write("\n\tallowFreeStandingZoneFaces\tfalse;")
oFW.write("\n\tfeatures")
oFW.write("\n\t(")

for ii in xrange(0,len(wSTL)):
	oFW.write("\n\t\t{")
	oFW.write('\n\t\t\tfile\t"'+wSTL[ii]+'.eMesh";')
	oFW.write('\n\t\t\tlevel\t2;')
	oFW.write("\n\t\t}")

oFW.write("\n\t);")

oFW.write("\n\trefinementSurfaces")
oFW.write("\n\t{")
namesSoFar = list()
for ii in xrange(0,len(wSTL)):
	name = wSTL[ii]
	refSurf = wRefSurf[ii]
	solid = wSolid[ii]
	if name not in namesSoFar:
		if ii!=0:
			oFW.write("\n\t\t\t}\n\t\t}\n")
		oFW.write("\n\t\t"+name)
		oFW.write("\n\t\t{")
		oFW.write("\n\t\t\tlevel\t"+refSurf+";")
		oFW.write("\n\t\t\tregions")
		oFW.write("\n\t\t\t{")
		oFW.write("\n\t\t\t\t"+name+"_"+solid)
		oFW.write("\n\t\t\t\t{")
		oFW.write("\n\t\t\t\t\tlevel\t"+refSurf+";")
		oFW.write("\n\t\t\t\t\tpatchInfo")
		oFW.write("\n\t\t\t\t\t{")
		oFW.write("\n\t\t\t\t\t\ttype\twall;")
		oFW.write("\n\t\t\t\t\t}")
		oFW.write("\n\t\t\t\t}")
	else:
		oFW.write("\n\t\t\t\t"+name+"_"+solid)
		oFW.write("\n\t\t\t\t{")
		oFW.write("\n\t\t\t\t\tlevel\t"+refSurf+";")
		oFW.write("\n\t\t\t\t\tpatchInfo")
		oFW.write("\n\t\t\t\t\t{")
		oFW.write("\n\t\t\t\t\t\ttype\twall;")
		oFW.write("\n\t\t\t\t\t}")
		oFW.write("\n\t\t\t\t}")
	if ii==len(wSTL)-1:
		oFW.write("\n\t\t\t}\n\t\t}\n")
oFW.write("\n\t}")

oFW.write("\n\trefinementRegions")
oFW.write("\n\t{")
namesSoFar = list()
for ii in xrange(0,len(RefinementList)):
	parts = RefinementList[ii]
	name = (parts.split("\t")[0]).split(".stl")[0]
	mode = parts.split("\t")[1]
	refVol = parts.split("\t")[2]
	oFW.write("\n\t\t"+ name + "_Ref\n\t\t{")
	oFW.write("\n\t\t\tmode\t"+mode+";")
	oFW.write("\n\t\t\tlevels\t"+refVol+";")
	oFW.write("\n\t\t}")

oFW.write("\n\t}")
oFW.write("\n}")

for ii in xrange(0,len(snappyAdd)):
	oFW.write(snappyAdd[ii])

for ii in xrange(0,len(wSTL)):
	oFW.write("\n\t\t"+wSTL[ii])
	oFW.write("\n\t\t{")
	oFW.write("\n\t\t\tnSurfaceLayers\t2;")
	oFW.write("\n\t\t}")

oFW.write("\n\t}")
oFW.write("\n}")

for ii in xrange(0,len(snappyQual)):
	oFW.write(snappyQual[ii])

oFW.write("\ndebug\t0;")
oFW.write("\nmergeTolerance\t1E-05;")
oFW.close()

## Create the files of the variables

variables = ["p","p_rgh","k","omega","epsilon","T","nut","U"]

for var in variables:
	oF = cwd + "/0/" + var
	oFW = open(oF,"w")

	for ii in xrange(0,len(header)):
		oFW.write(header[ii])

	oFW.write("\tclass\tvolScalarField;")
	oFW.write("\n\tobject\t"+var+";")
	oFW.write("\n}")
	
	if var == "p" or var == "p_rgh":
		oFW.write("\n\ndimensions\t[0 2 -2 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0;")
	elif var == "k":
		oFW.write("\n\ndimensions\t[0 2 -2 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0.1;")
	elif var == "omega":
		oFW.write("\n\ndimensions\t[0 0 -1 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0.1;")
	elif var == "epsilon":
		oFW.write("\n\ndimensions\t[0 2 -3 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0.1;")
	elif var == "T":
		oFW.write("\n\ndimensions\t[0 0 0 1 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t298.15;")
	elif var == "nut":	
		oFW.write("\n\ndimensions\t[0 2 -1 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t0;")
	elif var == "U":	
		oFW.write("\n\ndimensions\t[0 1 -1 0 0 0 0];")
		oFW.write("\n\ninternalField\tuniform\t(0 -2 0);")

	oFW.write("\n\nboundaryField\n{\n")

	s1 = ["Left","Right","Front","Top"]
	for s in s1:
		oFW.write("\n\t"+s)
		oFW.write("\n\t{")
		if var == "p" or var == "p_rgh" or var == "T":
			oFW.write("\n\t\ttype\tzeroGradient;")
		elif var == "k" or var == "omega" or var == "epsilon":
			oFW.write("\n\t\ttype\tfixedValue;")
			oFW.write("\n\t\tvalue\tuniform 0.1;")
		elif var == "nut":
			oFW.write("\n\t\ttype\tcalculated;")
			oFW.write("\n\t\tvalue\tuniform 0;")
		elif var == "U":
			oFW.write("\n\t\ttype\tfixedValue;")
			oFW.write("\n\t\tvalue\tuniform (0 -2 0);")	
		oFW.write("\n\t}")
	
	oFW.write("\n\tBack")
	oFW.write("\n\t{")
	if var == "p" or var == "p_rgh":
		oFW.write("\n\t\ttype\tfixedValue;")
		oFW.write("\n\t\tvalue\tuniform 0;")
	elif var == "k" or var == "omega" or var == "epsilon":
		oFW.write("\n\t\ttype\tzeroGradient;")
	elif var == "nut":
		oFW.write("\n\t\ttype\tcalculated;")
		oFW.write("\n\t\tvalue\tuniform 0;")
	elif var=="T":
		oFW.write("\n\t\ttype\tzeroGradient;")
	elif var == "U":
		oFW.write("\n\t\ttype\tzeroGradient;")
	oFW.write("\n\t}")

	for jj in xrange(0,len(wSTL)):
		name = wSTL[jj]
		solid = wSolid[jj]
		oFW.write("\n\t"+name+"_"+solid)
		oFW.write("\n\t{")
		if var == "p" or var == "p_rgh" or var == "T":
			oFW.write("\n\t\ttype\tzeroGradient;")
		elif var == "k":
			oFW.write("\n\t\ttype\tkqRWallFunction;")
			oFW.write("\n\t\tvalue\tuniform 0.1;") 
		elif var == "omega":
			oFW.write("\n\t\ttype\tomegaWallFunction;")
			oFW.write("\n\t\tvalue\tuniform 0.1;") 
		elif var == "epsilon":
			oFW.write("\n\t\ttype\tepsilonWallFunction;")
			oFW.write("\n\t\tvalue\tuniform 0.1;")
		elif var == "nut":
			oFW.write("\n\t\ttype\tnutkWallFunction;")
			oFW.write("\n\t\tvalue\tuniform 0;")
		elif var == "U":
			oFW.write("\n\t\ttype\tfixedValue;")
			oFW.write("\n\t\tvalue\tuniform (0 0 0);")	
		oFW.write("\n\t}")

	oFW.write("\n}")
	oFW.close()

## Create controlDict, decomposeParDict, fvSchemes and fvSolutions

files = ["controlDict","decomposeParDict","fvSchemes","fvSolution"]

for f in files:
	oF = cwd + "/system/" + f
	oFW = open(oF,"w")
	
	for ii in xrange(0,len(header)):
		oFW.write(header[ii])
	
	oFW.write("\tclass\tdictionary;")
	oFW.write("\n\tlocation\tsystem;")
	oFW.write("\n\tobject\t"+f+";")
	oFW.write("\n}\n")

	if f == "controlDict":
		oFW.write("\napplication\tsimpleFoam;")
		oFW.write("\nstartFrom\tstartTime;")
		oFW.write("\nstartTime\t0;")
		oFW.write("\nstopAt\tendTime;")
		oFW.write("\nendTime\t1000;")
		oFW.write("\ndeltaT\t1;")
		oFW.write("\nwriteControl\ttimeStep;")
		oFW.write("\nwriteInterval\t100;")
		oFW.write("\npurgeWrite\t0;")
		oFW.write("\nwriteFormat\tascii;")
		oFW.write("\nwritePrecision\t6;")
		oFW.write("\nwriteCompression\tcompressed;")
		oFW.write("\ntimeFormat\tgeneral;")
		oFW.write("\ntimePrecision\t6;")
		oFW.write("\nrunTimeModifiable\ttrue;")
	elif f == "decomposeParDict":
		oFW.write("\nnumberOfSubdomains\t1;")
		oFW.write("\n\nmethod\thierarchical;\n")
		oFW.write("\nsimpleCoeffs")
		oFW.write("\n{")
		oFW.write("\n\tn\t(1 1 1);")
		oFW.write("\n\tdelta\t0.001;")
		oFW.write("\n}\n")
		oFW.write("\nhierarchicalCoeffs")
		oFW.write("\n{")
		oFW.write("\n\tn\t(1 1 1);")
		oFW.write("\n\tdelta\t0.001;")
		oFW.write("\n\torder\txyz;")
		oFW.write("\n}\n")
		oFW.write("\nmanualCoeffs")
		oFW.write("\n{")
		oFW.write('\n\tdataFile\t"";')
		oFW.write("\n}")
	elif f == "fvSchemes":
		oFW.write("\n\nddtSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tsteadyState;")
		oFW.write("\n}")

		oFW.write("\n\ngradSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tGauss linear;")
		oFW.write("\n}")

		oFW.write("\n\ndivSchemes")
		oFW.write("\n{")	
		oFW.write("\n\tdefault\tnone;")
		oFW.write("\n\tdiv(phi,U)\tGauss upwind;")
		oFW.write("\n\tdiv(phi,T)\tGauss upwind;")
		oFW.write("\n\tdiv(phi,k)\tGauss upwind;")
		oFW.write("\n\tdiv(phi,epsilon)\tGauss upwind;")
		oFW.write("\n\tdiv(phi,omega)\tGauss upwind;")
		oFW.write("\n\tdiv((nuEff*dev(T(grad(U)))))\tGauss linear;")
		oFW.write("\n}")

		oFW.write("\n\nlaplacianSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tnone;")
		oFW.write("\n\tlaplacian(nuEff,U)\tGauss linear corrected;")
		oFW.write("\n\tlaplacian(kappaEff,T)\tGauss linear corrected;")
		oFW.write("\n\tlaplacian(DkEff,k)\tGauss linear corrected;")
		oFW.write("\n\tlaplacian(DepsilonEff,epsilon)\tGauss linear corrected;")
		oFW.write("\n\tlaplacian(DomegaEff,omega)\tGauss linear corrected;")
		oFW.write("\n\tlaplacian((1|A(U)),p)\tGauss linear corrected;")
		oFW.write("\n\tlaplacian((1|A(U)),p_rgh)\tGauss linear corrected;")
		oFW.write("\n}")

		oFW.write("\n\ninterpolationSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tlinear;")
		oFW.write("\n}")

		oFW.write("\n\nsnGradSchemes")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tcorrected;")
		oFW.write("\n}")

		oFW.write("\n\nfluxRequired")
		oFW.write("\n{")
		oFW.write("\n\tdefault\tno;")
		oFW.write("\n\tp	;")
		oFW.write("\n\tp_rgh	;")
		oFW.write("\n}")
	
	elif f == "fvSolution":
		oFW.write("\n\nsolvers")
		oFW.write("\n{")
		oFW.write('\n\t"p|p_rgh"')
		oFW.write("\n\t{")
		oFW.write("\n\t\tsolver\tGAMG;")
		oFW.write("\n\t\tsmoother\tGaussSeidel;")
		oFW.write("\n\t\ttolerance\t1e-08;")
		oFW.write("\n\t\trelTol\t0.05;")
		oFW.write("\n\t\tcacheAgglomeration\toff;")
		oFW.write("\n\t\tnCellsInCoarsestLevel\t20;")
		oFW.write("\n\t\tagglomerator\tfaceAreaPair;")
		oFW.write("\n\t\tmergeLevels\t1;")
		oFW.write("\n\t}")

		oFW.write('\n"k|omega|epsilon|U|T"')
		oFW.write("\n\t{")
		oFW.write("\n\t\tsolver\tPBiCG;")
		oFW.write("\n\t\tpreconditioner\tDILU;")
		oFW.write("\n\t\ttolerance\t1e-05;")
		oFW.write("\n\t\trelTol\t0.1;")
		oFW.write("\n\t}")
		oFW.write("\n}")

		oFW.write("\nSIMPLE")
		oFW.write("\n{")
		oFW.write("\n\tnNonOrthogonalCorrectors\t0;")
		oFW.write("\n\tpRefCell\t0;")
		oFW.write("\n\tpRefValue\t0;")
		oFW.write("\n\tresidualControl")
		oFW.write("\n\t{")
		oFW.write('\n\t\t"k|omega|epsilon|U|T"\t1e-4;')
		oFW.write('\n\t\t"p|p_rgh"\t1e-3;')
		oFW.write("\n\t}")
		oFW.write("\n}")

		oFW.write("\nrelaxationFactors")
		oFW.write("\n{")
		oFW.write("\n\tequations")
		oFW.write("\n\t{")
		oFW.write('\n\t\t"k|omega|epsilon|U|T"\t0.3;')
		oFW.write("\n\t}")
	
		oFW.write("\n\tfields")
		oFW.write("\n\t{")
		oFW.write('\n\t\t"p|p_rgh"\t0.7;')
		oFW.write("\n\t}")
		
		oFW.write("\n}")
	oFW.close()

## Create RASProperties and transportProperties files

oF = cwd + "/constant/RASProperties"
oFW = open(oF,"w")

for ii in xrange(0,len(header)):
	oFW.write(header[ii])

oFW.write("\tclass\tdictionary;")
oFW.write("\n\tlocation\tconstant;")
oFW.write("\n\tobject\tRASProperties;")
oFW.write("\n}\n")

oFW.write("\n\nRASModel\tkEpsilon;")
oFW.write("\n\nturbulence\ton;")
oFW.write("\n\nprintCoeffs\ton;")

oFW.close()

oF = cwd + "/constant/transportProperties"
oFW = open(oF,"w")

for ii in xrange(0,len(header)):
	oFW.write(header[ii])

oFW.write("\tclass\tdictionary;")
oFW.write("\n\tlocation\tconstant;")
oFW.write("\n\tobject\ttransportProperties;")
oFW.write("\n}\n")

oFW.write("\n\ntransportModel\tNewtonian;")
oFW.write("\n\nnu\tnu\t[0 2 -1 0 0 0 0]\t1.5E-05;")

oFW.close()
