## Purpose: Python code to create topoSetDict with rotatedBoxToCell option
## Usage: Run the code from command line wherever the STL file is located
## Usage: python rotBoxToCellCreate.py <STLFileName> 
## Notes: Needs STL file (one file with multiple solids is fine) of the bounding box
## Notes: If using SketchUp, refer to "skp_to_grpstl.rb" to export multiple solids in single STL
## Notes: Requires "numpy" library to be installed! Tested with Python 2.7
## Notes: No warranty on results. Use at your own risk/discretion
## Notes: Code is free. Appreciate feedback/acknowledging when using it
## Created by: Venugopalan Raghavan

import os as os
import numpy as np
import math as math
import sys as sys

cwd = str(os.getcwd()).replace("\\","/")

ins = sys.argv

fileName = ins[len(ins)-1]

iF = cwd + "/" + fileName
iFR = open(iF,"r")
currline = iFR.readline()

namedir = list()

while currline!='':
	if "solid" in currline:
		name = currline.strip().split()[1]
		if name not in namedir:
			namedir.append(name)
	currline = iFR.readline()

iFR.close()

oF = cwd + "/topoSetDict"
oFW = open(oF,"w")

oFW.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
oFW.write("| =========                 |                                                 |\n")
oFW.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
oFW.write("|  \\    /   O peration     | Version:  2.1.x                                 |\n")
oFW.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
oFW.write("|    \\/     M anipulation  |                                                 |\n")
oFW.write("\\*---------------------------------------------------------------------------*/\n")
oFW.write("FoamFile\n")
oFW.write("{\n")
oFW.write("    version     2.0;\n")
oFW.write("    format      ascii;\n")
oFW.write("    class       dictionary;\n")
oFW.write("    object      topoSetDict;\n")
oFW.write("}\n")
oFW.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")

iFR = open(iF,"r")
currline = iFR.readline()

ii = 1
solidcnt = 0

while currline!='':
	yes1 = "solid" in currline
	yes2 = "endsolid" in currline
	xlist = list()
	ylist = list()
	zlist = list()
	pts = list()
	
	xMax = -1e10
	yMax = -1e10
	zMax = -1e10
	
	xMin = 1e10
	yMin = 1e10
	zMin = 1e10
	if yes1 and ~yes2:
		name = currline.strip().split()[1]
		currline = iFR.readline()
		yes2 = "endsolid" in currline
		while ~yes2:
			if "vertex" in currline:
				s1 = currline[currline.index("x")+2:len(currline)]
				s2 = s1.strip().split()
				x = float(s2[0])
				y = float(s2[1])
				z = float(s2[2])
				if x > xMax:
					i1 = [x,y]
					xMax = x
				if x < xMin:
					i2 = [x,y]
					xMin = x
				if y > yMax:
					j1 = [x,y]
					yMax = y
				if y < yMin:
					j2 = [x,y]
					yMin = y
				if z > zMax:
					zMax = z
				if z < zMin:
					zMin = z
			currline = iFR.readline()
			yes2 = "endsolid" in currline
			if yes2:
				break
			#print currline,yes2
		
		if ii==1:
			oFW.write("actions\n(\n")
			ii = ii + 1
	
		oFW.write("\t{\n")
		oFW.write("\t\tname "+name+"C;\n")
		oFW.write("\t\ttype cellSet;\n")
		oFW.write("\t\taction new;\n")
		oFW.write("\t\tsource rotatedBoxToCell;\n")
		oFW.write("\t\tsourceInfo\n")
		oFW.write("\t\t{\n")
	
		o = j2
		vec1 = [i1[0]-o[0],i1[1]-o[1]]
		vec2 = [j1[0]-o[0],j1[1]-o[1]]
		vec3 = [i2[0]-o[0],i2[1]-o[1]]
	
		magVec1 = math.sqrt(vec1[0]**2 + vec1[1]**2)
		magVec2 = math.sqrt(vec2[0]**2 + vec2[1]**2)
		magVec3 = math.sqrt(vec3[0]**2 + vec3[1]**2)
	
		dummy = list()
		dummy.append(magVec1)
		dummy.append(magVec2)
		dummy.append(magVec3)
		tt = dummy.index(max(dummy))
	
		if tt==0:
			ang2 = math.degrees(math.acos(np.dot(vec2,[1,0])/magVec2))
			ang3 = math.degrees(math.acos(np.dot(vec3,[1,0])/magVec3))
			if ang2 < ang3:
				iVec = vec2
				jVec = vec3
			else:
				iVec = vec3
				jVec = vec2
		elif tt==1:
			ang1 = math.degrees(math.acos(np.dot(vec1,[1,0])/magVec1))
			ang3 = math.degrees(math.acos(np.dot(vec3,[1,0])/magVec3))
			if ang1 < ang3:
				iVec = vec1
				jVec = vec3
			else:
				iVec = vec3
				jVec = vec1
		else:
			ang1 = math.degrees(math.acos(np.dot(vec1,[1,0])/magVec1))
			ang2 = math.degrees(math.acos(np.dot(vec2,[1,0])/magVec2))
			if ang1 < ang2:
				iVec = vec1
				jVec = vec2
			else:
				iVec = vec2
				jVec = vec1
	
		origin = "(" + str(o[0]) + " " + str(o[1]) + " " + str(zMin) + ")"
		iPts = "(" + str(iVec[0]) + " " + str(iVec[1]) + " " + "0)"
		jPts = "(" + str(jVec[0]) + " " + str(jVec[1]) + " " + "0)"
		kPts = "(0 0 " + str(zMax-zMin) + ")"
		
		oFW.write("\t\t\torigin " + origin + ";\n")
		oFW.write("\t\t\ti " + iPts + ";\n")
		oFW.write("\t\t\tj " + jPts + ";\n")
		oFW.write("\t\t\tk " + kPts + ";\n")
		oFW.write("\t\t}\n")
		oFW.write("\t}\n")
	
		oFW.write("\t{\n")
		oFW.write("\t\tname "+name+";\n")
		oFW.write("\t\ttype cellZoneSet;\n")
		oFW.write("\t\taction new;\n")
		oFW.write("\t\tsource setToCellZone;\n")
		oFW.write("\t\tsourceInfo\n")
		oFW.write("\t\t{\n")
		oFW.write("\t\t\tset " + name + "C;\n")
		oFW.write("\t\t}\n")
		oFW.write("\t}\n")
	
		if "NS" in name:
			angle = "(" + str(math.cos(math.radians(22))) + " " + str(-math.sin(math.radians(22))) + " 0)"
		elif "EW" in name:
			angle = "(" + str(math.cos(math.radians(68))) + " " + str(math.sin(math.radians(68))) + " 0)"
	
			
		ii = ii+1

		print "Done ",name
		solidcnt = solidcnt + 1
		
	currline = iFR.readline()

print "Number of solids = ",solidcnt

oFW.write(");")

iFR.close()
oFW.close()
