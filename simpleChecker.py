# Usage: python simpleChecker.py <case name/location>
# Input: No input
# Purpose: Performs basic check on OpenFOAM case to ensure that bare minimal files are present; Also flags basic errors if present
# Notes: The <> argument is optional. If not given, assumes the code is called from the case folder itself 
# Notes: Specifying only case name implies that the current working directory has more than one folder 
# Notes: At the moment limited to mainly incompressible solvers such as simpleFoam, pimpleFoam, buoyantSimpleFoam, buoyantBoussinesqSimpleFoam
# Notes: Future versions may have greater variety
# Notes: No warranty on results. Use at your own risk/discretion
# Notes: Code is free. Appreciate feedback/acknowledging when using it
# Created by: Venugopalan Raghavan

import os 
import sys

# Get the input: python simpleChecker.py <case name/location>
argv_cnt = len(sys.argv)

if argv_cnt == 1:
	s = sys.argv[0]
	cwd = os.getcwd().replace("\\","/") 
else:
	s = sys.argv[1]
	if "\\" or "/" in s:
		cwd = str(s)
	else:
		cwd = os.getcwd().replace("\\","/") + str(s)

dirList = list()

######################### basic folders	#########################
# Check if 0, constant and system folders are present
print "Checking basic folders"
dirList = os.listdir(cwd)

if "0" not in dirList:
	print "The 0 folder is missing!"
	quit()
elif "constant" not in dirList:
	print "The constant folder is missing!"
	quit()
elif "system" not in dirList:
	print "The system folder is missing!"
	quit()

print "0 folder present....\tOK"
print "constant folder present....\tOK"
print "system folder present....\tOK"	

# Check if triSurface and polyMesh folders present in /constant folder
dirList = os.listdir(cwd+"/constant")	
if "triSurface" not in dirList:
	print "The constant/triSurface folder is missing!"
	quit()
elif "polyMesh" not in dirList:
	print "The constant/polyMesh folder is missing!"
	quit()

print "constant/polyMesh folder present....\tOK"
print "constant/triSurface folder present....\tOK"
print "All basic directories for OpenFOAM are present!"

######################### /constant folder	#########################
print "\nChecking files in constant folder"
# Check if RASProperties or LESProperties in /constant folder
if ("RASProperties" not in dirList) and ("LESProperties" not in dirList):
	print "RASProperties/LESProperties file missing!"
	quit()

# Check if blockMeshDict present in constant/polyMesh directory	
dirList = os.listdir(cwd + "/constant/polyMesh")
constList = ["blockMeshDict"]
c = filter(lambda x: (not(x in constList)),dirList)

if len(c)==0:
	print "No files missing in constant/polyMesh directory!"
else:
	print "File: " + str(c) + " missing in constant/polyMesh directory!"
	quit()

# Read the contents of blockMeshDict
iF = cwd + "/constant/polyMesh/blockMeshDict"
iFR = open(iF,"r")
currline = iFR.readline()

filesCheck = list()

while "boundary" not in currline:
	currline = iFR.readline()
	if "boundary" in currline:
		break

if "(" in currline:
	currline = iFR.readline()
	cnt = 1
else:
	currline = iFR.readline()
	currline = iFR.readline()
	cnt = 1

bNameList = list()
bTypeList = list()	

# Get boundary names and types
while currline!='':
	if (cnt == 1) and (");" not in currline):
		bNameList.append(currline.strip())
		cnt = 0
	if "}" in currline:
		cnt = 1
	elif "type" in currline:
		bTypeList.append(((currline.strip().split())[1].split(";"))[0])
	currline = iFR.readline()
	if "mergePatchPairs" in currline:
		break
		
iFR.close()

# Get turbulence model and if turbulence has been switched on
if "RASProperties" in os.listdir(cwd+"/constant"):
	iF = cwd + "/constant/RASProperties"
	iFR = open(iF,"r")
	currline = iFR.readline()
	while currline!='':
		currline = iFR.readline()
		if "RASModel" in currline:
			break
	turbModel = (currline.strip().split()[1]).split(";")[0]
	currline = iFR.readline()
	while currline!='':
		if "turbulence" in currline:
			break
		currline = iFR.readline()
	turbSwitch = (currline.strip().split()[1]).split(";")[0]
	if turbSwitch == "off":
		print "Turbulence has been switched off!"
	iFR.close()
else:
	turbModel = ""

# Get variables associated with turbulence model
if turbModel in ["kEpsilon","RNGkEpsilon","NonlinearKEShih","LienCubicKE","LaunderSharmaKE","LamBremhorstKE","LienCubicKELowRe","LienLeschzinerLowRe","realizableKE"]:
	turbulence = "kEpsilon"
	filesCheck.append("k")
	filesCheck.append("epsilon")
elif turbModel in ["kOmega","kOmegaSST"]:
	turbulence = "kOmega"
	filesCheck.append("k")
	filesCheck.append("omega")
elif turbModel in "SpalartAllmaras":
	turbulence = "SA"
	filesCheck.append("nut")
	
xList = list()
		
if turbulence == "kEpsilon":
	xList = ["k","epsilon"]
elif turbulence == "kOmega":
	xList = ["k","omega"]
elif turbulence == "SA":
	xList = ["nut"]

######################### /system folder	#########################
print "\nChecking files in system folder"	
# Get solver name
iF = cwd + "/system/controlDict"
iFR = open(iF,"r")
currline = iFR.readline()

while currline!='':
	if "application" in currline:
		break
	currline = iFR.readline();

# Add variables particular to solver to later check in 0 directory
solver = ((currline.strip().split()[1]).split(";"))[0]
if solver in ["simpleFoam","pimpleFoam","MRFSimpleFoam"]:
	filesCheck.append("U")
	filesCheck.append("p")
elif solver in ["buoyantSimpleFoam","buoyantPimpleFoam","buoyantBoussinesqSimpleFoam","buoyantBoussinesqPimpleFoam"]:
	filesCheck.append("U")
	filesCheck.append("p_rgh")
	filesCheck.append("T")

# Get the name & type of the STL geometries in snappyHexMeshDict
gNameList = list()
gTypeList = list()
STLNameList = list()

iF = cwd + "/system/snappyHexMeshDict"
iFR = open(iF,"r")
currline = iFR.readline()
cnt = 0

while currline!='':
	if "name" in currline:
		if cnt == 1:
			name = (currline.strip().split()[1]).split(";")[0]
			gNameList.append(name)
			cnt = 0
		else:
			name = (currline.strip().split()[1]).split(";")[0]
			STLNameList.append(name+".stl")
			cnt = 1
	elif "castellatedMeshControls" in currline:
		break
	currline = iFR.readline()

#print "Geometry Names = ",gNameList
	
while currline!='':
	for x in gNameList:
		if x in currline:
			while "type" not in currline:
				currline = iFR.readline()
			typ = (currline.strip().split()[1]).split(";")[0]
			gTypeList.append(typ)
	currline = iFR.readline()

iFR.close()

#print "Geometry Types = ",gTypeList	

# Check if all basic files present in /system directory
dirList = os.listdir(cwd + "/system")
sysList = ["decomposeParDict","fvSolution","fvSchemes","snappyHexMeshDict","controlDict"]
c = filter(lambda x: (not(x in sysList)),dirList)

if len(c)==0:
	print "No files missing in system directory!"
else:
	print "File: " + str(c) + " missing in system directory!"
	quit()

# Check that the number of subdomains in decomposePar is fine
iF = cwd + "/system/decomposeParDict"
iFR = open(iF,"r")
currline = iFR.readline()

while currline!='':
	if "numberOfSubdomains" in currline:
		break
	currline = iFR.readline()
	
nDF = float((currline.strip().split()[1]).split(";")[0])
nDI = int(nDF)#((currline.strip().split()[1]).split(";")[0])

if (abs(nDI-nDF)>1e-5):
	print "Number of sub-domains is not an integer!"

if nDI < 0:
	print "Number of sub-domains cannot be negative!"
	
while currline!='':
	if "method" in currline:
		break
	currline = iFR.readline()

methName = (currline.strip().split()[1]).split(";")[0]

yes1 = 0

while currline!='':
	if methName!="scotch":
		if methName+"Coeffs" in currline:
			yes1 = 1
		if ("n" in currline) and (yes1==1):
			break
	currline = iFR.readline()

nz = int((currline.strip().split()[3]).split(")")[0])
ny = int(currline.strip().split()[2])
nx = int((currline.strip().split()[1]).split("(")[1])

if (nx*ny*nz != nDI):
	print "Error in specified splitting of domains!"
	print "\tnumberOfSubdomains =",nDI,"but computed =",(nx*ny*nz)

######################### /0 folder	#########################
print "\nChecking files in constant/triSurface folder"
# Check that the stl files mentioned tally with those in constant/triSurface

dirList = os.listdir(cwd+"/constant/triSurface")
c = filter(lambda x: (not(x in dirList)),STLNameList)

if len(c)==0:
	print "No files missing in constant/triSurface directory!"
else:
	print "File: " + str(c) + " missing in constant/triSurface directory!"
	
######################### /0 folder	#########################
print "\nChecking files in 0 folder"

# Check if relevant files are present in 0 directory
dirList = os.listdir(cwd+"/0")
c = filter(lambda x: (not(x in dirList)),filesCheck)

if len(c)==0:
	print "No files missing in 0 directory!"
else:
	print "File: " + str(c) + " missing in 0 directory!"
	quit()

# Check the contents of each of the files in 0 to ensure BC present
for x in filesCheck:
	iF = cwd + "/0/" + x
	iFR = open(iF,"r")
	combinedObjList = list(bNameList+gNameList)
	currline = iFR.readline()
	while currline!='':
		for y in combinedObjList:
			if y in currline:
				combinedObjList.remove(y)
		currline = iFR.readline()
	if len(combinedObjList) > 0:
		print "Missing objects: ",combinedObjList,"in file: ",x
	iFR.close()
	
for x in xList:
	combinedObjList = list(bNameList+gNameList)
	combinedTypList = list(bTypeList+gTypeList)
	iF = cwd + "/0/" + x
	iFR = open(iF,"r")
	currline = iFR.readline()
	while currline!='':
		for y in combinedObjList:
			ind = combinedObjList.index(y)
			if y in currline:
				while "type" not in currline:
					currline = iFR.readline()
				typ = (currline.strip().split()[1]).split(";")[0]
				yes1 = "WallFunction" in typ
				yes2 = (combinedTypList[ind] == "patch")
				yes3 = (solver in ["buoyantSimpleFoam","buoyantPimpleFoam"])
				yes4 = "compressible::" in typ
				if (yes1 and yes2):
					print "Cannot assign wall function to type patch!","Patch name = ",y
				if (yes3 and not(yes4)):
					print "Compressible solvers need 'compressible::' before WallFunction!"
				yes5 = (x if (x in ["epsilon","omega"]) else (x+"qR" if x in "k" else (x+"k"))) in typ
				if (not(yes5) and yes1):
					print "Incorrect wall function for",x,"!","Obj = ",y
		currline = iFR.readline()
	iFR.close()

for x in filesCheck:
	iF = cwd + "/0/" + x
	iFR = open(iF,"r")
	currline = iFR.readline()
	
	while currline!='':
		if "dimensions" in currline:
			break
		currline = iFR.readline()
	
	s = currline.strip().split()
	dims = s[1]+" "+s[2]+" "+s[3]+" "+s[4]+" "+s[5]+" "+s[6]+" "+s[7].split(";")[0]
	iFR.close()
	
	if ((x=="k") and (dims != "[0 2 -2 0 0 0 0]")):
		print "Error in units of k!"
	elif ((x=="epsilon") and (dims != "[0 2 -3 0 0 0 0]")):
		print "Error in units of epsilon!"
	elif ((x=="omega") and (dims != "[0 0 -1 0 0 0 0]")):
		print "Error in units of omega!"
	elif ((x=="U") and (dims != "[0 1 -1 0 0 0 0]")):
		print "Error in units of U!"
	elif ((x=="T") and (dims != "[0 0 0 1 0 0 0]")):
		print "Error in units of T!"
	elif ((x=="nut") and (dims != "[0 2 -1 0 0 0 0]")):
		print "Error in units of nut!"
	elif ((x=="p") or (x=="p_rgh")):
		if solver in ["buoyantSimpleFoam","buoyantPimpleFoam"]:
			if (dims!="[1 -1 -2 0 0 0 0]"):
				print "Error in units of p!"
				print "Compressible case should be: Pa"
		else:
			if (dims!="[0 2 -2 0 0 0 0]"):
				print "Error in units of p!"
				print "Incompressible case should be: m2/s2"
