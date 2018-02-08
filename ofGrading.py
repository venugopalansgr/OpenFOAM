# Purpose: Simple tool to compute either the number of cells or the first cell length
# Purpose: For use with graded blockMesh in OpenFOAM
# Usage: python ofGrading.py -s <length> -a <first cell size> -R <size ratio between last cell and first cell>
# Usage: python ofGrading.py -s <length> -R <size ratio between last cell and first cell> -n <number of cells>
# Output: Prints on screen the "missing value" - either number of cells or first cell length
# Notes: For multiple fields for the same plane, please copy paste the normals and point on the plane into a new line with new field/variable
# Notes: Can use long form of options as well: --dist= (instead of -s), --numcell= (instead of -n), --ratio= (instead of -R) and --first= (instead of -a)
# Notes: No warranty on results. Use at your own risk/discretion
# Notes: Code is free. Appreciate feedback/acknowledging when using it

# Created by: Venugopalan Raghavan

# Version 0: 08 February 2018
#	Release Notes: Initial Release


# Here is the schematic of the geometry used:
# |-|--|----|--------|-----------|----------------|
#  a ar ar^2   ar^3                   ar^(n-1)
# ratio = R = ar^(n-1) / a = r^(n-1)
# s = a + ar + ar^2 + ar^3 + ... + ar^(n-1) = (a(r^n - 1))/(r - 1)


import math
import getopt
import sys

eps = 1e-4

def getA(S,R,n):
	if (math.fabs(R - eps) < eps):
		a = S
	elif (math.fabs(R - 1) < eps):
		a = S
	elif n!=1:
		rec = 1./(n-1)
		recn = n*rec
		a = (S*(R**rec - 1))/(R**recn - 1)
	return a

def getN(S,A,R):
	if (math.fabs(R - eps) < eps):
		n = 1
	elif (math.fabs(R - 1) < eps):
		n = math.ceil(S/A)
	elif (math.fabs(S - A) < eps):
		n = 1
	else:
		t1 = (S - A)/(S - A*R)
		n = int(math.ceil(1 + math.log(R)/math.log(t1)))
		
	return n		

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:],"hs:R:a:n:",["help","dist=","ratio=","first=","numcell="])
	except:
		print "ofGrading can be used in one of two modes:"
		print "1. Get number of cells"
		print "python ofGrading.py -s <length> -a <first cell size> -R <size ratio between last cell and first cell>"
		print "2. Get initial cell size"
		print "python ofGrading.py -s <length> -R <size ratio between last cell and first cell> -n <number of cells>"
		sys.exit(2)
	
	sFlag = False
	aFlag = False
	nFlag = False
	rFlag = False
	
	getAFlag = False
	getNFlag = False
	
	for opt, arg in opts:
		if opt == "-h":
			print "ofGrading can be used in one of two modes:"
			print "1. Get number of cells"
			print "python ofGrading.py -s <length> -a <first cell size> -r <size ratio between last cell and first cell>"
			print "2. Get initial cell size"
			print "python ofGrading.py -s <length> -r <size ratio between last cell and first cell> -n <number of cells>"			
			sys.exit()
		elif opt in ("-s", "--dist"):
			S = float(arg)
			if S <= 0:
				print "Can't have negative or zero length mate!"
				sys.exit()
			sFlag = True
		elif opt in ("-a", "--first"):
			A = float(arg)
			if A <= 0:
				print "Can't have negative or zero first cell size boss!"
				sys.exit()
			aFlag = True
		elif opt in ("-R", "--ratio"):
			R = float(arg)
			if R < 0:
				print "Erm. Negative ratios don't make sense here"
				sys.exit()
			rFlag = True
		elif opt in ("-n", "--numcell"):
			n = int(arg)
			print "This will be rounded to the nearest integer"
			if n < 0:
				print "How is it even possible to have a negative number of cells?!"
				sys.exit()
			nFlag = True
			
	if sFlag == False:
		print "Wah! Give me a length to work with, buddy!"
		print "Type -h to get the available options"
		sys.exit()

	if rFlag == False:
		print "Come on! I do need the ratio you know!"
		print "Type -h to get the available options"
		sys.exit()
		
	if ((aFlag == True) and (nFlag == True)):
		print "You think this funny? :\\"
		print "Give me first division size OR total number of cells. NOT BOTH!"
		print "I will assume that you need number of cells"
		getNFlag = True
	elif nFlag==False:
		print "Requested for number of cells"
		getNFlag = True
	elif aFlag==False:
		print "Requested for initial division size"
		getAFlag = True
		
	if getAFlag==True:
		A = getA(S,R,n)
		print "\nFirst cell length is: ", A
	elif getNFlag==True:
		n = getN(S,A,R)
		print "\nNumber of cells is: ", n
		
if __name__ == "__main__":
	main()
