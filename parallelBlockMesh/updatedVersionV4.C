/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockMesh

Description
    Parallel version of blockMesh.
    Current version supports single-block 
    Multi-grading in all directions works
    Mesh generation can be parallelized in one direction only
    This version does not support curved edges

    Uses the block mesh description found in
    \a system/blockMeshDict

Usage

    - mpirun -np <nprocs> parallelblockMesh [OPTION] -parallel

Notes:
	Following modifications required in blockMeshDict:
		- In patches, the patches with normal same as direction of cutting mesh must be specified at the end
		- Entries after simpleGrading entry must be given in manner similar to multi-grading i.e. (length num_of_cells expansion_ratio)
			- Length entries must be specified in absolute lengths

Created by:
	Venugopalan S.G. Raghavan (raghavanvsg@ihpc.a-star.edu.sg; venugopalan.raghavan@gmail.com)
	Dominic Chander
	Harish Gopalan

Date:
	v0: 24 Jun 2018
	v1: 27 Aug 2018
	v2: 11 Sep 2018
	v3: 11 Oct 2018
	v4: 17 Oct 2018

Comments added: 11 Oct 2019

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOdictionary.H"
#include "IOPtrList.H"

#include "blockMesh.H"
#include "attachPolyTopoChanger.H"
#include "emptyPolyPatch.H"
#include "cellSet.H"

#include "argList.H"
#include "OSspecific.H"
#include "OFstream.H"

#include "Pair.H"
#include "slidingInterface.H"
#include "mpi.h"
#include "mergePolyMesh.H"

using namespace Foam;


class multiMesh
{
    public:

    word regionName ;
    IOdictionary meshDict;
    Time & runTime;
    int & myRank;
    polyMesh* mesh;

    multiMesh
    (
	word & regionName,
	IOdictionary & meshDict,
	Time & runTime,
	int & myrank
    );

    ~multiMesh()
    {}

    void createMesh();
    void writeMesh();
};

multiMesh::multiMesh
(
	word & regionName_,
	IOdictionary & meshDict_,
	Time & runTime_,
	int & myRank_
)
:
regionName(regionName_),
meshDict(meshDict_),
runTime(runTime_),
myRank(myRank_)
{

}

// // // // // 
void* threadAccess(void* argAccess)
{
	reinterpret_cast<multiMesh*>(argAccess)->createMesh();
	return NULL;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void multiMesh::createMesh()
{
    blockMesh blocks(meshDict, regionName);

    word defaultFacesName = "defaultFaces";
    word defaultFacesType = emptyPolyPatch::typeName;
    word procName = "processor" + name(myRank);

    IOobject myob
    (
	".",
	runTime.constant(),
        runTime
    );

    this->mesh = new polyMesh
    (
	myob,
        xferCopy(blocks.points()),           // could we re-use space?
        blocks.cells(),
        blocks.patches(),
        blocks.patchNames(),
        blocks.patchDicts(),
        defaultFacesName,
        defaultFacesType,
	false
    );


    // Read in a list of dictionaries for the merge patch pairs
    if (meshDict.found("mergePatchPairs"))
    {
        List<Pair<word> > mergePatchPairs
        (
            meshDict.lookup("mergePatchPairs")
        );

        #include "mergePatchPairs.H"
    }
    else
    {
        if (myRank==0) Info<< nl << "There are no merge patch pairs edges" << endl;
    }


    // Set any cellZones (note: cell labelling unaffected by above
    // mergePatchPairs)

    label nZones = blocks.numZonedBlocks();

    if (nZones > 0)
    {
        Info<< nl << "Adding cell zones" << endl;

        // Map from zoneName to cellZone index
        HashTable<label> zoneMap(nZones);

        // Cells per zone.
        List<DynamicList<label> > zoneCells(nZones);

        // Running cell counter
        label cellI = 0;

        // Largest zone so far
        label freeZoneI = 0;

        forAll(blocks, blockI)
        {
            const block& b = blocks[blockI];
	
	    #ifdef _v5
	            const List<FixedList<label, 8>> blockCells = b.cells();
	    #endif
	    #ifdef _v4	
		    const labelListList blockCells = b.cells();
	    #endif

            const word& zoneName = b.zoneName();

            if (zoneName.size())
            {
                HashTable<label>::const_iterator iter = zoneMap.find(zoneName);

                label zoneI;

                if (iter == zoneMap.end())
                {
                    zoneI = freeZoneI++;

                    Info<< "    " << zoneI << '\t' << zoneName << endl;

                    zoneMap.insert(zoneName, zoneI);
                }
                else
                {
                    zoneI = iter();
                }

                forAll(blockCells, i)
                {
                    zoneCells[zoneI].append(cellI++);
                }
            }
            else
            {
                cellI += b.cells().size();
            }
        }


        List<cellZone*> cz(zoneMap.size());

        if (!myRank)	Info<< nl << "Writing cell zones as cellSets" << endl;

        forAllConstIter(HashTable<label>, zoneMap, iter)
        {
            label zoneI = iter();

            cz[zoneI] = new cellZone
            (
                iter.key(),
                zoneCells[zoneI].shrink(),
                zoneI,
                mesh->cellZones()
            );

            // Write as cellSet for ease of processing
            cellSet cset(*mesh, iter.key(), zoneCells[zoneI].shrink());
            cset.write();
        }

        mesh->pointZones().setSize(0);
        mesh->faceZones().setSize(0);
        mesh->cellZones().setSize(0);
        mesh->addZones(List<pointZone*>(0), List<faceZone*>(0), cz);
    }
}


void multiMesh::writeMesh()
{
    // Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< nl << "Writing polyMesh" << endl;

    mesh->removeFiles();
    if (!mesh->write())
    {
        FatalErrorIn("blockMesh")
            << "Failed writing polyMesh."
            << exit(FatalError);
    }
}


int main(int argc, char *argv[])
{
    string dateString = clock::date();
    string timeString = clock::clockTime();

    argList::addBoolOption
    (
        "blockTopology",
        "write block edges and centres as .obj files"
    );
    argList::addOption
    (
        "dict",
        "file",
        "specify alternative dictionary for the blockMesh description"
    );
    argList::addOption
    (	
	"vector",
	"vector",
	"specify direction to split the meshes"
    );

    Foam::mkDir(Foam::cwd() + "/processor0");
    int myRank, nThreads;

    #include "setRootCase.H"

   if (!Pstream::master() && !isDir(args.path()))
    {
        mkDir(args.path());
    }

    MPI_Comm_rank(MPI_COMM_WORLD,&myRank);
    MPI_Comm_size(MPI_COMM_WORLD,&nThreads);

    double t1 = MPI_Wtime();

    #include "createTimeMPI.H"

    const word dictName("blockMeshDict");

    word regionName = "region";
    word regionPath;

    // Search for the appropriate blockMesh dictionary....
    fileName dictPath;

    // Check if the dictionary is specified on the command-line
    if (args.optionFound("dict"))
    {
        dictPath = args["dict"];

        dictPath =
        (
            isDir(dictPath)
          ? dictPath/dictName
          : dictPath
        );
    }
    // Check if dictionary is present in the constant directory
    else if
    (
        exists
        (
            runTime.path()/runTime.constant()
           /regionPath/polyMesh::meshSubDir/dictName
        )
    )
    {
        dictPath =
            runTime.constant()
           /regionPath/polyMesh::meshSubDir/dictName;
    }
    // Otherwise assume the dictionary is present in the system directory
    else
    {
        dictPath = runTime.system()/regionPath/dictName;
    }

    
    vector v = vector(0,0,1);

    args.optionReadIfPresent("vector",v);  

//    blockMesh::verbose(true);

    multiMesh** meshObject = new multiMesh* [nThreads];

    word regionNameLocal = regionName+name(myRank);

    dictPath = runTime.system()/"blockMeshDict";

    IOobject meshDictIO
    (
	dictPath,
	runTime,
	IOobject::MUST_READ,
	IOobject::NO_WRITE,
	false
    );

    IOdictionary meshDict(meshDictIO);

    ITstream it = (meshDict.lookup("vertices")); 

    scalar startX = it[2].number();
    scalar startY = it[3].number();
    scalar startZ = it[4].number();
    scalar endX = it[7].number();
    scalar endY = it[13].number();
    scalar endZ = it[24].number();

    scalar zExtent = abs(endZ - startZ);
    scalar yExtent = abs(endY - startY);
    scalar xExtent = abs(endX - startX);

    vector x = vector(1,0,0);
    vector y = vector(0,1,0);
    vector z = vector(0,0,1);

    scalar dotx = v & x;
    scalar doty = v & y;
    scalar dotz = v & z;

    int numX, numY, numZ;
  
    numX = max(nThreads*dotx,1);
    numY = max(nThreads*doty,1);
    numZ = max(nThreads*dotz,1);


    // Need to find out if there is grading in the direction of mesh cutting
    // Need to be able to handle such grading

    // First open the dictionary and add spaces to the end of the entries

    ITstream isx = meshDict.lookup("blocks");
    DynamicList<token> tx;
    tx.setSize(isx.size()+100);
//    tx.setSize(isx.size());
    for(int i=0;i<isx.size();i++)
	tx[i] = isx[i];
    for(int i=isx.size();i<tx.size();i++)
	tx[i] = token(token::SPACE);
    ITstream isx2("dummy",tx);
    meshDict.remove("blocks");
    meshDict.add("blocks",isx2);

    // Now lookup for the block information to determine cell count & grading
	
    ITstream is = meshDict.lookup("blocks");
    label perNX1, perNY1, perNZ1;

    labelList cellCounts;
    cellCounts.setSize(nThreads);

    for(int i=0;i<nThreads;i++)
	cellCounts[i] = 0;

    if(is[13].type()==6)
    {
	label nx = is[13].labelToken();
	label ny = is[14].labelToken(); 
	label nz = is[15].labelToken();

	label perNX = nx/numX;
	label perNY = ny/numY;
	label perNZ = nz/numZ;

	perNX1 = 0;
	perNY1 = 0;
	perNZ1 = 0;

	MPI_Allreduce(&perNX,&perNX1,1, MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&perNY,&perNY1,1, MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&perNZ,&perNZ1,1, MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	perNX1 = (nx - perNX1);
	perNY1 = (ny - perNY1);
	perNZ1 = (nz - perNZ1);

	if (perNX1!=0)
		perNX = myRank<perNX1 ? perNX + 1 : perNX;
	
	if (perNY1!=0)
		perNY = myRank<perNY1 ? perNY + 1 : perNY;

	if (perNZ1!=0)
		perNZ = myRank<perNZ1 ? perNZ + 1 : perNZ;
			
	is[13] = perNX;
    	is[14] = perNY;
    	is[15] = perNZ;

	if (dotx)
		cellCounts[myRank] = perNX;
	else if (doty)
		cellCounts[myRank] = perNY;
	else
		cellCounts[myRank] = perNZ;
	
	Pstream::gatherList(cellCounts);
	Pstream::scatterList(cellCounts);
    }
    else if(is[13].type()==1)
    {
	label nx = is[14].labelToken();
	label ny = is[15].labelToken(); 
	label nz = is[16].labelToken();

	label perNX = nx/numX;
	label perNY = ny/numY;
	label perNZ = nz/numZ;

	perNX1 = 0;
	perNY1 = 0;
	perNZ1 = 0;

	MPI_Allreduce(&perNX,&perNX1,1, MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&perNY,&perNY1,1, MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&perNZ,&perNZ1,1, MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	perNX1 = (nx - perNX1);
	perNY1 = (ny - perNY1);
	perNZ1 = (nz - perNZ1);

	if (perNX1!=0)
		perNX = myRank<perNX1 ? perNX + 1 : perNX;
	
	if (perNY1!=0)
		perNY = myRank<perNY1 ? perNY + 1 : perNY;

	if (perNZ1!=0)
		perNZ = myRank<perNZ1 ? perNZ + 1 : perNZ;

	is[14] = perNX;
    	is[15] = perNY;
    	is[16] = perNZ;

	if (dotx)
		cellCounts[myRank] = perNX;
	else if (doty)
		cellCounts[myRank] = perNY;
	else
		cellCounts[myRank] = perNZ;

	Pstream::gatherList(cellCounts);
	Pstream::scatterList(cellCounts);
    }

    // We use (( to determine where to locally replace dictionary
    // Not the most elegant but is a hack that works
    // Currently we can support only one block 

    int val = is.size()-100;
    label locstart;
    if(is[12].isWord())
	locstart = 20;
    else
	locstart = 19;
    SubList<token> sbk(is,val-locstart-2,locstart);
    labelList L1 = findIndices(sbk,token(token::BEGIN_LIST));
    labelList L2 = findIndices(sbk,token(token::END_LIST));

    int numBracket = 0;
    DynamicList<label> bracketStartLocation;
    for(int i=0;i<L1.size()-1;i++)
	if((L1[i+1] - L1[i]) == 1)
	{
		numBracket++;
		bracketStartLocation.setSize(bracketStartLocation.size() + 1);
		bracketStartLocation[bracketStartLocation.size() -1] = i;
	}

    DynamicList<label> bracketEndLocation;
//    Info << "sbk = " << sbk << endl; 
//    Info << "List L2 is: " << L2 << endl;
//    Info << "List L1 is: " << L1 << endl;
    for(int i=1;i<L2.size();i++)
	if((L2[i] - L2[i-1]) == 1)
	{
		bracketEndLocation.setSize(bracketEndLocation.size() + 1);
		bracketEndLocation[bracketEndLocation.size() -1] = i;
	}

    Info << "Number of such brackets = " << numBracket << endl;
    Info << "Bracket start locations = " << bracketStartLocation << endl;
    Info << "Bracket end locations = " << bracketEndLocation << endl;

//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Finalize();
//    std::exit(0);

    DynamicList<scalar> zvalues; 
    zvalues.setSize(1);

    // Assign where in the list of (( we are supposed to start
    // If direction of cutting along X, then start at 0
    // If direction of cutting along Y, then start at 1
    // If direction of cutting along Z (or left default), then start at 2
    int iStart, iEnd, iFinal;
    if (dotx)
    {
	iStart = bracketStartLocation[0];
//	iEnd = bracketEndLocation[0];
	iEnd = bracketStartLocation[1];
    	zvalues[0] = startX;
    }
    else if (doty)
    {	
	iStart = bracketStartLocation[1];
//	iEnd = bracketEndLocation[1];
	iEnd = bracketStartLocation[2];
	zvalues[0] = startY;
    }
    else
    {
	iStart = bracketStartLocation[2];
//	iEnd = bracketEndLocation[2];
	iEnd = L1.size();
	zvalues[0] = startZ;
    }

    iFinal = bracketEndLocation[2];

    Info << "iStart = " << iStart << endl;    
    Info << "iEnd = " << iEnd << endl;
    Info << "iFinal = " << iFinal << endl;

//    for (int i=iStart+1;i<L1.size();i++)
    for(int i=iStart+1;i<iEnd;i++)
    {
	scalar locDt = sbk[L1[i]+1].number();
	label locNumCells = sbk[L1[i]+2].labelToken();
	scalar locRatio = sbk[L1[i]+3].number();

	scalar locr = (::fabs(locRatio - 1.0) < 1e-4) ? 1.0 : ::pow(locRatio,(1.0/(locNumCells - 1)));
	label oldSize = zvalues.size();
	scalar locA = (::fabs(locRatio - 1.0) < 1e-4) ? locDt/locNumCells : (locDt*(1 - locr))/(1 - pow(locr,locNumCells));

	label startPow = 0;

	zvalues.setSize(zvalues.size()+locNumCells);

	for (int j=oldSize;j<zvalues.size();j++)
	{
		zvalues[j] = zvalues[j-1] + locA*pow(locr,startPow++);
	}
    }

    scalar localStartX, localStartY, localStartZ, localEndX, localEndY, localEndZ;

    Info << "Initializing the local variables " << endl;

    label skip = 0;

    if (dotx)
    {
    	localStartZ = startZ;
    	localStartY = startY;

    	localEndZ = endZ;
    	localEndY = endY;

    	for(int i=myRank-1;i>=0;i--)
    	{
		skip = skip + cellCounts[i];
    	}

    	localStartX = zvalues[skip];
    	localEndX = zvalues[skip+cellCounts[myRank]];
    }
    else if (doty)
    {
    	localStartX = startX;
    	localStartZ = startZ;

    	localEndX = endX;
    	localEndZ = endZ;

    	for(int i=myRank-1;i>=0;i--)
    	{
		skip = skip + cellCounts[i];
    	}

    	localStartY = zvalues[skip];
    	localEndY = zvalues[skip+cellCounts[myRank]];
    }
    else
    {
    	localStartX = startX;
    	localStartY = startY;

    	localEndX = endX;
    	localEndY = endY;

    	for(int i=myRank-1;i>=0;i--)
    	{
		skip = skip + cellCounts[i];
    	}

    	localStartZ = zvalues[skip];
    	localEndZ = zvalues[skip+cellCounts[myRank]];
    }

    SubList<scalar> zLocal(zvalues,cellCounts[myRank]+1,skip);

    DynamicList<scalar> localLength;
    DynamicList<scalar> localFullR;
    DynamicList<label> localLengthCellCount;

    scalar locLen = 0.0;
    scalar locLenR = 1.0;
    label locCellCount = 0;

    scalar startr = (zLocal[2] - zLocal[1])/(zLocal[1] - zLocal[0]);
    locLen+=::fabs(zLocal[2] - zLocal[0]);
    locLenR*=startr;
    locCellCount+=2;
    
    int indCnt = 3;

    while(indCnt < zLocal.size())
    {
	scalar evalr = (zLocal[indCnt] - zLocal[indCnt-1])/(zLocal[indCnt-1] - zLocal[indCnt-2]);
	
	if (::fabs(startr - evalr) < 2e-4)
	{
		locLen+=zLocal[indCnt] - zLocal[indCnt-1];
		locLenR*=evalr;
		locCellCount++;
	}
	else
	{
		localLength.setSize(localLength.size() + 1);
		localFullR.setSize(localFullR.size() + 1);
		localLengthCellCount.setSize(localLengthCellCount.size() + 1);

		localLength[localLength.size()-1] = locLen;
		localFullR[localFullR.size()-1] = locLenR;
		localLengthCellCount[localLengthCellCount.size()-1] = locCellCount;

//		locLen = ::fabs(zLocal[indCnt] - zLocal[indCnt-1]);
//		startr = evalr;
//		locLenR = evalr;
//		locCellCount = 1;

		locLen = ::fabs(zLocal[indCnt + 1] - zLocal[indCnt - 1]);
		startr = (zLocal[indCnt + 1] - zLocal[indCnt])/(zLocal[indCnt] - zLocal[indCnt-1]);
		locLenR = startr;
		locCellCount = 2;
		indCnt++;
	}
	indCnt++;
    }

    localLength.setSize(localLength.size() + 1);
    localFullR.setSize(localFullR.size() + 1);
    localLengthCellCount.setSize(localLengthCellCount.size() + 1);

    localLength[localLength.size()-1] = locLen;
    localFullR[localFullR.size()-1] = locLenR;
    localLengthCellCount[localLengthCellCount.size()-1] = locCellCount;

//    Pout << "localLength = " << localLength << endl;
//    Pout << "localFullR = " << localFullR << endl;
//    Pout << "localLengthCellCount = " << localLengthCellCount << endl;

    label locCounterAdd = L1[iStart] + locstart;
    is[locCounterAdd++] = token::BEGIN_LIST;
    
    for(int i=0;i<localLength.size();i++)
    {
	is[locCounterAdd++] = token::BEGIN_LIST;
	is[locCounterAdd++] = token(localLength[i]);
	is[locCounterAdd++] = token(localLengthCellCount[i]);
	is[locCounterAdd++] = token(localFullR[i]);
	is[locCounterAdd++] = token::END_LIST;
    }

    if ((dotx) || (doty))
    {
    	is[locCounterAdd++] = token::END_LIST;

    	Info << "Adding final bits to stream " << endl;

    	for(int i=L1[iEnd];i<L2[iFinal];i++)
    	{
		is[locCounterAdd++] = token(sbk[i]);
   	}
    }
    
    is[locCounterAdd++] = token::END_LIST;
    is[locCounterAdd++] = token::END_LIST;
    is[locCounterAdd++] = token::END_LIST;
    is[locCounterAdd++] = token::END_STATEMENT;
    
    if (locCounterAdd < is.size())
    {
	for(int i=locCounterAdd;i<is.size();i++)
		is[i] = token::SPACE;
    }


//    Info << "Finally stream is: " << is << endl;

//    MPI_Barrier(MPI_COMM_WORLD);
//    MPI_Finalize();
//    std::exit(0);

/*    is[locCounterAdd++] = token::END_LIST;
    is[locCounterAdd++] = token::END_LIST;
    is[locCounterAdd++] = token::END_STATEMENT;

    if (locCounterAdd < is.size())
    {
	for(int i=locCounterAdd;i<is.size();i++)
		is[i] = token::SPACE;
    }*/

//    Info << "Amended stream now is = " << is << endl;

    vectorField pts(8);
    
    pts[0] = Vector<double>(localStartX,localStartY,localStartZ);
    pts[1] = Vector<double>(localEndX,localStartY,localStartZ);
    pts[2] = Vector<double>(localEndX,localEndY,localStartZ);
    pts[3] = Vector<double>(localStartX,localEndY,localStartZ);

    pts[4] = Vector<double>(localStartX,localStartY,localEndZ);
    pts[5] = Vector<double>(localEndX,localStartY,localEndZ);
    pts[6] = Vector<double>(localEndX,localEndY,localEndZ);
    pts[7] = Vector<double>(localStartX,localEndY,localEndZ);

    meshDict.remove("vertices");
    meshDict.add("vertices",pts);

  
    meshDict.remove("blocks");
    meshDict.add("blocks",is);
    
    ITstream it2 = (meshDict.lookup("boundary"));  
    

    DynamicList<token> t;
  
    if ((myRank==0) || (myRank==nThreads-1))
    {
	t.setSize(it2.size() + 22);
    }
    else
    {
	t.setSize(it2.size() + 44);
    }

    if (myRank==0)
    {
	labelList l1 = findIndices(it2,token(token::BEGIN_BLOCK));
	labelList l2 = findIndices(it2,token(token::END_BLOCK));
	label indStart = l1[l1.size()-1] - 1;
	label indEnd = l2[l2.size()-1];

	SubList<token> s(it2,indEnd - indStart,indStart);

	labelList l3 = findIndices(s,token(token::BEGIN_LIST));
	labelList l4 = findIndices(s,token(token::END_LIST));
	label subIndStart = l3[0];
	label subIndEnd = l4[l4.size()-1];

    	int pp=0;
    	for (int i=0;i<(indStart + subIndStart);i++)
		t[pp++] = it2[i];

	t[pp++] = token(token::BEGIN_LIST);
	t[pp++] = token(token::END_LIST);
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(token::END_BLOCK);

	t[pp++] = token(word("proc0to1"));
	t[pp++] = token(token::BEGIN_BLOCK);
	t[pp++] = token(word("type"));
	t[pp++] = token(word("processor"));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("transform"));
	t[pp++] = token(word("unknown"));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("myProcNo"));
	t[pp++] = token(label(0));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("neighbProcNo"));
	t[pp++] = token(label(1));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("matchTolerance"));
	t[pp++] = token(scalar(0.0001));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("faces"));
	
	for(int i=subIndStart;i<=subIndEnd;i++)
		t[pp++] = s[i];
	
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(token::END_BLOCK);
	t[pp++] = token(token::END_LIST);

	ITstream it3("temp",t);
	meshDict.remove("boundary");
	meshDict.add("boundary",it3);
    }
    else if (myRank==nThreads-1)
    {
	labelList l1 = findIndices(it2,token(token::BEGIN_BLOCK));
	labelList l2 = findIndices(it2,token(token::END_BLOCK));
	label indStart = l1[l1.size()-2] - 1;
	label indEnd = l2[l2.size()-2];

	SubList<token> s(it2,indEnd - indStart,indStart);

	labelList l3 = findIndices(s,token(token::BEGIN_LIST));
	labelList l4 = findIndices(s,token(token::END_LIST));
	label subIndStart = l3[0];
	label subIndEnd = l4[l4.size()-1];

    	int pp=0;
    	for (int i=0;i<(indStart + subIndStart);i++)
		t[pp++] = it2[i];

	t[pp++] = token(token::BEGIN_LIST);
	t[pp++] = token(token::END_LIST);
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(token::END_BLOCK);

	label nextIndStart = l1[l1.size()-1] - 1;
	label nextIndEnd = l2[l2.size()-1];

    	for (int i=nextIndStart;i<nextIndEnd;i++)
		t[pp++] = it2[i];

	t[pp++] = token(token::END_BLOCK);

	word procID = "proc" + name(nThreads-1) + "to" + name(nThreads-2);
	t[pp++] = token(procID);
	t[pp++] = token(token::BEGIN_BLOCK);
	t[pp++] = token(word("type"));
	t[pp++] = token(word("processor"));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("transform"));
	t[pp++] = token(word("unknown"));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("myProcNo"));
	t[pp++] = token(label(nThreads-1));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("neighbProcNo"));
	t[pp++] = token(label(nThreads-2));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("matchTolerance"));
	t[pp++] = token(scalar(0.0001));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("faces"));
	
	for(int i=subIndStart;i<=subIndEnd;i++)
		t[pp++] = s[i];
	
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(token::END_BLOCK);
	t[pp++] = token(token::END_LIST);

	ITstream it3("temp",t);
	meshDict.remove("boundary");
	meshDict.add("boundary",it3);
    }
    else
    {
	labelList l1 = findIndices(it2,token(token::BEGIN_BLOCK));
	labelList l2 = findIndices(it2,token(token::END_BLOCK));
	label indStart = l1[l1.size()-2] - 1;
	label indEnd = l2[l2.size()-2];

	SubList<token> s(it2,indEnd - indStart,indStart);

	labelList l3 = findIndices(s,token(token::BEGIN_LIST));
	labelList l4 = findIndices(s,token(token::END_LIST));
	label subIndStart = l3[0];
	label subIndEnd = l4[l4.size()-1];

    	int pp=0;
    	for (int i=0;i<(indStart + subIndStart);i++)
		t[pp++] = it2[i];

	t[pp++] = token(token::BEGIN_LIST);
	t[pp++] = token(token::END_LIST);
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(token::END_BLOCK);

	label nextIndStart = l1[l1.size()-1] - 1;
	label nextIndEnd = l2[l2.size()-1];

	SubList<token> sNext(it2,nextIndEnd - nextIndStart,nextIndStart);

	labelList l5 = findIndices(sNext,token(token::BEGIN_LIST));
	labelList l6 = findIndices(sNext,token(token::END_LIST));
	label nextSubIndStart = l5[0];
	label nextSubIndEnd = l6[l6.size()-1];

    	for (int i=nextIndStart;i<(nextIndStart + nextSubIndStart);i++)
		t[pp++] = it2[i];

	t[pp++] = token(token::BEGIN_LIST);
	t[pp++] = token(token::END_LIST);
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(token::END_BLOCK);

	word procID1 = "proc" + name(myRank) + "to" + name(myRank-1);
	t[pp++] = token(procID1);
	t[pp++] = token(token::BEGIN_BLOCK);
	t[pp++] = token(word("type"));
	t[pp++] = token(word("processor"));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("transform"));
	t[pp++] = token(word("unknown"));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("myProcNo"));
	t[pp++] = token(label(myRank));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("neighbProcNo"));
	t[pp++] = token(label(myRank-1));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("matchTolerance"));
	t[pp++] = token(scalar(0.0001));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("faces"));
	
	for(int i=subIndStart;i<=subIndEnd;i++)
		t[pp++] = s[i];

	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(token::END_BLOCK);

	word procID2 = "proc" + name(myRank) + "to" + name(myRank+1);
	t[pp++] = token(procID2);
	t[pp++] = token(token::BEGIN_BLOCK);
	t[pp++] = token(word("type"));
	t[pp++] = token(word("processor"));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("transform"));
	t[pp++] = token(word("unknown"));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("myProcNo"));
	t[pp++] = token(label(myRank));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("neighbProcNo"));
	t[pp++] = token(label(myRank+1));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("matchTolerance"));
	t[pp++] = token(scalar(0.0001));
	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(word("faces"));
	
	for(int i=nextSubIndStart;i<=nextSubIndEnd;i++)
		t[pp++] = sNext[i];

	t[pp++] = token(token::END_STATEMENT);
	t[pp++] = token(token::END_BLOCK);
	t[pp++] = token(token::END_LIST);
	
	ITstream it3("temp",t);
	meshDict.remove("boundary");
	meshDict.add("boundary",it3);
    }

    meshObject[myRank] = new multiMesh(regionNameLocal, meshDict,runTime, myRank);
    meshObject[myRank]->createMesh();

    MPI_Barrier(MPI_COMM_WORLD);

    double t2 = MPI_Wtime();
    Info << "Total time taken for meshing = " << (t2 - t1) << "s" << endl;

    meshObject[myRank]->writeMesh();

    MPI_Barrier(MPI_COMM_WORLD);

    double t3 = MPI_Wtime();
    Info << "Time taken for writing mesh = " << (t3 - t2) << "s" << endl;

    label perCellCount = meshObject[myRank]->mesh->nCells();
    label perFaceCount = meshObject[myRank]->mesh->nFaces();
    label perPointCount = meshObject[myRank]->mesh->nPoints();

    label totalCellCount = 0;
    label totalFaceCount = 0;
    label totalPointCount = 0;

    MPI_Allreduce(&perCellCount,&totalCellCount,1, MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&perFaceCount,&totalFaceCount,1, MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&perPointCount,&totalPointCount,1, MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

    Info<< "----------------" << nl
        << "Mesh Information" << nl
	<< "----------------" << nl
        << "  " << "boundingBox: " << boundBox(vector(startX,startY,startZ),vector(endX,endY,endZ)) << nl
	<< "  " << "nPoints: " << totalPointCount << nl
	<< "  " << "nCells: " << totalCellCount << nl
	<< "  " << "nFaces: " << totalFaceCount << endl;
    
    delete meshObject[myRank];    

    delete [] meshObject ;

    MPI_Barrier(MPI_COMM_WORLD);
    
    double t4 = MPI_Wtime();

    Info << "Total time taken for the entire process = " << (t4 - t1) << "s" << endl;
    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
