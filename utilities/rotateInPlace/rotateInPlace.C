/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    rotateInPlace

Description
    rotates a mesh about its centre.

Usage
    rotateInPlace -angle '(theta_x theta_y theta_z)' 

    Angles theta_x, theta_y and theta_z to be specified in degrees 

    Optionally a point set can be specified with the -pointSet flag
    When -pointSet is specified, only the pointSet will be rotated

    Useful when a mesh needs to be rotated before initializing/starting runs

Original
    Adapted from rotateMesh

Version History
    v0: 20 Feb 2017 - initial version for OF 2.4.x
    v1: 10 Feb 2020 - code clean-up, comments & port to OF 5.x

Adapted by
    Venugopalan Raghavan (IHPC)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IStringStream.H"
#include "mathematicalConstants.H"

#include "pointSet.H"

using namespace Foam;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "angle",
        "vector",
        "rotate by specified <vector> - eg, '(1 0 0)'"
    );

    argList::addOption
    (
        "pointSet",
        "word",
        "name of set to be rotated"
    );

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
    word regionName = polyMesh::defaultRegion;
    fileName meshDir;

    meshDir = polyMesh::meshSubDir;

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshDir, "points"),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    vector v;

    args.optionReadIfPresent("angle",v);

    Info<< "Rotating points by" << nl
        << "    angleX   " << v.x() << nl
        << "    angleY " << v.y() << nl
        << "    angleZ  " << v.z() << nl;

     // Convert to radians
     v *= pi/180.0;

     quaternion R = quaternion(vector(1, 0, 0), v.x());
     R *= quaternion(vector(0, 1, 0), v.y());
     R *= quaternion(vector(0, 0, 1), v.z());

     word nameOfSet = polyMesh::defaultRegion;     
     args.optionReadIfPresent("pointSet", nameOfSet);

     vector t1;
     vector t2;

     Info<< "Rotating all points by quaternion " << R << endl;
     vectorField temp = transform(R,points);

Info << "Initial Domain center = " << 0.5*(max(points.component(vector::X)) + min(points.component(vector::X))) << " " << 0.5*(max(points.component(vector::Y)) + min(points.component(vector::Y))) << " " << 0.5*(max(points.component(vector::Z)) + min(points.component(vector::Z))) << endl;

     if (nameOfSet!=polyMesh::defaultRegion)	
     {
	double preMinX,preMinY,preMinZ,preMaxX,preMaxY,preMaxZ;
	double postMinX,postMinY,postMinZ,postMaxX,postMaxY,postMaxZ;

	pointSet setToUse(mesh,nameOfSet);
	labelList pointLabel = setToUse.toc();

	preMinX=preMinY=preMinZ=postMinX=postMinY=postMinZ=ROOTVGREAT;
	preMaxX=preMaxY=preMaxZ=postMaxX=postMaxY=postMaxZ=-ROOTVGREAT;
	
	for(int i=0;i<setToUse.size();i++)
	{
		int index = pointLabel[i];

		double X = points[index].component(vector::X);
		double Y = points[index].component(vector::Y);
		double Z = points[index].component(vector::Z);

		preMinX = preMinX < X ? preMinX : X;
		preMinY = preMinY < Y ? preMinY : Y;
		preMinZ = preMinZ < Z ? preMinZ : Z;

		preMaxX = preMaxX > X ? preMaxX : X;
		preMaxY = preMaxY > Y ? preMaxY : Y;
		preMaxZ = preMaxZ > Z ? preMaxZ : Z;

		X = temp[index].component(vector::X);
		Y = temp[index].component(vector::Y);
		Z = temp[index].component(vector::Z);

		postMinX = postMinX < X ? postMinX : X;
		postMinY = postMinY < Y ? postMinY : Y;
		postMinZ = postMinZ < Z ? postMinZ : Z;

		postMaxX = postMaxX > X ? postMaxX : X;
		postMaxY = postMaxY > Y ? postMaxY : Y;
		postMaxZ = postMaxZ > Z ? postMaxZ : Z;

		points[index].component(vector::X) = temp[index].component(vector::X);
		points[index].component(vector::Y) = temp[index].component(vector::Y);
		points[index].component(vector::Z) = temp[index].component(vector::Z);
	}


     	t1.x() = 0.5*(preMinX + preMaxX);
    	t1.y() = 0.5*(preMinY + preMaxY);
     	t1.z() = 0.5*(preMinZ + preMaxZ);

	Info << "Original center = " << t1.x() << " " << t1.y() << " " << t1.z() << endl;

	Info<< "Restricting the transformation to pointSet only " << endl;

     	t2.x() = 0.5*(postMinX + postMaxX);
    	t2.y() = 0.5*(postMinY + postMaxY);
     	t2.z() = 0.5*(postMinZ + postMaxZ);

     	Info << "New center = " << t2.x() << " " << t2.y() << " " << t2.z() << endl;

     	Info << "Translation Vector = " << t1.x()-t2.x() << " " << t1.y()-t2.y() << " " << t1.z()-t2.z() << endl;

	for(int i=0;i<setToUse.size();i++)
		points[pointLabel[i]]+=(t1-t2);

     }

     else
     {
	t1.x() = 0.5*(max(points.component(vector::X)) + min(points.component(vector::X)));
     	t1.y() = 0.5*(max(points.component(vector::Y)) + min(points.component(vector::Y)));
     	t1.z() = 0.5*(max(points.component(vector::Z)) + min(points.component(vector::Z)));

     	Info << "Original center = " << t1.x() << " " << t1.y() << " " << t1.z() << endl;

	for(int i=0;i<points.size();i++)
	{
		int index = i;

		points[index].component(vector::X) = temp[index].component(vector::X);
		points[index].component(vector::Y) = temp[index].component(vector::Y);
		points[index].component(vector::Z) = temp[index].component(vector::Z);
	}

     	t2.x() = 0.5*(max(points.component(vector::X)) + min(points.component(vector::X)));
     	t2.y() = 0.5*(max(points.component(vector::Y)) + min(points.component(vector::Y)));
     	t2.z() = 0.5*(max(points.component(vector::Z)) + min(points.component(vector::Z)));

	Info << "New center = " << t2.x() << " " << t2.y() << " " << t2.z() << endl;

     	Info << "Translation Vector = " << t1.x()-t2.x() << " " << t1.y()-t2.y() << " " << t1.z()-t2.z() << endl;

     	points+=(t1-t2);
     }

    Info << "Final Domain center = " << 0.5*(max(points.component(vector::X)) + min(points.component(vector::X))) << " " << 0.5*(max(points.component(vector::Y)) + min(points.component(vector::Y))) << " " << 0.5*(max(points.component(vector::Z)) + min(points.component(vector::Z))) << endl;

	// Set the precision of the points data to 10
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

    Info<< "Writing points into directory " << points.path() << nl << endl;
    points.write();

    return 0;
}


// ************************************************************************* //
