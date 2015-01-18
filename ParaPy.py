## Purpose: Automate some of the post-processing tasks in Paraview like taking snapshots
## Usage: Can be run from command line from the case folder
## Usage: paraview --script=<Location of this file>
## Status: Experimental
## Notes: File called "sliceList" that provides the normal vector and origin
## Notes: No warranty on results. Use at your own risk/discretion
## Notes: Code is free. Appreciate feedback/acknowledging when using it
## Notes: Tested with Python 2.7 and Paraview 4.2
## Created by: Venugopalan Raghavan

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

import os as os

cwd = os.getcwd().replace("\\","/")

iF = cwd + "/sliceList"
iFR = open(iF,"r")
currline = iFR.readline()

normalsList = list()
originsList = list()

while(currline!=''):
	parts = currline.strip("\n").split()
	normalsList.append(parts[0])
	originsList.append(parts[1])
	currline = iFR.readline()
	
iFR.close()

print normalsList

caseFileName = cwd + "/paraview.foam"

paraview_foam = OpenFOAMReader(FileName=caseFileName,CaseType=1)
arr = paraview_foam.TimestepValues

RenderView1 = GetRenderView()
a1_p_PVLookupTable = GetLookupTableForArray( "p", 1, RGBPoints=[0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.498039215686275, 0.498039215686275, 0.498039215686275], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

a1_p_PiecewiseFunction = CreatePiecewiseFunction( Points=[-3256.28002929688, 0.0, 0.5, 0.0, 1429.8701171875, 1.0, 0.5, 0.0] )

a1_p_PVLookupTable.ScalarOpacityFunction = a1_p_PiecewiseFunction

for ii in xrange(0,len(normalsList)):	
	Slice1 = Slice(paraview_foam)
	Slice1.SliceOffsetValues = [0.0]
	Slice1.SliceType.Origin = [0.0, 0.0, 1.5]
	Slice1.SliceType = "Plane"
	Slice1.SliceType.Normal = [0.0, 1.0, 0.0]
	
	if normalsList[ii]=="[1,0,0]":
		Slice1.SliceType.Normal=[1,0,0]
		Slice1.SliceType.Origin=[float(originsList[ii]),0,0]
		outputFile=cwd+"Velocity_X_"+originsList[ii]+"_m.png"
		RenderView1.CameraViewUp = [0.0, 0.0, 1.0]
		RenderView1.CameraPosition = [-9.971049470053016, 0.0, 1.5]
	elif normalsList[ii]=="[0,1,0]":
		Slice1.SliceType.Normal=[0,1,0]
		Slice1.SliceType.Origin=[0,float(originsList[ii]),0]
		outputFile=cwd+"Velocity_Y_"+originsList[ii]+"_m.png"
		RenderView1.CameraViewUp = [0.0, 0.0, 1.0]
		RenderView1.CameraPosition = [0.0, -9.971049470053016, 1.5]
	elif normalsList[ii]=="[0,0,1]":
		Slice1.SliceType.Normal=[0,0,1]
		Slice1.SliceType.Origin=[0,0,float(originsList[ii])]
		outputFile=cwd+"Velocity_Z_"+originsList[ii]+"_m.png"
		RenderView1.CameraViewUp = [0.0, 1.0, 0.0]
		RenderView1.CameraPosition = [1.5, 0.0, -9.971049470053016]
	else:
		p=str(normalsList[ii])
		q=p.strip("[")
		r=q.strip("]")
		s=r.split(",")
		Slice1.SliceType.Normal=[float(s[0]),float(s[1]),float(s[2])]
		Slice1.SliceType.Origin=[0,0,0]
		outputFile=cwd+"Velocity_Misc+"+originsList[ii]+"_m.png"
		RenderView1.CameraViewUp = [0.0, 0.0, 1.0]
		RenderView1.CameraPosition = [0.0, -9.971049470053016, 1.5]

	#RenderView1.CameraPosition = [0.0, 0.0, 13.656858288132726]
	RenderView1.CameraFocalPoint = [0.0, 0.0, 1.5]
	RenderView1.CameraClippingRange = [9.0502897052514, 16.096711162454717]
	RenderView1.CameraParallelScale = 3.146426453581177

	active_objects.source.SMProxy.InvokeEvent('UserEvent', 'ShowWidget')

	DataRepresentation2 = Show()
	DataRepresentation2.EdgeColor = [0.0, 0.0, 0.5000076295109483]
	DataRepresentation2.SelectionPointFieldDataArrayName = 'p'
	DataRepresentation2.SelectionCellFieldDataArrayName = 'p'
	DataRepresentation2.ColorArrayName = ('POINT_DATA', 'p')
	DataRepresentation2.LookupTable = a1_p_PVLookupTable
	DataRepresentation2.ScaleFactor = 0.4199999809265137

	RenderView1.CameraClippingRange = [8.41925257370347, 16.891056831009024]

	active_objects.source.SMProxy.InvokeEvent('UserEvent', 'HideWidget')

	qa = Slice1.PointData.GetArray('U').GetRange()[1]

	a3_U_PVLookupTable = GetLookupTableForArray( "U", 3, RGBPoints=[0.0, 0.0, 0.0, 1.0, float(qa), 1.0, 0.0, 0.0], VectorMode='Magnitude', NanColor=[0.498039215686275, 0.498039215686275, 0.498039215686275], ColorSpace='HSV', ScalarRangeInitialized=1.0, AllowDuplicateScalars=1 )

	a3_U_PiecewiseFunction = CreatePiecewiseFunction( Points=[-3256.28002929688, 0.0, 0.5, 0.0, 1429.8701171875, 1.0, 0.5, 0.0] )

	ScalarBarWidgetRepresentation1 = CreateScalarBar(ComponentTitle='Magnitude', Title='U', Enabled=1, LabelFontSize=10, LookupTable=a3_U_PVLookupTable, TitleFontSize=10)

	GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)
	
	RenderView1.CacheKey = 1000.0
	RenderView1.CameraClippingRange = [9.871338975352487, 10.12061521210381]
	RenderView1.ViewTime = 1000.0
	RenderView1.UseCache = 0
	RenderView1.CameraFocalPoint = [0.0, 0.0, 1.5]
	RenderView1.CameraParallelScale = 2.5806975025091172

	DataRepresentation2.ColorArrayName = ('POINT_DATA', 'U')
	DataRepresentation2.LookupTable = a3_U_PVLookupTable

	a3_U_PVLookupTable.ScalarOpacityFunction = a3_U_PiecewiseFunction

	WriteImage(outputFile)
	Show()
	Render()
	Hide()
