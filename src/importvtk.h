#ifndef _importvtk
#define _importvtk
//
// Copyright (C) Remik Ziemlinski
// 
// MEL Command: importvtk
//
// $Id: importvtk.h,v 1.3 2008/05/26 13:07:10 rsz Exp $

#include <maya/MIOStream.h>
#include <maya/MPxCommand.h>
#include <maya/MArgList.h>
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>
#include <maya/MArgParser.h>
#include <maya/MFnParticleSystem.h>
#include <maya/MPointArray.h>
#include <maya/MObject.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatArray.h>
#include <maya/MFloatPointArray.h>
#include <maya/MFloatVectorArray.h>
#include <maya/MFnFluid.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MColorArray.h>
#include <maya/MFnNurbsCurve.h>

#include "vtkStructuredGrid.h"
#include "vtkDataSetReader.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkStructuredPoints.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkLookupTable.h"

#define kDebugFlag "-d"
#define kDebugFlagLong "-debug"
#define kFilenameFlag		"-f"
#define kFilenameFlagLong	"-filename"
#define kNameFlag			"-n"
#define kNameFlagLong		"-name"
#define kNormalsFlag "-N"
#define kNormalsFlagLong "-normals"
#define kReferenceFlag		"-r"
#define kReferenceFlagLong	"-reference"
#define kScalarsAttFlag		"-s"
#define kScalarsAttFlagLong	"-scalarsAttribute"
#define kVersionFlag "-v"
#define kVersionFlagLong "-version"

#define VERSION "1.1.0"

#define CHECKRESULT(stat,msg)     \
	if ( MS::kSuccess != stat ) {		\
		displayError( msg );					\
	}

class MArgList;

class importvtk : public MPxCommand
{

public:
	importvtk();
	virtual	~importvtk();

	MStatus	doIt( const MArgList& );
	MStatus	redoIt();
	MStatus	undoIt();
	bool isUndoable() const;

	static void* creator();
	static MSyntax newSyntax();
	MStatus parseArgs(const MArgList&);
	
	MStatus GetReferenceStringPointer(const MString&, void**, MString&);
	MStatus vtkStructuredGridToParticles(vtkStructuredGrid*);
	MStatus createParticles(MPointArray&, MDoubleArray*);
	MStatus vtkStructuredPointsToFluid(vtkStructuredPoints*);
	MStatus vtkPolyDataImport(vtkPolyData*);
	MStatus vtkPolyDataToMesh(vtkPolyData*);
	MStatus vtkPolyDataToLines(vtkPolyData*);
	MStatus doImport();

private:
	// Store the data you will need to undo the command here
	//	

	// Reference of VTK object in memory.
	MString reference;
	// Name of input VTK file.
	MString filename;
	// Name for the particle system.
	MString name;
	// Name of particle system node attribute that will have scalars attached to it.
	MString scalarsAttributeName;
	// Debug level. If 0, nothing will print.
	int debuglevel;
	// Let user decide if normals should be imported, which can be slow with current setVertexNormal implementation.
	int importnormals;
};

#endif
