//
// Copyright (C) Remik Ziemlinski
// 
// MEL Command: importvtk
//
// $Id: importvtk.cpp,v 1.4 2008/05/26 13:06:48 rsz Exp $
// Todo:
// - Add selection of existing particle system for redefinition. Low priority.

#include "importvtk.h"
#include <maya/MGlobal.h>

MStatus importvtk::doIt( const MArgList& args)
//
//	Description:
//		implements the MEL importvtk command.
//
//	Arguments:
//		args - the argument list that was passes to the command from MEL
//
//	Return Value:
//		MS::kSuccess - command succeeded
//		MS::kFailure - command failed (returning this value will cause the 
//                     MEL script that is being run to terminate unless the
//                     error is caught using a "catch" statement.
//
{
	MStatus stat = MS::kSuccess;

	stat = parseArgs(args);
	if (stat != MS::kSuccess)
		return stat;

	// Typically, the doIt() method only collects the infomation required
	// to do/undo the action and then stores it in class members.  The 
	// redo method is then called to do the actuall work.  This prevents
	// code duplication.
	//
	return redoIt();
}

MStatus importvtk::redoIt()
//
//	Description:
//		implements redo for the MEL importvtk command. 
//
//		This method is called when the user has undone a command of this type
//		and then redoes it.  No arguments are passed in as all of the necessary
//		information is cached by the doIt method.
//
//	Return Value:
//		MS::kSuccess - command succeeded
//		MS::kFailure - redoIt failed.  this is a serious problem that will
//                     likely cause the undo queue to be purged
//
{
	return doImport();
}
MStatus importvtk::doImport()
{
	MStatus status = MS::kFailure;
	int datasetType = -999;
	vtkStructuredGrid* sgrid = NULL;
	vtkStructuredPoints* spts = NULL;
	vtkPolyData* poly = NULL;
	vtkRectilinearGrid* rgrid = NULL;
	vtkUnstructuredGrid* ugrid = NULL;

	if (reference != "") {
		void* voidptr;
		MString classname;
		status = GetReferenceStringPointer(reference, &voidptr, classname);

		if (status == MS::kSuccess) {
			// Find the object type.
			if (classname == "vtkStructuredGrid") {
				datasetType = VTK_STRUCTURED_GRID;
				sgrid = (vtkStructuredGrid*)voidptr;
			} else if (classname == "vtkStructuredPoints") {
				datasetType = VTK_STRUCTURED_POINTS;
				spts = (vtkStructuredPoints*)voidptr;
			} else if (classname == "vtkPolyData") {
				datasetType = VTK_POLY_DATA;
				poly = (vtkPolyData*)voidptr;
			} else if (classname == "vtkRectilinearGrid") {
				datasetType = VTK_RECTILINEAR_GRID;
				rgrid = (vtkRectilinearGrid*)voidptr;
			} else if (classname == "vtkUnstructuredGrid") {
				datasetType = VTK_UNSTRUCTURED_GRID;
				ugrid = (vtkUnstructuredGrid*)voidptr;
			}
		}
	} else if (filename != "") {
		// Read the VTK file.
		vtkDataSetReader* r = vtkDataSetReader::New();
		r->SetFileName(filename.asChar());
		r->Update();
		// Check data set type.
		datasetType = r->ReadOutputType();

		sgrid = r->GetStructuredGridOutput();
		spts = r->GetStructuredPointsOutput();
		poly = r->GetPolyDataOutput();
		ugrid = r->GetUnstructuredGridOutput();
		rgrid = r->GetRectilinearGridOutput();

		status = MS::kSuccess;
	}

	if (status == MS::kSuccess) {
		switch (datasetType)
			{
			case VTK_STRUCTURED_GRID:
				if (debuglevel) 
					MGlobal::displayInfo( "\tImporting VTK_STRUCTURED_GRID.\n" );
				
				status = vtkStructuredGridToParticles(sgrid);
				break;
			case VTK_STRUCTURED_POINTS: // Wrapper for vtkImageData.
				if (debuglevel) 
					MGlobal::displayInfo( "\tImporting VTK_STRUCTURED_POINTS.\n" );
				
				status = vtkStructuredPointsToFluid(spts);
				break;
			case VTK_POLY_DATA:			
				if (debuglevel) 
					MGlobal::displayInfo( "\tImporting VTK_POLY_DATA.\n" );
				
				status = vtkPolyDataImport(poly);
				break;
			case VTK_RECTILINEAR_GRID:
			case VTK_UNSTRUCTURED_GRID:
			default:			
				MString msg = "importvtk: Object type currently unsupported.";
				displayError( msg );	
				status = MS::kFailure;		
				break;			
			}
	}

	// Since this class is derived off of MPxCommand, you can use the 
	// inherited methods to return values and set error messages
	//
	setResult( status );
	
	return status;
}

MStatus importvtk::undoIt()
//
//	Description:
//		implements undo for the MEL importvtk command.  
//
//		This method is called to undo a previous command of this type.  The 
//		system should be returned to the exact state that it was it previous 
//		to this command being executed.  That includes the selection state.
//
//	Return Value:
//		MS::kSuccess - command succeeded
//		MS::kFailure - redoIt failed.  this is a serious problem that will
//                     likely cause the undo queue to be purged
//
{

	// You can also display information to the command window via MGlobal
	//
    MGlobal::displayInfo( "importvtk command undone!\n" );

	return MS::kSuccess;
}

void* importvtk::creator()
//
//	Description:
//		this method exists to give Maya a way to create new objects
//      of this type. 
//
//	Return Value:
//		a new object of this type
//
{
	return new importvtk();
}

importvtk::importvtk() : reference(""), filename(""), name(""), scalarsAttributeName(""), debuglevel(0), importnormals(0)
//
//	Description:
//		importvtk constructor
//
{}

importvtk::~importvtk()
//
//	Description:
//		importvtk destructor
//
{
}

bool importvtk::isUndoable() const
//
//	Description:
//		this method tells Maya this command is undoable.  It is added to the 
//		undo queue if it is.  
//      This will not support undo to avoid caching large amounts of 
//      extra particle memory.
//
//	Return Value:
//		true if this command is undoable.
//
{
	return false;
}

MSyntax importvtk::newSyntax()
{
    MSyntax syntax;
		
    syntax.addFlag(kReferenceFlag, kReferenceFlagLong, MSyntax::kString);
    syntax.addFlag(kFilenameFlag, kFilenameFlagLong, MSyntax::kString);
    syntax.addFlag(kNameFlag, kNameFlagLong, MSyntax::kString);
    syntax.addFlag(kScalarsAttFlag, kScalarsAttFlagLong, MSyntax::kString); 
    syntax.addFlag(kDebugFlag, kDebugFlagLong, MSyntax::kLong); 
    syntax.addFlag(kVersionFlag, kVersionFlagLong, MSyntax::kNoArg); 
    syntax.addFlag(kNormalsFlag, kNormalsFlagLong, MSyntax::kNoArg);
 
    return syntax;
}
MStatus importvtk::parseArgs( const MArgList& args )
{    
	MStatus status = MS::kSuccess;	
	int displayversion = 0;

	//MArgDatabase results were screwy, so doing bruteforce parsing.
	for ( unsigned int i = 0; i < args.length(); i++ )
		{
			if ( ( MString( kFilenameFlag ) == args.asString( i, &status )
						 && MS::kSuccess == status ) || 
					 ( MString( kFilenameFlagLong ) == args.asString( i, &status )
						 && MS::kSuccess == status ) )
				{
					filename = args.asString( ++i, &status );
				}	else if ( ( MString( kNameFlag ) == args.asString( i, &status )
										 && MS::kSuccess == status ) || 
									 ( MString( kNameFlagLong ) == args.asString( i, &status )
										 && MS::kSuccess == status ) )
				{
					name = args.asString( ++i, &status );
				} else if ( ( MString( kScalarsAttFlag ) == args.asString( i, &status )
										 && MS::kSuccess == status ) || 
									 ( MString( kScalarsAttFlagLong ) == args.asString( i, &status )
										 && MS::kSuccess == status ) )
				{
					scalarsAttributeName = args.asString( ++i, &status );
				} else if ( ( MString( kDebugFlag ) == args.asString( i, &status )
										 && MS::kSuccess == status ) || 
									 ( MString( kDebugFlagLong ) == args.asString( i, &status )
										 && MS::kSuccess == status ) )
				{
					debuglevel = args.asInt( ++i, &status );
				} else if ( ( MString( kVersionFlag ) == args.asString( i, &status )
										 && MS::kSuccess == status ) || 
									 ( MString( kVersionFlagLong ) == args.asString( i, &status )
										 && MS::kSuccess == status ) )
				{
					MString msg = "importvtk version ";
					displayInfo(msg + VERSION);		
					displayversion = 1;
				} else if ( ( MString( kNormalsFlag ) == args.asString( i, &status )
											&& MS::kSuccess == status ) || 
										( MString( kNormalsFlagLong ) == args.asString( i, &status )
											&& MS::kSuccess == status ) )
				{
					importnormals = 1;
				} else if ( ( MString( kReferenceFlag ) == args.asString( i, &status )
											&& MS::kSuccess == status ) || 
										( MString( kReferenceFlagLong ) == args.asString( i, &status )
											&& MS::kSuccess == status ) )
				{
					reference = args.asString( ++i, &status );
				} else
				{
					MString msg = "Invalid flag: ";
					msg += args.asString( i );
					displayError( msg );
					return MS::kFailure;
				}
		}

	if ((!displayversion) && (filename == "") && (reference == ""))
		{
			MString msg = "Required filename or reference missing. Use the ";
			msg += kFilenameFlag;
			msg += " or ";
			msg += kReferenceFlag;
			msg += " flag.";
			displayError( msg );
			return MS::kFailure;
		}

	return status;
}
/*
Converts a pointer string like '_03929f92_p_vtkObject' into a void* pointer that can be cast to an object pointer.

In:
text - The encoded reference string.
Out:
ptr - A void pointer to the object reference.
classname - The VTK object class name of the instance being pointed to (alternative to ptr->GetClassName()).
*/
MStatus importvtk::GetReferenceStringPointer(const MString & text, void** ptr, MString & classname)
{
  int n;
	char buf[256];
  
  n = sscanf(text.asChar(), "_%lx_p_%s", (long*)ptr, buf);
  
  if (n != 2) {
		*ptr = NULL;
		return MS::kFailure;
	}

	classname.set(buf);

  return MS::kSuccess;
}

MStatus importvtk::vtkStructuredGridToParticles(vtkStructuredGrid* g)
{
	MStatus status = MS::kSuccess;
	// MPointArray needs an array of 4-tuples.
	double (*pts4)[4];
	unsigned int npts = ((vtkPointSet*)g)->GetNumberOfPoints();
	vtkPoints* pts = ((vtkPointSet*)g)->GetPoints();

	double *scalars = NULL;
	MDoubleArray *mscalars = NULL;

	pts4 = new double[npts][4];

	for(unsigned int i=0; i < npts; ++i)
		{
			// Should only populate the first 3 elements.
			// Can't use points memory directly because MPoint is a 4tuple, not 3tuple.
			pts->GetPoint(i, &pts4[i][0]);			
		}
	
	if (scalarsAttributeName != "")
		{	
			// Access the scalar memory directly instead of copying values to new array.
			vtkDataArray *vtkscalars = ((vtkDataSet*)g)->GetPointData()->GetScalars();
			switch (vtkscalars->GetDataType())
				{
				case VTK_DOUBLE:
					mscalars = new MDoubleArray(((vtkDoubleArray*)vtkscalars)->GetPointer(0), npts);
					break;
				case VTK_FLOAT:
					mscalars = new MDoubleArray(((vtkFloatArray*)vtkscalars)->GetPointer(0), npts);
					break;
				default:
					displayInfo("Only double and float type scalar arrays supported.\n");
					break;
				}
		}

	MPointArray mpts(pts4, npts);
	status = createParticles(mpts, mscalars);

	delete [] pts4;
	if (mscalars != NULL)	delete mscalars;

	return status;
}
MStatus importvtk::createParticles(MPointArray& pa, MDoubleArray* scalars)
{
 	MStatus stat = MS::kSuccess;
	
	/* Need to use dummy particle system to create an MObject, and then have
     another particle system do the emit. I got burned using a single 
     particle system for everything resulting in no particles. 
		 I realized this in the plugin example "devkit/particleSystemInfoCmd".
	*/
	MFnParticleSystem dummy;
	MObject pobj = dummy.create(&stat);
	CHECKRESULT(stat,"MFnParticleSystem::create(status) failed!");

	MFnParticleSystem ps( pobj, &stat );
	CHECKRESULT(stat,"MFnParticleSystem::MFnParticleSystem(MObject,status) failed!");

	stat = ps.emit(pa);
	CHECKRESULT(stat,"MFnParticleSystem::emit(MPointArray) failed!");

	if (scalars == NULL)
		/* Nothing more to do. */
		return stat;

	MFnDoubleArrayData dData;
	MObject dObj = dData.create(*scalars, &stat);
	CHECKRESULT(stat,"MFnDoubleArrayData::create() failed!");
	
	// Don't test stat, because stat!= success if test returns false.
	bool hasattr;
	MObject attr;

	// Create a new attribute.
	MFnTypedAttribute fnattr;	
	attr = fnattr.create(scalarsAttributeName, 
											 scalarsAttributeName.substring(0,2),
											 MFnData::kDoubleArray,
											 dObj, &stat);
	CHECKRESULT(stat,"MFnTypedAttribute::create() failed!");

	stat = ps.addAttribute(attr, MFnDependencyNode::kLocalDynamicAttr);
	CHECKRESULT(stat,"MFnParticleSystem::addAttribute() failed!");

	hasattr = ps.hasAttribute(scalarsAttributeName, &stat);
	CHECKRESULT(stat,"MFnParticleSystem doesn't have new attribute!");
	
	ps.setPerParticleAttribute(scalarsAttributeName, *scalars, &stat);
	CHECKRESULT(stat,"MFnParticleSystem::setPerParticleAttribute() failed!");

	stat = ps.saveInitialState();
	CHECKRESULT(stat,"MFnParticleSystem::saveInitialState() failed!");

	return stat;
}
MStatus importvtk::vtkStructuredPointsToFluid(vtkStructuredPoints* spts)
{
	MStatus stat = MS::kSuccess;

	if (spts == NULL)
		return stat;

	int dims[3];
	spts->GetDimensions(dims);
	double dxyz[3];
	spts->GetSpacing(dxyz);

	MFnFluid fn;
	MObject node = fn.create3D(dims[0], dims[1], dims[2],
														 dxyz[0], dxyz[1], dxyz[2],
														 MObject::kNullObj, &stat);
	CHECKRESULT(stat,"MFnFluid::create3D() failed!");

	if (scalarsAttributeName == "")
		return stat;

	vtkDataArray *vtkscalars = ((vtkDataSet*)spts)->GetPointData()->GetScalars();
	/* Pointer to fluid data, such as density. */
	float* dest = NULL;
	double* dsrc = NULL;
	float* fsrc = NULL;
	int npts = dims[0]*dims[1]*dims[2];

	stat = fn.setVelocityMode(MFnFluid::kZero, MFnFluid::kConstant);
	CHECKRESULT(stat,"MFnFluid::setDensityMode() failed!");

	if (scalarsAttributeName == "density")
		{
			// Note: Second argument (kConstant) is ignored in StaticGrid context.
			stat = fn.setDensityMode(MFnFluid::kStaticGrid, MFnFluid::kConstant);
			CHECKRESULT(stat,"MFnFluid::setDensityMode() failed!");

			dest = fn.density(&stat);
			CHECKRESULT(stat,"MFnFluid::density() failed!");
			// Exit because destination pointer is invalid so nothing to assign scalars to.
			if (stat != MS::kSuccess)
				return stat;
		}

	switch (vtkscalars->GetDataType())
		{
		case VTK_DOUBLE:
			dsrc = ((vtkDoubleArray*)vtkscalars)->GetPointer(0);
			for(unsigned int i=0; i<npts; ++i)
				dest[i] = (float)dsrc[i];
			
			break;
		case VTK_FLOAT:
			fsrc = ((vtkFloatArray*)vtkscalars)->GetPointer(0);
			for(unsigned int i=0; i<npts; ++i)
				dest[i] = fsrc[i];

			break;
		default:
			displayInfo("Only double and float type scalar arrays supported.\n");
			break;
		}

	stat = fn.updateGrid();
	CHECKRESULT(stat,"MFnFluid::updateGrid() failed!");			
				
	return stat;
}
/*
	Imports the various geometries part of a vtkPolyData.
	For now only imports polygons.
 */
MStatus importvtk::vtkPolyDataImport(vtkPolyData* poly)
{
	MStatus stat = MS::kSuccess;

	if (poly == NULL) {
		if (debuglevel) {
			displayInfo("importvtk::vtkPolyDataImport: poly == NULL\n");
		}
		return stat;
	}

	if ( (poly->GetNumberOfPolys() > 0) ||
			 (poly->GetNumberOfStrips() > 0) )
		{
			if (debuglevel) {
				char buf[255];
				sprintf(buf, "vtkPolyData has %d polys, %d strips\n", 
								poly->GetNumberOfPolys(), poly->GetNumberOfStrips());

				MString msg(buf);
				displayInfo(msg);
			}

			stat = vtkPolyDataToMesh(poly);			
			CHECKRESULT(stat,"vtkPolyDataToMesh() failed!");
		}

	if (poly->GetNumberOfLines() > 0)
		stat = vtkPolyDataToLines(poly);

	return stat;
}
/*
	Only imports vertices (colormapped if scalars present),
	normals and texture coords.
	Todo: 
	- Cells (data,texture,normals).
 */
MStatus importvtk::vtkPolyDataToMesh(vtkPolyData* poly)
{
	MStatus stat = MS::kSuccess;

	if (poly == NULL)
		return stat;

	vtkCellArray *cells = NULL;

	if (poly->GetNumberOfPolys())
		cells = poly->GetPolys();
	else if (poly->GetNumberOfStrips())
		cells = poly->GetStrips();

	// MPointArray needs an array of 4-tuples.
	double (*pts4)[4];
	unsigned int npts = poly->GetNumberOfPoints();
	unsigned int ncells = cells->GetNumberOfCells();
	vtkPoints* pts = poly->GetPoints();
	int *polycounts = new int[ncells];
	int *polyconnects = new int[cells->GetNumberOfConnectivityEntries()];
	pts4 = new double[npts][4];
	vtkIdTypeArray* vtkcelldata = cells->GetData();
	int *celldata = vtkcelldata->GetPointer(0);
	unsigned int i, imax;

	for(i=0; i < npts; ++i)
		{
			// Should only populate the first 3 elements.
			// Can't use points memory directly because MPoint is a 4tuple, not 3tuple.
			pts->GetPoint(i, &pts4[i][0]);			
		}

	if (debuglevel) 
		MGlobal::displayInfo( "\tProcessed points.\n" );

	imax = vtkcelldata->GetNumberOfTuples();
	// Poly index.
	unsigned int k = 0;
	// Poly connect index.
	unsigned int m = 0;
	// Copy out the values in the celldata array: (n1,id1,id2,...,n2,id1,id2,...)
	for(i=0; i < imax;)
		{
			polycounts[k++] = celldata[i++];
			
			for(unsigned int c=0; c < polycounts[k-1]; ++i, ++c, ++m)
				{
					polyconnects[m] = celldata[i];
				}			
		}
	
	if (debuglevel) 
		MGlobal::displayInfo( "\tProcessed vertex connectivity.\n" );

	vtkDataArray* vtkuv = ((vtkDataSet*)poly)->GetPointData()->GetTCoords();
	float *us = NULL;
	float *vs = NULL;
	unsigned int nuvs = 0;

	// Not every mesh will have texture coords.
	if (vtkuv != NULL)
		{
			if (debuglevel) 
				MGlobal::displayInfo( "\tProcessing texture coords.\n" );

			nuvs = npts;
			us = new float[npts];
			vs = new float[npts];

			for(unsigned int i=0; i < npts; ++i)
				{
					us[i] = (float)vtkuv->GetComponent(i, 0);
					vs[i] = (float)vtkuv->GetComponent(i, 1);
				}

			if (debuglevel) 
				MGlobal::displayInfo( "\t\tDone.\n" );
		}

	MFloatArray uArray(us, nuvs);
	MFloatArray vArray(vs, nuvs);

	MFloatPointArray ptsarray(pts4, npts);
	MIntArray polycountarray(polycounts, k);
	MIntArray polyconnectsarray(polyconnects, m);

	if (debuglevel) 
		MGlobal::displayInfo( "\tCreating Mesh object.\n" );

	MFnMesh fn;
	MObject obj = fn.create(npts,
													k,
													ptsarray,
													polycountarray,
													polyconnectsarray,
													uArray, vArray,
													MObject::kNullObj,
													&stat);
	CHECKRESULT(stat,"MFnMesh::create() failed!");			
	
	// We must explicitly map the texture coords to mesh vertices. 
	if (nuvs > 0)	{
		if (debuglevel) 
			MGlobal::displayInfo( "\tAssigning texture coords.\n" );

		MString uvSet("map1"); 
		stat = fn.assignUVs(polycountarray, polyconnectsarray, &uvSet);
		CHECKRESULT(stat,"MFnMesh::assignUVs() failed!");			
	}

	delete [] pts4;
	if (us != NULL) delete [] us;
	if (vs != NULL)	delete [] vs;

	// Normals ///////////////////////////////////////////////////	
	vtkDataArray *vtknormals = ((vtkDataSet*)poly)->GetPointData()->GetNormals();

	if ((vtknormals != NULL) && importnormals){
		if (debuglevel) {
			MString msg("\tProcessing ");
			msg += m;
			msg += " normals.";
			MGlobal::displayInfo(msg);
		}

		stat = fn.unlockVertexNormals(polyconnectsarray);
		CHECKRESULT(stat,"MFnMesh::unlockVertexNormals() failed!");			

		double tmp[3];
		double (*normals)[3] = new double[m][3];
		for(i=0; i < m; ++i) { // For each (redundant) polygon vertex.
			// polyconnects[i] gives us the point id.
			// VTK uses the notion of shared vertex normals like Maya.
			// Todo: Support per polygon (exclusive) normals with command options or face normals.
			vtknormals->GetTuple(polyconnects[i], normals[i]);			
		}

		MVectorArray vecarray(normals, m);
		stat = fn.setVertexNormals(vecarray, polyconnectsarray, MSpace::kWorld);
		CHECKRESULT(stat,"MFnMesh::setVertexNormals() failed!");
 
		delete [] normals;
		
		if (debuglevel) 
			MGlobal::displayInfo( "\t\tDone.\n" );
	}

	delete [] polycounts;
	delete [] polyconnects;

	// Scalars ///////////////////////////////////////////////////
	// todo: celldata support.
	vtkDataArray *vtkscalars = ((vtkDataSet*)poly)->GetPointData()->GetScalars();
	
	// Do not proceed if scalars aren't present.
	if (vtkscalars == NULL)
		return stat;

	vtkLookupTable *lut = vtkscalars->GetLookupTable();
	if (lut == NULL)
		return stat;

	/* It seems that VTK file saving doesn't export the correct LUT table range,
     because this print out always shows the default 0,1 range.
		 Instead, use the scalar range.  In the future the range can be specified
     through the MEL arguments or to have an LUT to Maya Color Set feature.
	double tmp[2];
	lut->GetTableRange(tmp);
	displayInfo(MString("lut range min,max = ")+tmp[0]+MString(" ")+tmp[1]);
	*/
	lut->SetRange(vtkscalars->GetRange());

	// MColorArray needs an array of 4-tuples, rgba.
  double (*colors4)[4] = new double[npts][4];	
	int *indices = new int[npts];
	double scalar;

	if (debuglevel) 
		MGlobal::displayInfo( "\tProcessing scalars.\n" );

	for(i=0; i < npts; ++i)
		{
			indices[i] = i;
			scalar = vtkscalars->GetTuple1(i);
			lut->GetColor(scalar, &colors4[i][0]);
			colors4[i][3] = lut->GetOpacity(scalar);
		}	
	
	if (debuglevel) 
		MGlobal::displayInfo( "\tAssigning vertex colors.\n" );

	// Assign to the default color set.
	MColorArray mcolors(colors4, npts);
	MIntArray mindices(indices, npts);

	stat = fn.setVertexColors(mcolors,
														mindices,
														(MDGModifier*)NULL);									 
	CHECKRESULT(stat,"MFnMesh::setVertexColors() failed!");			
	
	delete [] colors4;
	delete [] indices;

	// todo: if cell setFaceColors

	return stat;	
}
MStatus importvtk::vtkPolyDataToLines(vtkPolyData* poly)
{
	MStatus stat = MS::kSuccess;

	if (poly == NULL)
		return stat;

	vtkCellArray *cells = poly->GetLines();
	
	if (cells == NULL)
		return stat;

	int nlines = cells->GetNumberOfCells();
	vtkIdType *celldata = cells->GetPointer();
	unsigned int i, j, k, npts;
	double (*pts4)[4];
	float *knots;

	for(i=0, j=0; i < nlines; ++i)
		{
			npts = celldata[j];
			pts4 = new double[npts][4];
			knots = new float[npts];

			for(k=0; k < npts; ++k)
				{
					poly->GetPoint(celldata[j+k+1], &pts4[k][0]);
					knots[k] = 1;
				}
			
			j += npts + 1 ;

			MFnNurbsCurve fn;
			MObject obj = fn.create(MPointArray(pts4, npts), 
															MDoubleArray(knots, npts),
															1,
															MFnNurbsCurve::kOpen,
															false,
															false,
															MObject::kNullObj,
															&stat);
			CHECKRESULT(stat,"MFnNurbsCurve::create() failed!");			
			
			delete [] pts4;
			delete [] knots;
		}

	return stat;
}
