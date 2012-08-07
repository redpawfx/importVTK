//
// Copyright (C) Remik Ziemlinski
// 
// # $Id: pluginMain.cpp,v 1.2 2008/04/23 20:48:47 rsz Exp $

#include "importvtk.h"

#include <maya/MFnPlugin.h>

MStatus initializePlugin( MObject obj )
//
//	Description:
//		this method is called when the plug-in is loaded into Maya.  It 
//		registers all of the services that this plug-in provides with 
//		Maya.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{ 
	MStatus   status;
	MFnPlugin plugin( obj, "Remik Ziemlinski", "7.0", "Any");

	status = plugin.registerCommand( "importvtk", 
										importvtk::creator,
										importvtk::newSyntax);
	if (!status) {
		status.perror("registerCommand");
		return status;
	}

	return status;
}

MStatus uninitializePlugin( MObject obj )
//
//	Description:
//		this method is called when the plug-in is unloaded from Maya. It 
//		deregisters all of the services that it was providing.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{
	MStatus   status;
	MFnPlugin plugin( obj );

	status = plugin.deregisterCommand( "importvtk" );
	if (!status) {
		status.perror("deregisterCommand");
		return status;
	}

	return status;
}
