#include "resetVtxRemapCmd.h"
#include "resetVtxRemapNode.h"

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
	MFnPlugin plugin( obj, PLUGIN_COMPANY, "4.0", "Any");

	status = plugin.registerCommand( "resetVtxRemap", resetVtxRemap::creator );
	if (!status) {
		status.perror("registerCommand");
		return status;
	}

	status = plugin.registerNode( "resetVtxRemapNode",
								  resetVtxRemapNode::id,
								  resetVtxRemapNode::creator,
								  resetVtxRemapNode::initialize );
	if (!status) {
		status.perror("registerNode");
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

	status = plugin.deregisterCommand( "resetVtxRemap" );
	if (!status) {
		status.perror("deregisterCommand");
		return status;
	}

	status = plugin.deregisterNode( resetVtxRemapNode::id );
	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

	return status;
}
