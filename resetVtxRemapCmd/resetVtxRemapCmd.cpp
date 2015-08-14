#include "resetVtxRemapCmd.h"
#include "resetVtxRemapNode.h"

// Function Sets
//
#include <maya/MFnDependencyNode.h>
#include <maya/MFnMesh.h>
#include <maya/MFnSingleIndexedComponent.h>

// Iterators
//
#include <maya/MItSelectionList.h>
#include <maya/MItMeshPolygon.h>

// General Includes
//
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MPlug.h>

#include <maya/MIOStream.h>

// Status Checking Macro - MCheckStatus (Debugging tool)
//
#define MCheckStatus(status,message)	\
	if( MS::kSuccess != status ) {		\
		cerr << message << "\n";		\
		return status;					\
	}


resetVtxRemap::resetVtxRemap()
//
//	Description:
//		resetVtxRemap constructor
//
{}

resetVtxRemap::~resetVtxRemap()
//
//	Description:
//		resetVtxRemap destructor
//
{}

void* resetVtxRemap::creator()
//
//	Description:
//		this method exists to give Maya a way to create new objects
//      of this type. 
//
//	Return Value:
//		a new object of this type
//
{
	return new resetVtxRemap();
}

bool resetVtxRemap::isUndoable() const
//
//	Description:
//		this method tells Maya this command is undoable.  It is added to the 
//		undo queue if it is.
//
//	Return Value:
//		true if this command is undoable.
//
{
	return true;
}

MStatus resetVtxRemap::doIt( const MArgList& )
//
//	Description:
//		implements the MEL resetVtxRemap command.
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
	MStatus status;

	// Parse the selection list for objects with selected UV components.
	// To simplify things, we only take the first object that we find with
	// selected UVs and operate on that object alone.
	//
	// All other objects are ignored and return warning messages indicating
	// this limitation.
	//
	MSelectionList selList;
	MGlobal::getActiveSelectionList( selList );
	MItSelectionList selListIter( selList );
	selListIter.setFilter( MFn::kMesh );

	// The resetVtxRemap node only accepts a component list input, so we build
	// a component list using MFnComponentListData.
	//
	// MIntArrays could also be passed into the node to represent the uvIds,
	// but are less storage efficient than component lists, since consecutive 
	// components are bundled into a single entry in component lists.
	//
	MFnComponentListData compListFn;
	compListFn.create();
	bool found = false;
	bool foundMultiple = false;

	for( ; !selListIter.isDone(); selListIter.next() )
	{
		MDagPath dagPath;
		selListIter.getDagPath( dagPath );

        if( !found )
        {
            dagPath.extendToShape();
            cout<<dagPath.fullPathName().asChar()<<endl;
            setMeshNode( dagPath );
            found = true;
        }
        else
        {
            foundMultiple = true;
            break;
		}
	}
    
	if( foundMultiple )
	{
		displayWarning("Found more than one object with selected UVs - Only operating on first found object.");
	}

	// Initialize the polyModifierCmd node type - mesh node already set
	//
	setModifierNodeType( resetVtxRemapNode::id );

	if( found )
	{
        status = doModifyPoly();
        
        if( status == MS::kSuccess )
        {
            setResult( "resetVtxRemap command succeeded!" );
        }
        else
        {
            displayError( "resetVtxRemap command failed!" );
        }
    }
	
	return status;
}

MStatus resetVtxRemap::redoIt()
//
//	Description:
//		Implements redo for the MEL resetVtxRemap command. 
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
	MStatus status;

	// Process the polyModifierCmd
	//
	status = redoModifyPoly();

	if( status == MS::kSuccess )
	{
		setResult( "resetVtxRemap command succeeded!" );
	}
	else
	{
		displayError( "resetVtxRemap command failed!" );
	}

	return status;
}

MStatus resetVtxRemap::undoIt()
//
//	Description:
//		implements undo for the MEL resetVtxRemap command.  
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
	MStatus status;

	status = undoModifyPoly();

	if( status == MS::kSuccess )
	{
		setResult( "resetVtxRemap undo succeeded!" );
	}
	else
	{
		setResult( "resetVtxRemap undo failed!" );
	}
    
	return status;
}

MStatus resetVtxRemap::initModifierNode( MObject modifierNode )
{
	MStatus status;

	MFnDependencyNode depNodeFn( modifierNode );

	return status;
}

MStatus resetVtxRemap::directModifier( MObject mesh )
{
	MStatus status;

	fresetVtxRemapFactory.setMesh( mesh );

	status = fresetVtxRemapFactory.doIt();

	return status;
}

