#include "resetVtxRemapNode.h"

// Function Sets
//
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnComponentListData.h>
#include <maya/MFnSingleIndexedComponent.h>

// General Includes
//
#include <maya/MGlobal.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MIOStream.h>

// Macros
//
#define MCheckStatus(status,message)	\
	if( MStatus::kSuccess != status ) {	\
		cerr << message << "\n";		\
		return status;					\
	}


// Unique Node TypeId
//
MTypeId     resetVtxRemapNode::id( 0x94567 );


resetVtxRemapNode::resetVtxRemapNode()
{}

resetVtxRemapNode::~resetVtxRemapNode()
{}

MStatus resetVtxRemapNode::compute( const MPlug& plug, MDataBlock& data )
//
//	Description:
//		This method computes the value of the given output plug based
//		on the values of the input attributes.
//
//	Arguments:
//		plug - the plug to compute
//		data - object that provides access to the attributes for this node
//
{
	MStatus status = MS::kSuccess;
 
	MDataHandle stateData = data.outputValue( state, &status );
	MCheckStatus( status, "ERROR getting state" );

	// Check for the HasNoEffect/PassThrough flag on the node.
	//
	// (stateData is an enumeration standard in all depend nodes - stored as short)
	// 
	// (0 = Normal)
	// (1 = HasNoEffect/PassThrough)
	// (2 = Blocking)
	// ...
	//
	if( stateData.asShort() == 1 )
	{
		MDataHandle inputData = data.inputValue( inMesh, &status );
		MCheckStatus(status,"ERROR getting inMesh");

		MDataHandle outputData = data.outputValue( outMesh, &status );
		MCheckStatus(status,"ERROR getting outMesh");

		// Simply redirect the inMesh to the outMesh for the PassThrough effect
		//
		outputData.set(inputData.asMesh());
	}
	else
	{
		// Check which output attribute we have been asked to 
		// compute. If this node doesn't know how to compute it, 
		// we must return MS::kUnknownParameter
		// 
		if (plug == outMesh)
		{
			MDataHandle inputData = data.inputValue( inMesh, &status );
			MCheckStatus(status,"ERROR getting inMesh");

			MDataHandle outputData = data.outputValue( outMesh, &status );
			MCheckStatus(status,"ERROR getting outMesh"); 

			// Copy the inMesh to the outMesh, and now you can
			// perform operations in-place on the outMesh
			//
			outputData.set(inputData.asMesh());
			MObject mesh = outputData.asMesh();

			fresetVtxRemapFactory.setMesh( mesh );

			status = fresetVtxRemapFactory.doIt();

			outputData.setClean();
		}
		else
		{
			status = MS::kUnknownParameter;
		}
	}

	return status;
}

void* resetVtxRemapNode::creator()
//
//	Description:
//		this method exists to give Maya a way to create new objects
//      of this type. 
//
//	Return Value:
//		a new object of this type
//
{
	return new resetVtxRemapNode();
}

MStatus resetVtxRemapNode::initialize()
//
//	Description:
//		This method is called to create and initialize all of the attributes
//      and attribute dependencies for this node type.  This is only called 
//		once when the node type is registered with Maya.
//
//	Return Values:
//		MS::kSuccess
//		MS::kFailure
//		
{
	MStatus				status;

	MFnTypedAttribute attrFn;

	inMesh = attrFn.create("inMesh", "im", MFnMeshData::kMesh);
	attrFn.setStorable(true);	// To be stored during file-save

	// Attribute is read-only because it is an output attribute
	//
	outMesh = attrFn.create("outMesh", "om", MFnMeshData::kMesh);
	attrFn.setStorable(false);
	attrFn.setWritable(false);

	status = addAttribute( inMesh );
		if (!status)
		{
			status.perror("addAttribute");
			return status;
		}
	status = addAttribute( outMesh);
		if (!status)
		{
			status.perror("addAttribute");
			return status;
		}

	// Set up a dependency between the input and the output.  This will cause
	// the output to be marked dirty when the input changes.  The output will
	// then be recomputed the next time the value of the output is requested.
	//
	status = attributeAffects( inMesh, outMesh );
		if (!status)
		{
			status.perror("attributeAffects");
			return status;
		}

	return MS::kSuccess;

}
