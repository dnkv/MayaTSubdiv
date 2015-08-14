#ifndef _updateTCCDataNode
#define _updateTCCDataNode

#include "polyModifierNode.h"
#include "updateTCCDataFty.h"

// General Includes
//
#include <maya/MTypeId.h>
 
class updateTCCDataNode : public polyModifierNode
{
public:
						updateTCCDataNode();
	virtual				~updateTCCDataNode(); 

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static  void*		creator();
	static  MStatus		initialize();

public:

	// There needs to be a MObject handle declared for each attribute that
	// the node will have.  These handles are needed for getting and setting
	// the values later.
	//
	// polyModifierNode has predefined the standard inMesh and outMesh attributes.
	//
	// We define an input attribute for our UV list input
	//
	static  MObject		polyOrder;  // order of old polygons
    static  MObject     cShift;
    static  MObject     vtxRemap;
    static  MObject     delta_F, delta_nFV;
	static  MObject		nV;
    

	// The typeid is a unique 32bit indentifier that describes this node.
	// It is used to save and retrieve nodes of this type from the binary
	// file format.  If it is not unique, it will cause file IO problems.
	//
	static	MTypeId		id;

	updateTCCDataFty			fupdateTCCDataFactory;
};

#endif
