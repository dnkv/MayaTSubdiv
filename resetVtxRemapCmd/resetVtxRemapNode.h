#ifndef _resetVtxRemapNode
#define _resetVtxRemapNode
// 
// File: resetVtxRemapNode.h
//
// Dependency Graph Node: resetVtxRemapNode
//
// Authors: Lonnie Li, Jeyprakash Michaelraj
//

#include "polyModifierNode.h"
#include "resetVtxRemapFty.h"

// General Includes
//
#include <maya/MTypeId.h>
 
class resetVtxRemapNode : public polyModifierNode
{
public:
						resetVtxRemapNode();
	virtual				~resetVtxRemapNode(); 

	virtual MStatus		compute( const MPlug& plug, MDataBlock& data );

	static  void*		creator();
	static  MStatus		initialize();

public:


	// The typeid is a unique 32bit indentifier that describes this node.
	// It is used to save and retrieve nodes of this type from the binary
	// file format.  If it is not unique, it will cause file IO problems.
	//
	static	MTypeId		id;

	resetVtxRemapFty			fresetVtxRemapFactory;
};

#endif
