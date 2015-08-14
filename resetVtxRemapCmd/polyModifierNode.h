#ifndef _polyModifierNode
#define _polyModifierNode

#include <maya/MPxNode.h>

class polyModifierNode : public MPxNode
{
public:
						polyModifierNode();
	virtual				~polyModifierNode(); 

public:

	// There needs to be a MObject handle declared for each attribute that
	// the node will have.  These handles are needed for getting and setting
	// the values later.
	//
	static  MObject		inMesh;
	static  MObject		outMesh;
};

#endif
