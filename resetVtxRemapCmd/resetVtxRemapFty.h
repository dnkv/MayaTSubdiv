#ifndef _resetVtxRemapFty
#define _resetVtxRemapFty

#include "polyModifierFty.h"

// General Includes
//
#include <maya/MObject.h>
#include <maya/MIntArray.h>
#include <maya/MString.h>

class resetVtxRemapFty : public polyModifierFty
{

public:
				resetVtxRemapFty();
	virtual		~resetVtxRemapFty();

	void		setMesh( MObject mesh );

	// polyModifierFty inherited methods
	//
	MStatus		doIt();

private:

	MObject		fMesh;
};

#endif
