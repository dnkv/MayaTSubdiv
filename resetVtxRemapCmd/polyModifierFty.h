#ifndef _polyModifierFty
#define _polyModifierFty

#include <maya/MStatus.h>

class polyModifierFty
{
public:
						polyModifierFty();
	virtual				~polyModifierFty();

	// Pure virtual doIt()
	//
	virtual MStatus		doIt() = 0;
};

#endif
