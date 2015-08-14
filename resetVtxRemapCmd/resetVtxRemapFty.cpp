#include "resetVtxRemapFty.h"

resetVtxRemapFty::resetVtxRemapFty()
//
//	Description:
//		resetVtxRemapFty constructor
//
{
}

resetVtxRemapFty::~resetVtxRemapFty()
//
//	Description:
//		resetVtxRemapFty destructor
//
{}

void resetVtxRemapFty::setMesh( MObject mesh )
//
//	Description:
//		Sets the mesh object that this factory will operate on
//
{
	fMesh = mesh;
}
