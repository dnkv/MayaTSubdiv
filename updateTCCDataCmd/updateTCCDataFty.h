#ifndef _updateTCCDataFty
#define _updateTCCDataFty

#include "polyModifierFty.h"

// General Includes
//
#include <maya/MObject.h>
#include <maya/MIntArray.h>
#include <maya/MString.h>
#include <maya/MFnMesh.h>

class updateTCCDataFty : public polyModifierFty
{

public:
				updateTCCDataFty();
	virtual		~updateTCCDataFty();

	void		setMesh( MObject mesh ) { fMesh = mesh; }
	void		setPolyOrder( MIntArray polyOrder ) { fPolyOrder = polyOrder; }
	void		setCShift( MIntArray cShift ) { fCShift = cShift; }
	void		setVtxRemap( MIntArray vtxRemap ) { fVtxRemap = vtxRemap; }
	void		setDelta_nFV( MIntArray delta_nFV ) { fDelta_nFV = delta_nFV; }
	void		setDelta_F( MIntArray delta_F ) { fDelta_F = delta_F; }

	void		setnV( size_t nV ) { fnV = nV; }

    MIntArray &getPolyOrder() { return fPolyOrder; }
    MIntArray &getCShift() { return fCShift; }
    MIntArray &getVtxRemap() { return fVtxRemap; }
    MIntArray &getDelta_nFV() { return fDelta_nFV; }
    MIntArray &getDelta_F() { return fDelta_F; }

	unsigned int &getnV() { return fnV; }
    
	// polyModifierFty inherited methods
	//
    
	MStatus		 doIt();

private:
    MStatus      remapMeshData(MFnMesh &meshFn);

	MObject		 fMesh;
    MIntArray    fVtxRemap;            // mapping old to new indices
	MIntArray	 fPolyOrder;           // indices into old polygons
	MIntArray	 fCShift;              // circular shift offset when remapping new to old polygons
    MIntArray    fDelta_F, fDelta_nFV; // delta per-face vertex indices
	unsigned int fnV;                  // number of new vertices
};

#endif
