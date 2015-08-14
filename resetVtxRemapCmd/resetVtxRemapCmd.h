#ifndef _resetVtxRemapCmd
#define _resetVtxRemapCmd

#include "polyModifierCmd.h"
#include "resetVtxRemapFty.h"

// Function Sets
//
#include <maya/MFnComponentListData.h>

// Forward Class Declarations
//
class MArgList;

class resetVtxRemap : public polyModifierCmd
{

public:
	////////////////////
	// Public Methods //
	////////////////////

				resetVtxRemap();
	virtual		~resetVtxRemap();

	static		void* creator();

	bool		isUndoable() const;

	MStatus		doIt( const MArgList& );
	MStatus		redoIt();
	MStatus		undoIt();

	/////////////////////////////
	// polyModifierCmd Methods //
	/////////////////////////////

	MStatus		initModifierNode( MObject modifierNode );
	MStatus		directModifier( MObject mesh );

private:

	// resetVtxRemap Factory
	//
	resetVtxRemapFty				fresetVtxRemapFactory;
};

#endif
