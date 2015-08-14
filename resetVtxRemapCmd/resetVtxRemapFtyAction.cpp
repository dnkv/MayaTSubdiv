#include "resetVtxRemapFty.h"

// General Includes
//
#include <maya/MGlobal.h>
#include <maya/MIOStream.h>

// Function Sets
//
#include <maya/MFnMesh.h>
#include <maya/MFnSingleIndexedComponent.h>

// Iterators
//
#include <maya/MItMeshPolygon.h>
#define VTXREMAP_BLINDDATA_ID 9005

MStatus initialize_mesh_vertex_blindData(MFnMesh &meshFn)
{
    MStatus stat;
    
    if (meshFn.isBlindDataTypeUsed(VTXREMAP_BLINDDATA_ID, &stat))
    {
        stat = meshFn.clearBlindData(MFn::kMeshVertComponent, VTXREMAP_BLINDDATA_ID, "vtxRemap");  // remove possible old blind data
    }
    else
    {
        MStringArray longNames, shortNames, formatNames;
        longNames.append("vtxRemap");
        shortNames.append("vR");
        formatNames.append("int");
        stat = meshFn.createBlindDataType(VTXREMAP_BLINDDATA_ID, longNames, shortNames, formatNames);
    }
    
    size_t nV = meshFn.numVertices(); 
    
    MIntArray idx(nV, 0), val(nV, 0);
    
    for (size_t k=0; k<nV; k++)
    {
        idx[k] = k;
        val[k] = k;
    }
    
    stat = meshFn.setIntBlindData(idx, MFn::kMeshVertComponent, VTXREMAP_BLINDDATA_ID, "vtxRemap", val);
    
    if (stat != MS::kSuccess)
    {
        stat.perror("");
    }
    return stat;
}


MStatus resetVtxRemapFty::doIt()
//
//	Description:
//		Performs the actual resetVtxRemap operation on the given object and UVs
//
{
	MStatus status = MS::kSuccess;

	MFnMesh meshFn( fMesh );
    
    initialize_mesh_vertex_blindData(meshFn);

	return status;
}
