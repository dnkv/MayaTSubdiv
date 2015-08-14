//
//  blinddata.h
//  updateTCCDataCmd
//
//  Created by Denis Kovacs on 7/30/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#ifndef _BLINDDATA_H
#define _BLINDDATA_H

#include <maya/MFnMesh.h>
#include <maya/MFnData.h>
#include <maya/MIntArray.h>
#include <maya/MFloatArray.h>
#include <maya/MDoubleArray.h>


void print_MIntArray(MIntArray &a);


bool get_mesh_int_blindData(MFnMesh &meshFn, MFn::Type compType, int typeId, MString bdName, MIntArray &array);
bool get_mesh_double_blindData(MFnMesh &meshFn, MFn::Type compType, int typeId, MString bdName, MDoubleArray &array);
MStatus set_mesh_int_blindData(MFnMesh &meshFn, MFn::Type compType, int typeId, MString bdName, MIntArray &val);
MStatus set_mesh_double_blindData(MFnMesh &meshFn, MFn::Type compType, int typeId, MString bdName, MDoubleArray &val);


#endif