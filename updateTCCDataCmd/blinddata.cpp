//
//  blinddata.cpp
//  updateTCCDataCmd
//
//  Created by Denis Kovacs on 7/30/12.
//  Copyright 2012 __MyCompanyName__. All rights reserved.
//

#include "blinddata.h"

void print_MIntArray(MIntArray &a)
{
    cout<<"["<<a.length()<<"] ";
    for (size_t k=0; k<a.length(); k++)
    {
        cout<<a[k]<<" ";
    }
    cout<<endl;
}

void print_MDoubleArray(MDoubleArray &a)
{
    cout<<"["<<a.length()<<"] ";
    for (size_t k=0; k<a.length(); k++)
    {
        cout<<a[k]<<" ";
    }
    cout<<endl;
}

bool get_mesh_double_blindData(MFnMesh &meshFn, MFn::Type compType, int typeId, MString bdName, MDoubleArray &array)
{
    if (!meshFn.isBlindDataTypeUsed(typeId))
    {
        return false;
    }
    
    MIntArray idx; MDoubleArray val;
    
    meshFn.getDoubleBlindData(compType, typeId, bdName, idx, val);
    
    for (size_t k=0; k<idx.length(); k++) array[idx[k]] = val[k];
    
    return true;
}

bool get_mesh_int_blindData(MFnMesh &meshFn, MFn::Type compType, int typeId, MString bdName, MIntArray &array)
{
    if (!meshFn.isBlindDataTypeUsed(typeId))
    {
        return false;
    }
    
    MIntArray idx; MIntArray val;
    
    meshFn.getIntBlindData(compType, typeId, bdName, idx, val);
    
    for (size_t k=0; k<idx.length(); k++) array[idx[k]] = val[k];
    
    return true;
}

MStatus set_mesh_double_blindData(MFnMesh &meshFn, MFn::Type compType, int typeId, MString bdName, MDoubleArray &val)
{
    if (meshFn.isBlindDataTypeUsed(typeId))
    {
        meshFn.clearBlindData(compType, typeId, bdName);
    }
    else
    {
        MStringArray longNames, shortNames, formatNames;
        longNames.append(bdName);
        shortNames.append(bdName);
        formatNames.append("double");
        meshFn.createBlindDataType(typeId, longNames, shortNames, formatNames);
    }
    
    MIntArray idx(val.length());
    
    for (size_t k=0; k<val.length(); k++) idx[k] = k;
    
    meshFn.setDoubleBlindData(idx, compType, typeId, bdName, val);
    
    return MS::kSuccess;
}

MStatus set_mesh_int_blindData(MFnMesh &meshFn, MFn::Type compType, int typeId, MString bdName, MIntArray &val)
{
    if (meshFn.isBlindDataTypeUsed(typeId))
    {
        meshFn.clearBlindData(compType, typeId, bdName);
    }
    else
    {
        MStringArray longNames, shortNames, formatNames;
        longNames.append(bdName);
        shortNames.append(bdName);
        formatNames.append("int");
        meshFn.createBlindDataType(typeId, longNames, shortNames, formatNames);
    }
    
    MIntArray idx(val.length());
    
    for (size_t k=0; k<val.length(); k++) idx[k] = k;
    
    meshFn.setIntBlindData(idx, compType, typeId, bdName, val);
    
    return MS::kSuccess;
}
