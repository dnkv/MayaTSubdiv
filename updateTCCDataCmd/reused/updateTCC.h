//
//  updateTCC.h
//  MatlabKVCHDS
//
//  Created by Denis Kovacs on 9/15/12.
//  Copyright 2012 Denis Kovacs. All rights reserved.
//

#ifndef MatlabKVCHDS_updateTCC_h
#define MatlabKVCHDS_updateTCC_h

#include "kvc.h"

struct TCCNodeData: public KVCStruct
{
    KVCMutableRealArray<uint64_t>  nFV;
    KVCMutableIndexArray<uint64_t> tip;
    KVCMutableRealArray<double>    V;
    
    KVCMutableRealArray<uint64_t>  pole;
    KVCMutableRealArray<uint64_t>  corner;
    
    KVCMutableRealArray<bool>      T;
    KVCMutableRealArray<uint64_t>  eqc;
    KVCMutableRealArray<double>    itv;
    KVCMutableRealArray<uint64_t>  err;
    
    KVCMutableIndexArray<uint64_t> selHE;
    KVCMutableIndexArray<uint64_t> selV;

    
    TCCNodeData()
    {
        std::string  k[] = { "nFV", "tip", "V", "pole", "T", "eqc", "itv", "err", "selHE", "selV", "corner"};
        KVCObject   *v[] = { &nFV,  &tip,  &V,  &pole,  &T,  &eqc,  &itv,  &err,  &selHE , &selV,  &corner};
        KVC(k, v, sizeof(v)/sizeof(KVCObject *));
    }
};

void updateTCC(TCCNodeData &nd);

#endif
