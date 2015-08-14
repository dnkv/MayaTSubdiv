#ifndef HDS_H
#define HDS_H

#include "kvc.h"
#include <complex>


// This is created automatically - DO NOT MODIFY!
//----------------------------------------------------------------------------

struct HDS_R: public KVCStruct
{
    // per face
    KVCRealArray<uint64_t>  nFV;
    KVCIndexArray<uint64_t> F2H;
    
    // per halfedge
    KVCIndexArray<uint64_t> H2F;
    KVCIndexArray<uint64_t> next, prev, tip, twin;
    
    KVCRealArray<bool>      T;
    KVCRealArray<double>    itv;
    
    // per vertex
    KVCIndexArray<uint64_t> V2H;

    KVCRealArray<double>     V;
    
    size_t nF() { return std::length(this->nFV); }
    size_t nHE() { return std::length(this->tip); }
    size_t nIHE() { return std::length(this->H2F); }
    size_t nV() { return std::length(this->V2H); }

    HDS_R()
    {
        std::string  k[] = { "nFV", "F2H", "H2F", "next", "prev", "tip", "twin", "T", "itv", "V2H"};
        KVCObject   *v[] = { &nFV,  &F2H,  &H2F,  &next,  &prev,  &tip,  &twin,  &T,  &itv,  &V2H };
        KVC(k, v, sizeof(v)/sizeof(KVCObject *));
    }
};

struct HDS_R_TOPO: public KVCStruct
{
    // per face
    KVCRealArray<uint64_t> nFV;
    KVCIndexArray<uint64_t> F2H;
    
    // per halfedge
    KVCIndexArray<uint64_t> H2F;
    KVCIndexArray<uint64_t> next, prev, tip, twin;
    
    KVCRealArray<bool>     T;
    KVCRealArray<double>   itv;
    
    // per vertex
    KVCIndexArray<uint64_t> V2H;
    
    size_t nF() { return length(this->nFV); }
    size_t nHE() { return length(this->tip); }
    size_t nIHE() { return length(this->H2F); }
    size_t nV() { return std::length(this->V2H); }

    HDS_R_TOPO()
    {
        std::string  k[] = { "nFV", "F2H", "H2F", "next", "prev", "tip", "twin", "T", "itv", "V2H"};
        KVCObject   *v[] = { &nFV,  &F2H,  &H2F,  &next,  &prev,  &tip,  &twin,  &T,  &itv,  &V2H };
        KVC(k, v, sizeof(v)/sizeof(KVCObject *));
    }
};


struct HDS: public KVCStruct
{
    // per face
    KVCRealArray<uint64_t>  nFV;
    KVCIndexArray<uint64_t> F2H;
    
    // per halfedge
    KVCIndexArray<uint64_t> H2F;
    KVCIndexArray<uint64_t> next, prev, tip, tipK, twin;
    
    KVCRealArray<bool>      T;
    KVCRealArray<uint64_t>  pole;
    KVCRealArray<double>    itv;
    KVCRealArray<double>    val;
    
    // per vertex
    KVCIndexArray<uint64_t> V2H;
    KVCRealArray<uint64_t>      corner;
    
    KVCRealArray<double>     V;
    KVCRealArray<double>     VK, VFK;
//    KVCRealArray<double>     UV;

    size_t nF() { return std::length(this->nFV); }
    size_t nHE() { return std::length(this->tip); }
    size_t nIHE() { return std::length(this->H2F); }
    size_t nV() { return std::length(this->V2H); }
    
    HDS()
    {
        std::string  k[] = { "nFV", "F2H", "H2F", "next", "prev", "tip", "tipK", "twin", "T", "itv", "V2H", "V", "VK", "VFK", "pole", "corner", "val"};
        KVCObject   *v[] = { &nFV,  &F2H,  &H2F,  &next,  &prev,  &tip,  &tipK,  &twin,  &T,  &itv,  &V2H,  &V , &VK , &VFK , &pole, &corner, &val};
        KVC(k, v, sizeof(v)/sizeof(KVCObject *));
    }
    
};

/*
struct HDS: public KVCStruct
{
    // per face
    KVCMutableRealArray<uint64_t> nFV;
    KVCMutableIndexArray<uint64_t> F2H;
    
    // per halfedge
    KVCMutableIndexArray<uint64_t> H2F;
    KVCMutableIndexArray<uint64_t> next, prev, tip, twin;
    
    KVCMutableRealArray<bool>     T;
    KVCMutableRealArray<double>   itv;
    
    // per vertex
    KVCMutableIndexArray<uint64_t> V2H;
    
    HDS()
    {
        std::string  k[] = { "nFV", "F2H", "H2F", "next", "prev", "tip", "twin", "T", "itv", "V2H"};
        KVCObject   *v[] = { &nFV,  &F2H,  &H2F,  &next,  &prev,  &tip,  &twin,  &T,  &itv,  &V2H };
        KVC(k, v, sizeof(v)/sizeof(KVCObject *));
    }
};
*/

//----------------------------------------------------------------------------

bool get_face_vertices(HDS_R &hds, uint64_t fI, KVCIndexArray<uint64_t> &v);
bool get_face_halfedges(HDS_R &hds, uint64_t fI, KVCIndexArray<uint64_t> &he);
bool get_halfedge_onering(HDS_R &hds, uint64_t hI, KVCIndexArray<uint64_t> &v);
bool get_halfedge_halfedges(HDS_R &hds, uint64_t hI, KVCIndexArray<uint64_t> &he);
bool get_vertex_onering(HDS_R &hds, uint64_t vI, KVCIndexArray<uint64_t> &v);
bool get_vertex_halfedges(HDS_R &hds, uint64_t vI, KVCIndexArray<uint64_t> &he);


bool finalize_HDS(HDS &HDS);
bool subdivide_HDS(HDS &hds, HDS &hds_sd);

#endif