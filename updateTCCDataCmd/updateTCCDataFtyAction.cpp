#include "updateTCCDataFty.h"

// General Includes
//
#include <maya/MGlobal.h>
#include <maya/MIOStream.h>

// Function Sets
//
#include <maya/MFnMesh.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnData.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFloatPointArray.h>

// Iterators
//
#include <maya/MItMeshPolygon.h>

#include "kvcMatlabVar.h"

#include "blinddata.h"

// YOUR STUFF HERE -----------------------------------------------\

struct HalfedgeData
{
    MIntArray    T,   o_T;
    MIntArray    eqc, o_eqc;
    MDoubleArray itv, o_itv;
    
    HalfedgeData(size_t o_nHE, size_t nHE):
    T(nHE, 0),     o_T(o_nHE, 0),    
    eqc(nHE, 0),  o_eqc(o_nHE, 0), 
    itv(nHE, 0),   o_itv(o_nHE, 0)
    { }
};

inline void remap_HalfedgeData(HalfedgeData &d, size_t srcIdx, size_t tgtIdx)
{
    d.T[tgtIdx]   = d.o_T[srcIdx];
    d.itv[tgtIdx] = d.o_itv[srcIdx];
    d.eqc[tgtIdx] = d.o_eqc[srcIdx];
}


struct VertexData
{
    MIntArray pole, o_pole;
    
    VertexData(size_t o_nV, size_t nV): pole(nV, 2), o_pole(o_nV, 2) { }
};

inline void remap_VertexData(VertexData &d, size_t srcIdx, size_t tgtIdx)
{
    d.pole[tgtIdx] = d.o_pole[srcIdx];
}


void get_vertex_blindData(MFnMesh &meshFn, VertexData &v)
{
    get_mesh_int_blindData(meshFn, MFn::kMeshVertComponent, POLE_BLINDDATA_ID, "pole", v.o_pole);
}

void get_halfedge_blindData(MFnMesh &meshFn, HalfedgeData &he)
{
    get_mesh_int_blindData(meshFn, MFn::kMeshFaceVertComponent, T_BLINDDATA_ID, "T", he.o_T);
    get_mesh_int_blindData(meshFn, MFn::kMeshFaceVertComponent, EQC_BLINDDATA_ID, "eqc", he.o_eqc);
    get_mesh_double_blindData(meshFn, MFn::kMeshFaceVertComponent, ITV_BLINDDATA_ID, "itv", he.o_itv);
}

void set_vertex_blindData(MFnMesh &meshFn, VertexData &v)
{
    set_mesh_int_blindData(meshFn, MFn::kMeshFaceVertComponent, POLE_BLINDDATA_ID, "pole", v.pole);
}

void set_halfedge_blindData(MFnMesh &meshFn, HalfedgeData &he)
{
    set_mesh_int_blindData(meshFn, MFn::kMeshFaceVertComponent, T_BLINDDATA_ID, "T", he.T);
    set_mesh_int_blindData(meshFn, MFn::kMeshFaceVertComponent, EQC_BLINDDATA_ID, "eqc", he.eqc);
    set_mesh_double_blindData(meshFn, MFn::kMeshFaceVertComponent, ITV_BLINDDATA_ID, "itv", he.itv);
}


#include "kvc.h"

struct TCCNodeData: public KVCStruct
{
    KVCMutableRealArray<uint64_t>      nFV;
    KVCMutableIndexArray<uint64_t>     F;

    KVCMutableRealArray<unsigned char> pole;
    KVCMutableRealArray<bool>          T;
    KVCMutableRealArray<double>        itv;
    KVCMutableRealArray<uint64_t>      eqc;
    
    KVCMutableIndexArray<uint64_t>     bad;           // faces with bad topology: nondyadic T-mesh or non-analysis-suitable
    KVCMutableIndexArray<uint64_t>     open_TF;       // faces with unassigned T-joints
    KVCMutableIndexArray<uint64_t>     open_eqc;      // faces in the same eqc as openTF
    KVCMutableIndexArray<uint64_t>     mismatched;    // itv values not matched throughout eqc
    KVCMutableIndexArray<uint64_t>     inconsistent;  // inconsistent eqc: there is no way to assign matching itv values.
    
    TCCNodeData()
    {
        std::string  k[] = { "pole", "T", "itv", "eqc", "nFV", "F", "bad", "open_TF", "open_eqc", "mismatched", "inconsistent" };
        KVCObject   *v[] = { &pole,  &T,  &itv,  &eqc , &nFV,  &F , &bad,  &open_TF,  &open_eqc,  &mismatched,  &inconsistent };
        KVC(k, v, sizeof(v)/sizeof(KVCObject *));
    }
};


// YOUR STUFF UNTIL HERE -----------------------------------------------/

void write_to_TCCNodeData(MIntArray &nFV, MIntArray &F, VertexData &v, HalfedgeData &he, TCCNodeData &nd)
{
    // vertex data 
    nd.pole.setDims(1,v.pole.length());
    for (size_t k=0; k<v.pole.length(); k++) nd.pole[k] = v.pole[k];

    // halfedge data
    nd.T.setDims(1,he.T.length());
    for (size_t k=0; k<he.T.length(); k++) nd.T[k] = he.T[k];
    
    nd.itv.setDims(1,he.itv.length());
    for (size_t k=0; k<he.itv.length(); k++) nd.itv[k] = he.itv[k];
    
    nd.eqc.setDims(1,he.eqc.length());
    for (size_t k=0; k<he.eqc.length(); k++) nd.eqc[k] = he.eqc[k];

    
    nd.nFV.setDims(1,nFV.length());
    for (size_t k=0; k<nFV.length(); k++) nd.nFV[k] = nFV[k];

    nd.F.setDims(1,F.length());
    for (size_t k=0; k<F.length(); k++) nd.F[k] = F[k];
}

void read_from_TCCNodeData(TCCNodeData &nd, VertexData &v, HalfedgeData &he)
{
    for (size_t k=0; k<v.pole.length(); k++) v.pole[k] = nd.pole[k];
    for (size_t k=0; k<he.T.length(); k++) he.T[k] = nd.T[k];
    for (size_t k=0; k<he.itv.length(); k++) he.itv[k] = nd.itv[k];
    for (size_t k=0; k<he.eqc.length(); k++) he.eqc[k] = nd.eqc[k];
}

MStatus update_tags(MIntArray &nFV, MIntArray &F, VertexData &v, HalfedgeData &he)
{
    // export into TCCNodeData, send to MATLAB, read back results and store
    TCCNodeData nd;
    
    write_to_TCCNodeData(nFV, F, v, he, nd);
    
    char localhost[] = "localhost";
    int remote_socket = socketConnect(3003, localhost);
    cout<<"socket: "<<remote_socket<<endl;
    
    if (!send_KVCObject_to_Socket(&nd, remote_socket))
    {
        cout<<"sending failed!"<<endl;
        close(remote_socket);
        return MS::kUnknownParameter;
    }
    
    if (!init_KVCObject_from_Socket(&nd, remote_socket)) 
    {
        cout<<"receiving failed!"<<endl;
        close(remote_socket);
        return MS::kUnknownParameter;
    }
    close(remote_socket);
    
    read_from_TCCNodeData(nd, v, he);
    
    return MS::kSuccess;
}



size_t compute_nHE(MIntArray &polyOrder, MIntArray &o_nFV, MIntArray &delta_nFV)
{
    size_t kD = 0, nHE = 0;
    for (size_t k=0; k<polyOrder.length(); k++)
    {
        if (polyOrder[k]>=0) 
            nHE += o_nFV[polyOrder[k]];
        else 
            nHE += delta_nFV[kD++];
    }
    
    return nHE;
}


MStatus updateTCCDataFty::remapMeshData(MFnMesh &meshFn)
{
    MIntArray o_nFV, o_F;
    meshFn.getVertices(o_nFV, o_F);
    size_t o_nF  = o_nFV.length();
    size_t o_nHE = o_F.length();
    size_t o_nV  = meshFn.numVertices();
    
    MIntArray o_F2H(o_nF); // maps face index to halfedge index
    {
        size_t kS=0;
        for (size_t k=0; k<o_nF; k++)
        {
            o_F2H[k] = kS;
            kS+=o_nFV[k];
        }
    }
    
    size_t nF  = fPolyOrder.length();
    size_t nHE = compute_nHE(fPolyOrder, o_nFV, fDelta_nFV);
    size_t nV  = fnV;
    
    HalfedgeData he(o_nHE, nHE);
    VertexData   v(o_nV, nV);
    
    MFloatPointArray V(nV);
    MIntArray nFV(nF);
    MIntArray F(nHE);
    
    get_vertex_blindData(meshFn, v);
    get_halfedge_blindData(meshFn, he);
    
    // add old/new faces in the order given by fPolyOrder and remap their halfedge data (if exists)
    
    size_t delta_kF = 0, delta_kHE=0;
    size_t kF = 0, kHE = 0;
    for (size_t k=0; k<fPolyOrder.length(); k++)
    {
        if (fPolyOrder[k]>=0) 
        {
            size_t cF = fPolyOrder[k], cnFV = o_nFV[cF];
            
            size_t o_cF2H = o_F2H[cF];
            size_t cShift = fCShift[cF];
            for (size_t kFV=0; kFV<cnFV; kFV++) 
            {
                size_t o_kHE = o_cF2H + ((kFV + cShift) % cnFV);
                
                F[kHE] = fVtxRemap[o_F[o_kHE]];
                
                remap_HalfedgeData(he, o_kHE, kHE);
                
                kHE++;
            }

            nFV[kF] = cnFV; kF++;
        } 
        else // this is a new face
        { 
            size_t cnFV = fDelta_nFV[delta_kF]; delta_kF++;
            for (size_t kFV=0; kFV<cnFV; kFV++) 
            {
                F[kHE] = fDelta_F[delta_kHE];
                
                kHE++; delta_kHE++;
            }
            nFV[kF] = cnFV; kF++;
        }
    }    
    
    // remap vertex data
    for (size_t k=0; k<o_nV; k++)
    {
        if (fVtxRemap[k]>=0)
        {
            remap_VertexData(v, k, fVtxRemap[k]);
        }
    }
    
    MStatus stat = meshFn.createInPlace(nV, nF, V, nFV, F);
    if (stat != MS::kSuccess)
    {
        std::cerr<<"createInPlace failed"<<endl;
        std::cerr<<stat.errorString()<<endl;
    }
    
    
    update_tags(nFV, F, v, he);
    
    set_vertex_blindData(meshFn, v);
    set_halfedge_blindData(meshFn, he);
    
    return stat;
}

MStatus updateTCCDataFty::doIt()
//
//	Description:
//		Performs the actual updateTCCData operation on the given object and UVs
//
{
	MStatus status = MS::kSuccess;

    MFnMesh meshFn( fMesh );
    
    remapMeshData(meshFn);
	
    return status;
}
