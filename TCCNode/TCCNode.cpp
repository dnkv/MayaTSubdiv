#include <maya/MTime.h>
#include <maya/MFnMesh.h>
#include <maya/MPoint.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MDoubleArray.h>
#include <maya/MFloatArray.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnPlugin.h>
#include <maya/MStatus.h>

#include <maya/MPxNode.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MFnMeshData.h>

#include <maya/MIOStream.h>

#include "blinddata.h"
#include "hds.h"
#include "tcc_max.h"

#include "msockets.h"
#include "kvcMatlabVar.h"

#define McheckErr(stat,msg)            \
    if ( MS::kSuccess != stat ) {    \
        stat.perror(msg);           \
        return MS::kFailure;        \
    }


struct TCCData 
{
    MIntArray    &nFV, &F;
    MIntArray    &nFVc, &Fc; // cached versions

    MIntArray    &pole, &corner;
    MIntArray    &T;
    MIntArray    &eqc;
    MDoubleArray &itv;
    MIntArray    &err;
    
    TCCData(MIntArray &rnFV, MIntArray &rF, 
            MIntArray &rnFVc, MIntArray &rFc, 
            MIntArray &rP, MIntArray &rC, MIntArray &rT, MIntArray &rEqc, MDoubleArray &rItv, MIntArray &rErr)
    : nFV(rnFV), F(rF), nFVc(rnFVc), Fc(rFc), pole(rP), corner(rC), T(rT), eqc(rEqc), itv(rItv), err(rErr) { }
};

class TCC : public MPxNode
{
public:
    TCC() { };
    virtual         ~TCC() {};
    virtual MStatus compute(const MPlug& plug, MDataBlock& data);

    static  void*   creator();
    static  MStatus initialize();

    
    static MObject    aRes, aRefRes;
    static MObject    aLineThickness;
    static MObject    aInputMesh, aOutputMesh;

    static MObject    anFVc, aFc;
    static MObject    aPole, aCorner, aT, aEqc, aItv, aErr;
    
    static MTypeId    id;
    
protected:
    bool    validTopology(TCCData &tccData);
    MStatus setErrorColors(MFnMesh &meshFn, TCCData &tccData);
    MStatus createSubdividedMesh(int sdRes, int sdRefRes, MFnMesh &srcMesh, TCCData &tccData, MDataHandle outMeshHandle, float lineThickness, MStatus& stat);
    
private:
    void copy_HDS(HDS &a, HDS &b);
};

MObject TCC::aRes;
MObject TCC::aLineThickness;
MObject TCC::aRefRes;
MObject TCC::aInputMesh;

MObject TCC::anFVc;
MObject TCC::aFc;

MObject TCC::aPole;
MObject TCC::aCorner;
MObject TCC::aT;
MObject TCC::aEqc;
MObject TCC::aItv;
MObject TCC::aErr;

MObject TCC::aOutputMesh;



MTypeId TCC::id( 0x80004 );

void* TCC::creator()
{
    return new TCC;
}

MStatus TCC::initialize()
{
    MFnNumericAttribute numAttr;
    MFnTypedAttribute attrFn;

    MFnIntArrayData defaultIArrayDataFn; 
    MFnDoubleArrayData defaultDArrayDataFn; 
    
    MStatus stat;

    aRes = numAttr.create( "SubdivisionResolution", "res", MFnNumericData::kInt, 3 );
    aRefRes = numAttr.create( "SubdivisionRefinementResolution", "refres", MFnNumericData::kInt, 0 );
    aLineThickness = numAttr.create( "UVLineThickness", "th", MFnNumericData::kFloat, 0 );
    
    
    aInputMesh = attrFn.create("inputMesh", "inMesh", MFnMeshData::kMesh);
    
    defaultIArrayDataFn.create( );
    anFVc = attrFn.create("nFVcached", "nFVc", MFnData::kIntArray, defaultIArrayDataFn.object());
    attrFn.setConnectable(false);
    attrFn.setStorable(true);

    defaultIArrayDataFn.create( );
    aFc = attrFn.create("Fcached", "Fc", MFnData::kIntArray, defaultIArrayDataFn.object());
    attrFn.setConnectable(false);
    attrFn.setStorable(true);

    defaultIArrayDataFn.create( );
    aPole = attrFn.create("pole", "pole", MFnData::kIntArray, defaultIArrayDataFn.object());
    attrFn.setConnectable(false);
    attrFn.setStorable(true);
    
    defaultIArrayDataFn.create( );
    aCorner = attrFn.create("corner", "corner", MFnData::kIntArray, defaultIArrayDataFn.object());
    attrFn.setConnectable(false);
    attrFn.setStorable(true);
    
    defaultIArrayDataFn.create( );
    aT = attrFn.create("T", "T", MFnData::kIntArray, defaultIArrayDataFn.object());
    attrFn.setConnectable(false);
    attrFn.setStorable(true);
    
    defaultIArrayDataFn.create( );
    aEqc = attrFn.create("eqc", "eqc", MFnData::kIntArray, defaultIArrayDataFn.object());
    attrFn.setConnectable(false);
    attrFn.setStorable(true);

    defaultDArrayDataFn.create( );
    aItv = attrFn.create("itv", "itv", MFnData::kDoubleArray, defaultDArrayDataFn.object());
    attrFn.setConnectable(false);
    attrFn.setStorable(true);
    
    defaultIArrayDataFn.create( );
    
    aErr = attrFn.create("err", "err", MFnData::kIntArray, defaultIArrayDataFn.object());
    attrFn.setConnectable(false);
    attrFn.setStorable(true);
    
    aOutputMesh = attrFn.create( "outputMesh", "out", MFnData::kMesh);
    attrFn.setStorable(false);
    attrFn.setWritable(false);

    stat = addAttribute(aRes); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aRefRes); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aLineThickness); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aInputMesh); McheckErr(stat, "ERROR adding attribute\n");
    
    stat = addAttribute(anFVc); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aFc); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aPole); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aCorner); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aT); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aEqc); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aItv); McheckErr(stat, "ERROR adding attribute\n");
    stat = addAttribute(aErr); McheckErr(stat, "ERROR adding attribute\n");
    
    stat = addAttribute(aOutputMesh); McheckErr(stat, "ERROR adding attribute\n");

    stat = attributeAffects(aRes, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aRefRes, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aLineThickness, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aInputMesh, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    
    stat = attributeAffects(anFVc, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aFc, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aPole, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aCorner, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aT, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aEqc, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aItv, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");
    stat = attributeAffects(aErr, aOutputMesh); McheckErr(stat, "ERROR in attributeAffects\n");

    return MS::kSuccess;
}

void store_in_hds(HDS &hds, MFloatPointArray &points, MIntArray &nFV, MIntArray &F)
{
    size_t nV = points.length();
    hds.V.setDims(3, nV);
    for (size_t k=0; k<nV; k++) 
    {
        hds.V[3*k+0] = points[k](0);
        hds.V[3*k+1] = points[k](1);
        hds.V[3*k+2] = points[k](2);
    }
    
    hds.nFV.setDims(1,nFV.length());
    for (size_t k=0; k<nFV.length(); k++) hds.nFV[k] = nFV[k];
    
    hds.tip.setDims(1,F.length());
    for (size_t k=0; k<F.length(); k++) hds.tip[k] = F[k];
}                      

void load_from_hds(HDS &hds, MFloatPointArray &points, MIntArray &nFV, MIntArray &F)
{
    size_t nV   = hds.nV();
    size_t nF   = hds.nF();
    size_t nIHE = hds.nIHE();
    
    points.setLength(nV);
    for (size_t k=0; k<nV; k++) 
    {
        points[k](0) = hds.V[3*k+0];
        points[k](1) = hds.V[3*k+1];
        points[k](2) = hds.V[3*k+2];
    }
    
    nFV.setLength(nF);
    for (size_t k=0; k<nF; k++) nFV[k] = hds.nFV[k];
    
    F.setLength(nIHE);
    for (size_t k=0; k<nIHE; k++) F[k] = hds.tip[k];
}

void create_subdivided_face(int sdRes, int nV, double *uvs, double *itv, int Tidx, MFloatArray &uA, MFloatArray &vA, MIntArray &uvIdx)
{
    HDS hds;
    
    hds.V.setDims(2, nV); memcpy(&hds.V.v[0], uvs, 2*nV*sizeof(double));
    hds.nFV.setDims(1, 1); hds.nFV[0] = nV;

    hds.tip.setDims(1, nV); 
    for (size_t k=0; k<nV; k++) 
    {
        hds.tip[k] = k;
    }
    finalize_HDS(hds);
    
    size_t nHE = hds.nHE(), nIHE = hds.nIHE();⟵
    hds.T.setDims(1, nHE);
    hds.itv.setDims(1, nHE);

    memset(&hds.T.v[0], 0, nHE*sizeof(bool)); if (nV==5) hds.T[Tidx] = 1;
    memcpy(&hds.itv.v[0], itv, nV*sizeof(double));

    // border halfedge tags
    for (size_t k=nIHE; k<nHE; k++) 
    {
        hds.itv[k] = hds.itv[hds.twin[k]];
    }
    
    TCC_MAX::linear_subdivide(hds, sdRes);
    
    int sd_nV = hds.nV();
    int sd_nIHE = hds.nIHE();

    uA.setLength(sd_nV);
    vA.setLength(sd_nV);
    for (int k=0; k<sd_nV; k++)
    {
        uA[k] = hds.V[2*k+0];
        vA[k] = hds.V[2*k+1];
    }
    
    uvIdx.setLength(sd_nIHE);
    for (int k=0; k<sd_nIHE; k++)
    {
        uvIdx[k]=hds.tip[k];
    }
}


// uScale and vScale switched because first edge pointing in v direction
void copy_patch_uvs(int &kV, int &kHE, MFloatArray &u_src, MFloatArray &v_src, MIntArray &uvIdx_src, MFloatArray &u_dst, MFloatArray &v_dst, float vScale, float uScale, MFloatArray &sc_u_dst, MFloatArray &sc_v_dst, MIntArray &uvIdx_dst, float lineThickness)
{
    int nHE = uvIdx_src.length();
    for (int k=0; k<nHE; k++)
    {
        uvIdx_dst[kHE] = kV + uvIdx_src[k]; // offset by beginning of this block of UVs
        
        kHE++;
    }
    
//    float offset_u = (uScale - 1) * 0.5*lineThickness;
//    float offset_v = (vScale - 1) * 0.5*lineThickness;
    
    lineThickness = 0.5 * lineThickness; 
    
    float offset_u = lineThickness * (uScale - 1.0)/(uScale - 2*lineThickness);
    float offset_v = lineThickness * (vScale - 1.0)/(vScale - 2*lineThickness);

    int nV = u_src.length();
    for (int k=0; k<nV; k++)
    {
        u_dst[kV] = u_src[k] * (1.0 - 2.0*offset_u) + offset_u;
        v_dst[kV] = v_src[k] * (1.0 - 2.0*offset_v) + offset_v;
        
        sc_u_dst[kV] = uScale * u_src[k];
        sc_v_dst[kV] = vScale * v_src[k];
        
        kV++;
    }
}

void createUVset(TCCData &tccData, int sdRes, MFloatArray &uArray, MFloatArray &vArray, MFloatArray &sc_uArray, MFloatArray &sc_vArray, MIntArray &uvIdx, float lineThickness)
{
    MFloatArray u4, v4, u5T0, v5T0, u5T1, v5T1, u5T2, v5T2, u5T3, v5T3, u5T4, v5T4;
    MIntArray id4, id5T0, id5T1, id5T2, id5T3, id5T4;
    
    double hds4_uvs[]   = { 0,1, 1,1, 1,0, 0,0 },       hds4_itv[] = {1,1,1,1};       create_subdivided_face(sdRes, 4, hds4_uvs, hds4_itv, -1, u4, v4, id4);
    double hds5T0_uvs[] = { 0,.5, 0,1, 1,1, 1,0, 0,0 }, hds5T0_itv[] = {.5,.5,1,1,1}; create_subdivided_face(sdRes, 5, hds5T0_uvs, hds5T0_itv, 0, u5T0, v5T0, id5T0);
    double hds5T1_uvs[] = { 0,0, 0,.5, 0,1, 1,1, 1,0 }, hds5T1_itv[] = {1,.5,.5,1,1}; create_subdivided_face(sdRes, 5, hds5T1_uvs, hds5T1_itv, 1, u5T1, v5T1, id5T1);
    double hds5T2_uvs[] = { 1,0, 0,0, 0,.5, 0,1, 1,1 }, hds5T2_itv[] = {1,1,.5,.5,1}; create_subdivided_face(sdRes, 5, hds5T2_uvs, hds5T2_itv, 2, u5T2, v5T2, id5T2);
    double hds5T3_uvs[] = { 1,1, 1,0, 0,0, 0,.5, 0,1 }, hds5T3_itv[] = {1,1,1,.5,.5}; create_subdivided_face(sdRes, 5, hds5T3_uvs, hds5T3_itv, 3, u5T3, v5T3, id5T3);
    double hds5T4_uvs[] = { 0,1, 1,1, 1,0, 0,0, 0,.5 }, hds5T4_itv[] = {.5,1,1,1,.5}; create_subdivided_face(sdRes, 5, hds5T4_uvs, hds5T4_itv, 4, u5T4, v5T4, id5T4);
    
    
    int nV4 = u4.length(), nHE4 = id4.length();
    int nV5 = u5T0.length(), nHE5 = id5T0.length();
    
    int nF = tccData.nFV.length();

    int nHE = 0, nV = 0;
    for (int kF = 0; kF<nF; kF++)
    {
        switch (tccData.nFV[kF])
        {
            case 4: nV+=nV4; nHE+=nHE4; break;
            case 5: nV+=nV5; nHE+=nHE5; break;
        }
    }
    
    uArray.setLength(nV); sc_uArray.setLength(nV);
    vArray.setLength(nV); sc_vArray.setLength(nV);
    uvIdx.setLength(nHE);
    
    int kHE = 0, sd_kV = 0, sd_kHE = 0;
    for (int kF = 0; kF<nF; kF++)
    {
        int faceConfig = tccData.nFV[kF];

        if (faceConfig==5)
        {
            for (int kT=0; kT<5; kT++)
            {
                if (tccData.T[kHE+kT]) break;
                faceConfig++;
            }
        }

        // uScale and vScale are chosen to not use T-edges, but match the U/V directions of the patch (convention: first edge gets vScale)
        switch (faceConfig)
        {
            case 4: copy_patch_uvs(sd_kV, sd_kHE, u4, v4, id4, uArray, vArray, tccData.itv[kHE], tccData.itv[kHE+1], sc_uArray, sc_vArray, uvIdx, lineThickness); break;
            case 5: copy_patch_uvs(sd_kV, sd_kHE, u5T0, v5T0, id5T0, uArray, vArray, tccData.itv[kHE+3], tccData.itv[kHE+4], sc_uArray, sc_vArray, uvIdx, lineThickness); break;
            case 6: copy_patch_uvs(sd_kV, sd_kHE, u5T1, v5T1, id5T1, uArray, vArray, tccData.itv[kHE+4], tccData.itv[kHE], sc_uArray, sc_vArray, uvIdx, lineThickness); break;
            case 7: copy_patch_uvs(sd_kV, sd_kHE, u5T2, v5T2, id5T2, uArray, vArray, tccData.itv[kHE], tccData.itv[kHE+1], sc_uArray, sc_vArray, uvIdx, lineThickness); break;
            case 8: copy_patch_uvs(sd_kV, sd_kHE, u5T3, v5T3, id5T3, uArray, vArray, tccData.itv[kHE+1], tccData.itv[kHE+2], sc_uArray, sc_vArray, uvIdx, lineThickness); break;
            case 9: copy_patch_uvs(sd_kV, sd_kHE, u5T4, v5T4, id5T4, uArray, vArray, tccData.itv[kHE+2], tccData.itv[kHE+3], sc_uArray, sc_vArray, uvIdx, lineThickness); break;
        }
        
        
        kHE += tccData.nFV[kF];
    }
}


void TCC::copy_HDS(HDS &a, HDS &b) // copies from a to b
{
    b.nFV = a.nFV;
    b.F2H = a.F2H;
    
    b.H2F = a.H2F;
    b.next = a.next; 
    b.prev = a.prev; 
    b.twin = a.twin;
    
    b.tip = a.tip; 
    
    b.T = a.T;
    b.pole = a.pole;
    b.corner = a.corner;
    
    b.itv = a.itv;
    
    b.V2H = a.V2H;
    
    b.V = a.V;
}

MStatus TCC::createSubdividedMesh(int sdRes, int sdRefRes, MFnMesh &srcMesh, TCCData &tccData, MDataHandle outMeshHandle, float lineThickness, MStatus& stat)
{   
    HDS hds;
    bool shouldCreateUVs = true;
    
    size_t nV = srcMesh.numVertices();
    size_t nF = srcMesh.numPolygons();
    size_t nIHE = tccData.F.length();
    
    bool consistentSizes= (tccData.pole.length()==nV) && (tccData.T.length()==nIHE) && (tccData.itv.length()==nIHE) & (tccData.corner.length()==nV);
    
    if ((nV==0)||(nF==0)||(!consistentSizes)) return MS::kFailure;

    MFloatArray uArray, vArray, sc_uArray, sc_vArray;
    MIntArray uvIdx;
    if (shouldCreateUVs)
    {
        createUVset(tccData, sdRes, uArray, vArray, sc_uArray, sc_vArray, uvIdx, lineThickness);
    }
    
    
    MFloatPointArray points;
    srcMesh.getPoints(points);
    
    store_in_hds(hds, points, tccData.nFV, tccData.F);     // convert to HDS

    finalize_HDS(hds);
    size_t nHE = hds.nHE();

    hds.T.setDims(1, nHE);
    hds.itv.setDims(1, nHE);
    hds.corner.setDims(1, nV);
    
    // interior halfedge tags
    for (size_t k=0; k<nV; k++) 
    {
        hds.corner[k] = tccData.corner[k];
    }

    // interior halfedge tags
    for (size_t k=0; k<nIHE; k++) 
    {
        hds.T[k] = tccData.T[k];
        hds.itv[k] = tccData.itv[k];
    }
    
    // border halfedge tags
    for (size_t k=nIHE; k<nHE; k++) 
    {
        hds.T[k] = false;
        hds.itv[k] = hds.itv[hds.twin[k]];
    }
        
    TCC_MAX::subdivide(hds, sdRes);
    
    if (sdRefRes>0)
    {
        HDS hds2;
        copy_HDS(hds, hds2);
        TCC_MAX::subdivide(hds2, sdRefRes);
        memcpy(&hds.V[0], &hds2.V[0], hds.V.size() * sizeof(double));
    }
    

    
    MObject outMeshObj = outMeshHandle.asMesh();
    MFnMesh outMeshFn(outMeshObj);
    
    // if no topology change necessary, just update points!
    if ( (outMeshFn.numFaceVertices() == hds.nIHE()) && (outMeshFn.numPolygons() == hds.nF()) )
    {
        size_t nV   = hds.nV();
        points.setLength(nV);
        for (size_t k=0; k<nV; k++) 
        {
            points[k](0) = hds.V[3*k+0];
            points[k](1) = hds.V[3*k+1];
            points[k](2) = hds.V[3*k+2];
        }
        stat = outMeshFn.setPoints(points); McheckErr(stat, "ERROR creating outputData");
        
        if (shouldCreateUVs)
        {
            MString uvSet = "UnitPatchUVs";
            MString sc_uvSet = "ScaledPatchUVs";
            stat = outMeshFn.setUVs(uArray, vArray, &uvSet); McheckErr(stat, "ERROR setting UVs");
            stat = outMeshFn.setUVs(sc_uArray, sc_vArray, &sc_uvSet); McheckErr(stat, "ERROR setting UVs");
        }
        
        return MS::kSuccess;
    }

    
    // Have to update connectivity and geometry

    load_from_hds(hds, points, tccData.nFV, tccData.F);

    nV = points.length();
    nF = tccData.nFV.length();
    
    MFnMeshData dataCreator;
    MObject newOutputData = dataCreator.create(&stat); McheckErr(stat, "ERROR creating outputData");
    MFnMesh newOutMeshFn;
    
    MObject newMesh;
    
    newMesh = newOutMeshFn.create(nV, nF, points, tccData.nFV, tccData.F, newOutputData, &stat); McheckErr(stat, "ERROR in MFnMesh.create\n");

    if (shouldCreateUVs)
    {
        MString uvSet = "UnitPatchUVs";
        MString sc_uvSet = "ScaledPatchUVs";
        
        uvSet = newOutMeshFn.createUVSetDataMeshWithName(uvSet, &stat); McheckErr(stat, "ERROR creating UVset");
        stat = newOutMeshFn.clearUVs(&uvSet);
        stat = newOutMeshFn.setUVs(uArray, vArray, &uvSet); McheckErr(stat, "ERROR setting UVs");
        stat = newOutMeshFn.assignUVs(tccData.nFV, uvIdx, &uvSet); McheckErr(stat, "ERROR assigning UVs");
        
        sc_uvSet = newOutMeshFn.createUVSetDataMeshWithName(sc_uvSet, &stat); McheckErr(stat, "ERROR creating UVset");
        stat = newOutMeshFn.clearUVs(&sc_uvSet);
        stat = newOutMeshFn.setUVs(sc_uArray, sc_vArray, &sc_uvSet); McheckErr(stat, "ERROR setting UVs");
        stat = newOutMeshFn.assignUVs(tccData.nFV, uvIdx, &sc_uvSet); McheckErr(stat, "ERROR assigning UVs");
    }
    
    
    if (stat == MS::kSuccess)
    {
        outMeshHandle.set(newOutputData);
    }
    
    
    return MS::kSuccess;
}

bool TCC::validTopology(TCCData &tccData)
{
    if (tccData.nFV.length() != tccData.nFVc.length()) return false;
    if (tccData.F.length() != tccData.Fc.length()) return false;


    bool hasErr=false;
    for (size_t k=0; k<tccData.err.length(); k++) hasErr = hasErr | (tccData.err[k]>0);
    if (hasErr) return false;
    
    size_t nHE = tccData.F.length();
    for (size_t k=0; k<nHE; k++) 
    {
        if (tccData.F[k] != tccData.Fc[k]) return false;
    }

    return true; 
}


MStatus TCC::setErrorColors(MFnMesh &meshFn, TCCData &tccData)
{
    MStatus stat;
    if ((tccData.err.length()>0) && (tccData.err.length()==tccData.nFV.length()))
    {
        for (size_t k=0; k<tccData.err.length(); k++)
        {
            int err = tccData.err[k];
            
            if (err>0)
            {
                MColor col(1.0, 1.0, 1.0);
                switch(err)
                {
                    case 1: col = MColor(1.0, 0.0, 0.0); break; // bad topology / inconsistent eqc
                    case 2: col = MColor(1.0, 0.0, 1.0); break; // mismatched itv
                    case 3: col = MColor(0.0, 0.0, 1.0); break; // unassigned itvs
                    case 4: col = MColor(0.5, 0.5, 1.0); break; // unassigned T-joint
                    default: break;
                }
                stat = meshFn.setFaceColor(col, k);
                if (stat!=MS::kSuccess)
                {
                    char bla[] = "bla";
                    stat.perror(bla);
                }
            }
        }
    }
    
    return stat;
}


MStatus TCC::compute(const MPlug& plug, MDataBlock& data)

{
    MStatus stat;

    if (plug == aOutputMesh) {
        /* Get time */
        int subdivRes = data.inputValue(aRes, &stat).asInt();
        int subdivRefRes = data.inputValue(aRefRes, &stat).asInt();
        float lineThickness = data.inputValue(aLineThickness, &stat).asFloat();

        MDataHandle inMeshHandle = data.inputValue( aInputMesh, &stat ); McheckErr(stat,"ERROR getting attribute");
        MObject inMeshObj = inMeshHandle.asMesh();
        MFnMesh inMeshFn(inMeshObj);
        
        MIntArray nFV, F;
        inMeshFn.getVertices(nFV, F);
        
        MIntArray nFVc = MFnIntArrayData( data.inputValue( anFVc ).data() ).array(&stat); McheckErr(stat,"ERROR getting attr");
        MIntArray Fc   = MFnIntArrayData( data.inputValue( aFc ).data() ).array(&stat); McheckErr(stat,"ERROR getting attr");
        
        MIntArray pole = MFnIntArrayData( data.inputValue( aPole ).data() ).array(&stat); McheckErr(stat,"ERROR getting attr");
        MIntArray corner = MFnIntArrayData( data.inputValue( aCorner ).data() ).array(&stat); McheckErr(stat,"ERROR getting attr");
        MIntArray T    = MFnIntArrayData( data.inputValue( aT ).data() ).array(&stat); McheckErr(stat,"ERROR getting attr");
        MIntArray eqc  = MFnIntArrayData( data.inputValue( aEqc ).data() ).array(&stat); McheckErr(stat,"ERROR getting attr");
        MDoubleArray itv  = MFnDoubleArrayData( data.inputValue( aItv ).data() ).array(&stat); McheckErr(stat,"ERROR getting attr");
        MIntArray err  = MFnIntArrayData( data.inputValue( aErr ).data() ).array(&stat); McheckErr(stat,"ERROR getting attr");

        TCCData tccData(nFV, F, nFVc, Fc, pole, corner, T, eqc, itv, err);

        /* Get output object */
        MDataHandle outMeshHandle = data.outputValue(aOutputMesh, &stat); McheckErr(stat, "ERROR getting attribute\n");

        if (validTopology(tccData))
        {
            stat = createSubdividedMesh(subdivRes, subdivRefRes, inMeshFn, tccData, outMeshHandle, lineThickness, stat);
        } 
        else 
        {
            outMeshHandle.setMObject(inMeshObj);
            
            MFnMesh outMeshFn(outMeshHandle.asMesh());
            
            stat = setErrorColors(outMeshFn, tccData);
        }
        
        data.setClean( plug );
    } else
        return MS::kUnknownParameter;

    return stat;
}

MStatus initializePlugin(MObject obj)
{
    MStatus   status;
    MFnPlugin plugin(obj, PLUGIN_COMPANY, "3.0", "Any");

    status = plugin.registerNode("TCC", TCC::id, TCC::creator, TCC::initialize);
    if (!status) {
        status.perror("registerNode");
        return status;
    }

    return status;
}

MStatus uninitializePlugin(MObject obj)
{
    MStatus      status;
    MFnPlugin plugin(obj);

    status = plugin.deregisterNode(TCC::id);
    if (!status) {
        status.perror("deregisterNode");
        return status;
    }

    return status;
}
