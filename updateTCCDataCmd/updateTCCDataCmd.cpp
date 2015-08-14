#include "updateTCCDataCmd.h"

// Function Sets

#include <maya/MFnPlugin.h>
//
#include <maya/MFnDependencyNode.h>
#include <maya/MFnMesh.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFloatPointArray.h>
#include <maya/MDagPath.h>
#include <maya/MFnDoubleIndexedComponent.h>

// Iterators
//
#include <maya/MItSelectionList.h>
#include <maya/MItMeshPolygon.h>

// General Includes
//
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>

#include <maya/MIOStream.h>

#include <vector>
#include <map>
#include "blinddata.h"
#include "TCCNodeData.h"
#include "kvcMatlabVar.h"

#include "updateTCC.h"

using namespace std;

#define MCheckStatus(status,message)    \
if( MS::kSuccess != status ) {        \
cerr << message << "\n";        \
return status;                    \
}


updateTCCData::updateTCCData()
{
}

updateTCCData::~updateTCCData()
{}

void* updateTCCData::creator()
{
    return new updateTCCData();
}

bool updateTCCData::isUndoable() const
{
    return true;
}

MStatus updateTCCData::getTCCNode( MSelectionList &selList )
{
    MStatus status; 
    MItSelectionList selListIter( selList );
    selListIter.setFilter( MFn::kMesh );
    
    for( ; !selListIter.isDone(); selListIter.next() )
    {
        MDagPath dagPath;
        MObject component;
        selListIter.getDagPath( dagPath, component );
        status = dagPath.extendToShape();
    
        if (status == MS::kSuccess)
        {
            MFnDependencyNode meshShapeFn(dagPath.node());
            
            MObject outMeshAttr = meshShapeFn.attribute( "outMesh" );
            MPlug outMeshPlug(dagPath.node(), outMeshAttr );
//            worldMeshPlug = worldMeshPlug.elementByLogicalIndex( 0 );
            cout<<"inMesh-connected dagnode name: "<<endl<<" > "<<outMeshPlug.info().asChar()<<endl;
            
            MPlugArray plugArray;
            outMeshPlug.connectedTo(plugArray, true, true);
            cout<<"connected attrs: "<<plugArray.length()<<endl;
            
            for (size_t k=0; k<plugArray.length(); k++)
            {
                MObject TCCnode(plugArray[k].node());
                MFnDependencyNode TCCDepNode( TCCnode );
                
                cout<<"Dependency node type: "<<TCCDepNode.typeName().asChar()<<endl;
                if (TCCDepNode.typeName() == "TCC")
                {
                    fSrcDagPath = dagPath;
                    fTCCnode = TCCnode;
                    fComponent = component;
                    return MS::kSuccess;
                }
            }
        }
    }    
    
    return MS::kFailure;
}


MStatus updateTCCData::doIt( const MArgList& )
{
    MStatus status;
    
    MSelectionList selList;
    MGlobal::getActiveSelectionList( selList );
    
    status = getTCCNode(selList); if (status!=MS::kSuccess) return status;
    
    status = backup_TCCData(); fDoDgModifier = false;
    
    fDgModifier.commandToExecute("resetVtxRemap "+fSrcDagPath.fullPathName());
    
    status = redoIt();
    
    return status;
}

MStatus updateTCCData::redoIt()
{
    MStatus status;
    
    if( hasTopologyChanged() || hasUnassignedItvs())
    {
        status = update();
        
        fDgModifier.doIt(); fDoDgModifier = true;
        
        if (status==MS::kSuccess)
        {
            setResult( "updateTCCData command succeeded!" );
        } 
        else
        {
            setResult( "updateTCCData command failed!" );
        }
    }
    else
    {
        displayError( "updateTCCData skipped: No topology change detected" );
        status = MS::kSuccess;
    }
    
    return status;
}

MStatus updateTCCData::undoIt()
{
    MStatus status = MS::kSuccess;
    
    status = restore_TCCData();
    
    if( status == MS::kSuccess )
    {
        if (fDoDgModifier)
        {
            fDgModifier.undoIt(); fDoDgModifier = false;
        }
        setResult( "updateTCCData undo succeeded!" );
    }
    else
    {
        setResult( "updateTCCData undo failed!" );
    }
    
    return status;
}

bool updateTCCData::hasTopologyChanged()
{

    return true;
}

bool updateTCCData::hasUnassignedItvs()
{

    return true;
}

MStatus getMIntArray(MObject &node, MString attrString, MIntArray &arr)
{
    MStatus stat; 
    MFnDependencyNode nodeFn(node);
    
    MObject attr = nodeFn.attribute( attrString, &stat );  MCheckStatus(stat,"ERROR node.attribute()"); 
    MPlug plug( node, attr );
    MDataHandle handle; plug.getValue( handle );
    arr.copy(MFnIntArrayData(handle.data()).array());
    plug.destructHandle(handle);
    
    return MS::kSuccess;
}

MStatus getMDoubleArray(MObject &node, MString attrString, MDoubleArray &arr)
{
    MStatus stat; 
    MFnDependencyNode nodeFn(node);
    
    MObject attr = nodeFn.attribute( attrString, &stat );  MCheckStatus(stat,"ERROR node.attribute()"); 
    MPlug plug( node, attr );
    MDataHandle handle; plug.getValue( handle );
    arr.copy(MFnDoubleArrayData(handle.data()).array());
    plug.destructHandle(handle);
    
    return MS::kSuccess;
}


MStatus updateTCCData::backup_TCCData()
{
    MStatus stat;

    stat = getMIntArray(fTCCnode, "nFVc", fnFVc); MCheckStatus(stat,"ERROR getMIntArray"); 
    stat = getMIntArray(fTCCnode, "Fc", fFc); MCheckStatus(stat,"ERROR getMIntArray"); 
    stat = getMIntArray(fTCCnode, "pole", fPole); MCheckStatus(stat,"ERROR getMIntArray"); 
    stat = getMIntArray(fTCCnode, "corner", fCorner); MCheckStatus(stat,"ERROR getMIntArray");
    stat = getMIntArray(fTCCnode, "T", fT); MCheckStatus(stat,"ERROR getMIntArray"); 
    stat = getMIntArray(fTCCnode, "eqc", fEqc); MCheckStatus(stat,"ERROR getMIntArray"); 
    stat = getMDoubleArray(fTCCnode, "itv", fItv); MCheckStatus(stat,"ERROR getMIntArray"); 
    stat = getMIntArray(fTCCnode, "err", fErr); MCheckStatus(stat,"ERROR getMIntArray");
    
 
    /*
    print_MIntArray(fnFVc);
    print_MIntArray(fFc);
    print_MIntArray(fPole);
    print_MIntArray(fT);
    print_MIntArray(fEqc);
    print_MDoubleArray(fItv);
    print_MIntArray(fErr);
    */
    
    return MS::kSuccess;
}


MStatus setMIntArray(MObject &node, MString attrString, MIntArray &arr)
{
    MStatus stat; 
    MFnDependencyNode nodeFn(node);
    
    MObject attr = nodeFn.attribute( attrString, &stat );  MCheckStatus(stat,"ERROR node.attribute()"); 
    MPlug plug( node, attr );
    MDataHandle handle; plug.getValue( handle );
    MFnIntArrayData(handle.data()).array().copy(arr);
    plug.destructHandle(handle);
    
    return MS::kSuccess;
}

MStatus setMDoubleArray(MObject &node, MString attrString, MDoubleArray &arr)
{
    MStatus stat; 
    MFnDependencyNode nodeFn(node);
    
    MObject attr = nodeFn.attribute( attrString, &stat );  MCheckStatus(stat,"ERROR node.attribute()"); 
    MPlug plug( node, attr );
    MDataHandle handle; plug.getValue( handle );
    MFnDoubleArrayData(handle.data()).array().copy(arr);
    plug.destructHandle(handle);
    
    return MS::kSuccess;
}


MStatus updateTCCData::restore_TCCData()
{
    MStatus stat;
    
    stat = setMIntArray(fTCCnode, "nFVc", fnFVc); MCheckStatus(stat,"ERROR setMIntArray"); 
    stat = setMIntArray(fTCCnode, "Fc", fFc); MCheckStatus(stat,"ERROR setMIntArray"); 
    stat = setMIntArray(fTCCnode, "pole", fPole); MCheckStatus(stat,"ERROR setMIntArray"); 
    stat = setMIntArray(fTCCnode, "corner", fCorner); MCheckStatus(stat,"ERROR setMIntArray"); 
    stat = setMIntArray(fTCCnode, "T", fT); MCheckStatus(stat,"ERROR setMIntArray"); 
    stat = setMIntArray(fTCCnode, "eqc", fEqc); MCheckStatus(stat,"ERROR setMIntArray"); 
    stat = setMDoubleArray(fTCCnode, "itv", fItv); MCheckStatus(stat,"ERROR setMIntArray"); 
    stat = setMIntArray(fTCCnode, "err", fErr); MCheckStatus(stat,"ERROR setMIntArray"); 
    
    return MS::kSuccess;
}

#define VTXREMAP_BLINDDATA_ID 9005
#define NO_DST_VTXID -1

MStatus matlab_update_TCCData(TCCData &tccData, MFloatPointArray &V);


MStatus compute_halfedge_indices(MIntArray &nFV, MIntArray &F, MIntArray &selV, MIntArray &selF, MIntArray &selHE)
{
    MIntArray F2H(nFV.length(), 0);
    
    size_t cumsum = 0;
    for (size_t k=0; k<nFV.length(); k++)
    {
        F2H[k] = cumsum; 
        cumsum += nFV[k];
    }
    
    selHE.setLength(selF.length());
    for (size_t k=0; k<selF.length(); k++)
    {
        size_t cV = selV[k];
        size_t cF = selF[k];
        size_t cnFV = nFV[cF];
        size_t cF2H = F2H[cF];
        for (size_t kFV=0; kFV<cnFV; kFV++)
        {
            if (F[cF2H+kFV]==cV)
            {
                selHE[k] = cF2H + kFV; break;
            }
        }
    }
    
    return MS::kSuccess;
}


MStatus updateTCCData::update()
{
    MStatus stat;    
    
    MIntArray nFV, F;
    MIntArray polyOrder, circShift;
    MFloatPointArray V;

    MFnMesh srcOutMesh(fSrcDagPath);
    srcOutMesh.getVertices(nFV, F);
    srcOutMesh.getPoints(V);
    fnV = srcOutMesh.numVertices();
    size_t nHE = F.length();
    
    MIntArray selHE;
    if (fComponent.apiType() == MFn::kMeshVtxFaceComponent)
    {
        MIntArray selF, selV;
        MFnDoubleIndexedComponent compFn(fComponent);
        compFn.getElements(selV, selF);
        compute_halfedge_indices(nFV, F, selV, selF, selHE);
    }
    
    MIntArray vR(fnV, NO_DST_VTXID);
    get_mesh_int_blindData(srcOutMesh, MFn::kMeshVertComponent, VTXREMAP_BLINDDATA_ID, "vtxRemap", vR);
    
    MIntArray newCorner(fnV, 2);
    MIntArray newPole(fnV, 2);
    MIntArray newT(nHE, 0);
    MIntArray newEqc(nHE, 0);
    MDoubleArray newItv(nHE, 0);
    MIntArray newErr;
    
    TCCData tccData(nFV,     fnFVc,
                    F,       fFc,
                    newPole, fPole,
                    newCorner, fCorner, 
                    newT,    fT, 
                    newEqc,  fEqc, 
                    newItv,  fItv, 
                    selHE, newErr);

    compute_remap(tccData, vR, polyOrder, circShift);
    
    remapTCCData(tccData, vR, polyOrder, circShift);
    
    stat = matlab_update_TCCData(tccData, V);
    
    if (stat == MS::kSuccess)
    {
        stat = setMIntArray(fTCCnode, "nFVc", nFV); MCheckStatus(stat,"ERROR setMIntArray"); 
        stat = setMIntArray(fTCCnode, "Fc",   F); MCheckStatus(stat,"ERROR setMIntArray"); 
        stat = setMIntArray(fTCCnode, "pole", newPole); MCheckStatus(stat,"ERROR setMIntArray"); 
        stat = setMIntArray(fTCCnode, "T",    newT); MCheckStatus(stat,"ERROR setMIntArray"); 
        stat = setMIntArray(fTCCnode, "eqc",  newEqc); MCheckStatus(stat,"ERROR setMIntArray"); 
        stat = setMDoubleArray(fTCCnode, "itv",  newItv); MCheckStatus(stat,"ERROR setMIntArray"); 
        stat = setMIntArray(fTCCnode, "err",  newErr); MCheckStatus(stat,"ERROR setMIntArray"); 
        stat = setMIntArray(fTCCnode, "corner",  newCorner); MCheckStatus(stat,"ERROR setMIntArray"); 
    }
    return stat;
}



MStatus updateTCCData::compute_remap(TCCData &tccData, MIntArray &vR, MIntArray &pO, MIntArray &cS)
{
    pair<multimap<int, int>::iterator,multimap<int, int>::iterator> fRange;
    multimap<int, int> dst_V2F; // maps smallest face vertex index to all its faces
    
    MIntArray &src_nFV(tccData.nFV), &dst_nFV(tccData.o_nFV);
    MIntArray &src_F(tccData.F),     &dst_F(tccData.o_F);
    
    // get source and destination F and nFV arrays
    size_t src_nHE = src_F.length(), src_nF = src_nFV.length();
    size_t dst_nHE = dst_F.length(), dst_nF = dst_nFV.length();
    
    pO.setLength(src_nF);
    cS.setLength(src_nF);
    
    MIntArray dst_F2H(dst_nF, 0);
    MIntArray dst_cS(dst_nF, 0);
    
    // build lookup map
    size_t kHE = 0;
    for (size_t kF=0; kF<dst_nF; kF++)
    {
        dst_F2H[kF] = kHE;
        size_t cnFV = dst_nFV[kF];
        size_t minVI = dst_nHE;     // big enough...
        size_t ccS = 0;
        for (size_t kFV=0; kFV<cnFV; kFV++)
        {
            if (dst_F[kHE]<minVI) 
            { 
                minVI = dst_F[kHE];
                ccS = kFV;
            }
            kHE++;
        }
        dst_V2F.insert(pair<int,int>(minVI, kF));
        dst_cS[kF] = ccS;
    }
    
    size_t src_F2H = 0;
    
    // look up new faces in old mesh
    for (size_t kF=0; kF<src_nF; kF++)
    {
        size_t cnFV = src_nFV[kF];
        
        // find minimum remapped vtx index, remap vertices
        size_t minVI = src_nHE;     // start with smth big enough...
        size_t src_ccS = 0;
        bool newFace = false;
        for (size_t kFV=0; kFV<cnFV; kFV++)
        {
            int vRI = vR[src_F[src_F2H+kFV]];
            
            if (vRI == NO_DST_VTXID) { newFace = true; break; }
            if (vRI < minVI) {
                minVI = vRI;
                src_ccS = kFV;
            }
        }
        

        // no new vertices? ok, then let's see if we can find this face in dst
        if (!newFace)
        {
            bool faceMismatch = true;
            
            fRange = dst_V2F.equal_range(minVI);
            for (multimap<int, int>::iterator it=fRange.first; it!=fRange.second; ++it)
            {
                int dst_rF = (*it).second; // remapped face candidate
                if (dst_nFV[dst_rF] != cnFV) continue; // face vertex count mismatch -> not the right face
                
                faceMismatch = false;
                
                size_t ccS = dst_cS[dst_rF] + cnFV - src_ccS; // additional cnFV makes ccS always positive.
                
                for (size_t kFV=0; kFV<cnFV; kFV++)
                {
                    // index in dst polygon, starting with the same vertex offset (relative to smallest index) as src polygon
                    size_t o_kHE = dst_F2H[dst_rF] + ((kFV + ccS) % cnFV);
                    
                    int src_vI = vR[ src_F[ src_F2H + kFV ] ];
                    int dst_vI = dst_F[ o_kHE ];
                    
                    if ( src_vI != dst_vI ) { faceMismatch = true; break; }
                }
                if (!faceMismatch) {
                    pO[kF] = dst_rF;
                    cS[kF] = ccS;
                    break;
                }
            }
            
            if (faceMismatch) newFace = true;
        }
        
        // if new face, add it to delta_F and delta_nFV
        if (newFace)
        {
            pO[kF] = NO_DST_VTXID;
        }
        
        src_F2H += cnFV;
    }
    
    return MS::kSuccess;
}


void remap_HalfedgeData(TCCData &tccData, size_t o_kHE, size_t kHE)
{
    if (tccData.o_T.length()>o_kHE)   tccData.T[kHE] = tccData.o_T[o_kHE];
    if (tccData.o_eqc.length()>o_kHE) tccData.eqc[kHE] = tccData.o_eqc[o_kHE];
    if (tccData.o_itv.length()>o_kHE) tccData.itv[kHE] = tccData.o_itv[o_kHE];
}

void remap_VertexData(TCCData &tccData, size_t o_kV, size_t kV)
{
    if (tccData.o_pole.length()>o_kV) tccData.pole[kV] = tccData.o_pole[o_kV];
    if (tccData.o_corner.length()>o_kV) tccData.corner[kV] = tccData.o_corner[o_kV];
}

MStatus updateTCCData::remapTCCData(TCCData &tccData, MIntArray &vR, MIntArray &polyOrder, MIntArray &circShift)
{
    MIntArray &o_nFV=tccData.o_nFV, &nFV=tccData.nFV;
    size_t o_nF  = o_nFV.length();
    
    MIntArray o_F2H(o_nF); // maps face index to halfedge index
    {
        size_t kS=0;
        for (size_t k=0; k<o_nF; k++)
        {
            o_F2H[k] = kS;
            kS+=o_nFV[k];
        }
    }
    
    size_t nF  = nFV.length();
    size_t nV  = fnV;
    
    // add old/new faces in the order given by fPolyOrder and remap their halfedge data (if exists)
    
    size_t kHE = 0;
    for (size_t k=0; k<nF; k++)
    {
        if (polyOrder[k]>=0) 
        {
            size_t cF = polyOrder[k], cnFV = o_nFV[cF];
            
            size_t o_cF2H = o_F2H[cF];
            size_t cShift = circShift[k];  // cShift is stored for new faces
            for (size_t kFV=0; kFV<cnFV; kFV++) 
            {
                size_t o_kHE = o_cF2H + ((kFV + cShift) % cnFV);
                
                remap_HalfedgeData(tccData, o_kHE, kHE);
                
                kHE++;
            }
        } 
        else // this is a new face
        { 
            kHE += nFV[k];
        }
    }    
    
    
    // remap vertex data
    for (size_t k=0; k<nV; k++)
    {
        if (vR[k]>=0)
        {
            remap_VertexData(tccData, vR[k], k);
        }
    }
    
    return MS::kSuccess;
}



void write_to_TCCNodeData(TCCData &tcc, MFloatPointArray &V, TCCNodeData &nd)
{
    size_t nV = V.length();
    nd.V.setDims(3, nV);
    for (size_t k=0; k<nV; k++) 
    {
        nd.V[3*k+0] = V[k](0);
        nd.V[3*k+1] = V[k](1);
        nd.V[3*k+2] = V[k](2);
    }
    
    // vertex data 
    nd.pole.setDims(1,tcc.pole.length());
    for (size_t k=0; k<tcc.pole.length(); k++) nd.pole[k] = tcc.pole[k]==1;

    nd.corner.setDims(1,tcc.corner.length());
    for (size_t k=0; k<tcc.corner.length(); k++) nd.corner[k] = tcc.corner[k];
    
    // halfedge data
    nd.T.setDims(1,tcc.T.length());
    for (size_t k=0; k<tcc.T.length(); k++) nd.T[k] = tcc.T[k];
    
    nd.itv.setDims(1,tcc.itv.length());
    for (size_t k=0; k<tcc.itv.length(); k++) nd.itv[k] = tcc.itv[k];
    
    nd.eqc.setDims(1,tcc.eqc.length());
    for (size_t k=0; k<tcc.eqc.length(); k++) nd.eqc[k] = tcc.eqc[k];

    // err is not needed!
    
    nd.nFV.setDims(1,tcc.nFV.length());
    for (size_t k=0; k<tcc.nFV.length(); k++) nd.nFV[k] = tcc.nFV[k];
    
    nd.tip.setDims(1,tcc.F.length());
    for (size_t k=0; k<tcc.F.length(); k++) nd.tip[k] = tcc.F[k];

    nd.selHE.setDims(1,tcc.selHE.length());
    for (size_t k=0; k<tcc.selHE.length(); k++) nd.selHE[k] = tcc.selHE[k];
}

void read_from_TCCNodeData(TCCNodeData &nd, TCCData &tcc)
{
    tcc.pole.setLength(length(nd.pole));
    for (size_t k=0; k<tcc.pole.length(); k++) tcc.pole[k] = nd.pole[k];
    tcc.corner.setLength(length(nd.corner));
    for (size_t k=0; k<tcc.corner.length(); k++) tcc.corner[k] = nd.corner[k];
    tcc.T.setLength(length(nd.T));
    for (size_t k=0; k<tcc.T.length(); k++)    tcc.T[k]    = nd.T[k];
    tcc.itv.setLength(length(nd.itv));
    for (size_t k=0; k<tcc.itv.length(); k++)  tcc.itv[k]  = nd.itv[k];
    tcc.eqc.setLength(length(nd.eqc));
    for (size_t k=0; k<tcc.eqc.length(); k++)  tcc.eqc[k]  = nd.eqc[k];
    tcc.err.setLength(length(nd.err));
    for (size_t k=0; k<tcc.err.length(); k++)  tcc.err[k]  = nd.err[k];
}

MStatus matlab_update_TCCData_MATLAB(TCCData &tccData, MFloatPointArray &V)
{
    // export into TCCNodeData, send to MATLAB, read back results and store
    TCCNodeData nd;
    
    write_to_TCCNodeData(tccData, V, nd);
    
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
    
    read_from_TCCNodeData(nd, tccData);
    
    return MS::kSuccess;
}


MStatus matlab_update_TCCData(TCCData &tccData, MFloatPointArray &V)
{
    // export into TCCNodeData, send to MATLAB, read back results and store
    TCCNodeData nd;
    
    write_to_TCCNodeData(tccData, V, nd);
    
    updateTCC(nd);    
    
    read_from_TCCNodeData(nd, tccData);
    
    return MS::kSuccess;
}


MStatus initializePlugin( MObject obj )
//
//	Description:
//		this method is called when the plug-in is loaded into Maya.  It 
//		registers all of the services that this plug-in provides with 
//		Maya.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{ 
	MStatus   status;
	MFnPlugin plugin( obj, PLUGIN_COMPANY, "4.0", "Any");
    
	status = plugin.registerCommand( "updateTCCData", updateTCCData::creator );
	if (!status) {
		status.perror("registerCommand");
		return status;
	}
    
	return status;
}

MStatus uninitializePlugin( MObject obj )
//
//	Description:
//		this method is called when the plug-in is unloaded from Maya. It 
//		deregisters all of the services that it was providing.
//
//	Arguments:
//		obj - a handle to the plug-in object (use MFnPlugin to access it)
//
{
	MStatus   status;
	MFnPlugin plugin( obj );
    
	status = plugin.deregisterCommand( "updateTCCData" );
	if (!status) {
		status.perror("deregisterCommand");
		return status;
	}
    
	return status;
}
