#include "updateTCCDataNode.h"

// Function Sets
//
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnData.h>
#include <maya/MFnNumericData.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMesh.h>

// General Includes
//
#include <maya/MGlobal.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MIOStream.h>

// Macros
//
#define MCheckStatus(status,message)    \
    if( MStatus::kSuccess != status ) {    \
        cerr << message << "\n";        \
        return status;                    \
    }


// Unique Node TypeId
//
MTypeId     updateTCCDataNode::id( 0x34500 );

// Node attributes (in addition to inMesh and outMesh defined by polyModifierNode)
//
MObject     updateTCCDataNode::polyOrder;
MObject     updateTCCDataNode::cShift;
MObject     updateTCCDataNode::vtxRemap;
MObject     updateTCCDataNode::delta_F;
MObject     updateTCCDataNode::delta_nFV;
MObject     updateTCCDataNode::nV;

updateTCCDataNode::updateTCCDataNode()
{
    cout<<"Node Constructor"<<endl;
}

updateTCCDataNode::~updateTCCDataNode()
{}

MStatus updateTCCDataNode::compute( const MPlug& plug, MDataBlock& data )
//
//    Description:
//        This method computes the value of the given output plug based
//        on the values of the input attributes.
//
//    Arguments:
//        plug - the plug to compute
//        data - object that provides access to the attributes for this node
//
{
    MStatus status = MS::kSuccess;
 
    MDataHandle stateData = data.outputValue( state, &status );
    MCheckStatus( status, "ERROR getting state" );

    // Check for the HasNoEffect/PassThrough flag on the node.
    //
    // (stateData is an enumeration standard in all depend nodes - stored as short)
    // 
    // (0 = Normal)
    // (1 = HasNoEffect/PassThrough)
    // (2 = Blocking)
    // ...
    //
    if( stateData.asShort() == 1 )
    {
        MDataHandle inputData = data.inputValue( inMesh, &status );
        MCheckStatus(status,"ERROR getting inMesh");

        MDataHandle outputData = data.outputValue( outMesh, &status );
        MCheckStatus(status,"ERROR getting outMesh");

        // Simply redirect the inMesh to the outMesh for the PassThrough effect
        //
        outputData.set(inputData.asMesh());
    }
    else
    {
        // Check which output attribute we have been asked to 
        // compute. If this node doesn't know how to compute it, 
        // we must return MS::kUnknownParameter
        // 
        if (plug == outMesh)
        {
            MDataHandle inputData = data.inputValue( inMesh, &status );
            MCheckStatus(status,"ERROR getting inMesh");

            MDataHandle outputData = data.outputValue( outMesh, &status );
            MCheckStatus(status,"ERROR getting outMesh"); 

            MIntArray vR = MFnIntArrayData( data.inputValue( vtxRemap ).data() ).array(&status);
            MCheckStatus(status,"ERROR getting vtxRemap");
            
            MIntArray pO = MFnIntArrayData( data.inputValue( polyOrder ).data() ).array(&status);
            MCheckStatus(status,"ERROR getting polyOrder");

            MIntArray cS = MFnIntArrayData( data.inputValue( cShift ).data() ).array(&status);
            MCheckStatus(status,"ERROR getting cShift");
            
            MIntArray dnFV = MFnIntArrayData( data.inputValue( delta_nFV ).data() ).array(&status);
            MCheckStatus(status,"ERROR getting deltanFV");

            MIntArray dF = MFnIntArrayData( data.inputValue( delta_F ).data() ).array(&status);
            MCheckStatus(status,"ERROR getting deltaF");
            
            int nVtx = data.inputValue( nV ).asInt();
            MCheckStatus(status,"ERROR getting nV");
            
            // Copy the inMesh to the outMesh, and now you can
            // perform operations in-place on the outMesh
            //
            outputData.set(inputData.asMesh());
            MObject mesh = outputData.asMesh();
            
            fupdateTCCDataFactory.setMesh( mesh );
            fupdateTCCDataFactory.setVtxRemap( vR );
            fupdateTCCDataFactory.setPolyOrder( pO );
            fupdateTCCDataFactory.setCShift( cS );
            fupdateTCCDataFactory.setDelta_nFV( dnFV );
            fupdateTCCDataFactory.setDelta_F( dF );
            fupdateTCCDataFactory.setnV( nVtx );

            // Now, perform the updateTCCData
            //
            status = fupdateTCCDataFactory.doIt();

            // Mark the output mesh as clean
            //
            outputData.setClean();
        }
        else
        {
            status = MS::kUnknownParameter;
        }
    }

    return status;
}

void* updateTCCDataNode::creator()
//
//    Description:
//        this method exists to give Maya a way to create new objects
//      of this type. 
//
//    Return Value:
//        a new object of this type
//
{
    return new updateTCCDataNode();
}

MStatus updateTCCDataNode::initialize()
//
//    Description:
//        This method is called to create and initialize all of the attributes
//      and attribute dependencies for this node type.  This is only called 
//        once when the node type is registered with Maya.
//
//    Return Values:
//        MS::kSuccess
//        MS::kFailure
//        
{
    MStatus             status;

    MFnTypedAttribute   attrFn;
    MFnNumericAttribute attrNum;

    MIntArray defaultArray;
    MFnIntArrayData defaultArrayDataFn; 
    
    defaultArrayDataFn.create( defaultArray );
    vtxRemap = attrFn.create("vtxRemap", "vR", MFnData::kIntArray, defaultArrayDataFn.object());
    attrFn.setStorable(true);

    defaultArrayDataFn.create( defaultArray );
    polyOrder = attrFn.create("polyOrder", "pO", MFnData::kIntArray, defaultArrayDataFn.object());
    attrFn.setStorable(true);
    
    defaultArrayDataFn.create( defaultArray );
    cShift = attrFn.create("cShift", "cS", MFnData::kIntArray, defaultArrayDataFn.object());
    attrFn.setStorable(true);
    
    defaultArrayDataFn.create( defaultArray );
    delta_F = attrFn.create("delta_F", "dF", MFnData::kIntArray, defaultArrayDataFn.object());
    attrFn.setStorable(true);
    
    defaultArrayDataFn.create( defaultArray );
    delta_nFV = attrFn.create("delta_nFV", "dnFV", MFnData::kIntArray, defaultArrayDataFn.object());
    attrFn.setStorable(true);
    
    nV = attrNum.create("nV", "nV", MFnNumericData::kInt, 0.0);
    attrFn.setStorable(true);
    
    inMesh = attrFn.create("inMesh", "iM", MFnMeshData::kMesh);
    attrFn.setStorable(true);


    outMesh = attrFn.create("outMesh", "om", MFnMeshData::kMesh);
    attrFn.setStorable(false);
    attrFn.setWritable(false);


    
    
    status = addAttribute( vtxRemap );
    MCheckStatus(status,"ERROR addAttribute"); 

    status = addAttribute( polyOrder );
    MCheckStatus(status,"ERROR addAttribute"); 

    status = addAttribute( cShift );
    MCheckStatus(status,"ERROR addAttribute"); 
    
    status = addAttribute( delta_F );
    MCheckStatus(status,"ERROR addAttribute"); 

    status = addAttribute( delta_nFV );
    MCheckStatus(status,"ERROR addAttribute"); 

    status = addAttribute( nV );
    MCheckStatus(status,"ERROR addAttribute"); 

    status = addAttribute( inMesh );
    MCheckStatus(status,"ERROR addAttribute"); 

    status = addAttribute( outMesh);
    MCheckStatus(status,"ERROR addAttribute"); 



    
    status = attributeAffects( inMesh, outMesh );
    MCheckStatus(status,"ERROR attributeAffects"); 

    status = attributeAffects( vtxRemap, outMesh );
    MCheckStatus(status,"ERROR attributeAffects"); 

    status = attributeAffects( polyOrder, outMesh );
    MCheckStatus(status,"ERROR attributeAffects"); 

    status = attributeAffects( cShift, outMesh );
    MCheckStatus(status,"ERROR attributeAffects"); 
    
    status = attributeAffects( delta_F, outMesh );
    MCheckStatus(status,"ERROR attributeAffects"); 
    
    status = attributeAffects( delta_nFV, outMesh );
    MCheckStatus(status,"ERROR attributeAffects"); 
    
    status = attributeAffects( nV, outMesh );
    MCheckStatus(status,"ERROR attributeAffects"); 

    return MS::kSuccess;

}
