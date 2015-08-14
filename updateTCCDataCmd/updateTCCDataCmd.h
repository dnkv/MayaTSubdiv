#ifndef _updateTCCDataCmd
#define _updateTCCDataCmd

#include <maya/MPxCommand.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MDGModifier.h>
#include <maya/MDagPath.h>
#include <maya/MObject.h>
#include <maya/MSelectionList.h>

class MArgList;

struct TCCData 
{
    MIntArray    &nFV, &o_nFV;
    MIntArray    &F,   &o_F;
    
    MIntArray    &pole, &o_pole;
    MIntArray    &corner, &o_corner;

    MIntArray    &T,    &o_T;
    MIntArray    &eqc,  &o_eqc;
    MDoubleArray &itv,  &o_itv;

    MIntArray    &err;
    
    MIntArray    &selHE;
    
    TCCData (MIntArray &anFV,    MIntArray &ao_nFV,
             MIntArray &aF,      MIntArray &ao_F,
             MIntArray &apole,   MIntArray &ao_pole, 
             MIntArray &acorner, MIntArray &ao_corner, 
             MIntArray &aT,      MIntArray &ao_T, 
             MIntArray &aeqc,    MIntArray &ao_eqc, 
             MDoubleArray &aitv, MDoubleArray &ao_itv, 
             MIntArray &aselHE,
             MIntArray &aerr)
    
    :nFV(anFV), o_nFV(ao_nFV), F(aF), o_F(ao_F), pole(apole), o_pole(ao_pole), corner(acorner), o_corner(ao_corner), 
    T(aT), o_T(ao_T), eqc(aeqc), o_eqc(ao_eqc), itv(aitv), o_itv(ao_itv), selHE(aselHE), err(aerr) { }
};



class updateTCCData : public MPxCommand
{

public:
    ////////////////////
    // Public Methods //
    ////////////////////

                updateTCCData();
    virtual    ~updateTCCData();

    static      void* creator();

    bool        isUndoable() const;

    MStatus     doIt( const MArgList& );
    MStatus     redoIt();
    MStatus     undoIt();
    

    /////////////////////////////
    // polyModifierCmd Methods //
    /////////////////////////////

private:
    /////////////////////
    // Private Methods //
    /////////////////////

    MStatus     getTCCNode( MSelectionList &selList );
    MStatus     update();
    bool        hasTopologyChanged();
    bool        hasUnassignedItvs();

    MStatus     backup_TCCData();
    MStatus     restore_TCCData();

    MStatus     compute_remap(TCCData &tccData, MIntArray &vR, MIntArray &pO, MIntArray &cS);
    MStatus     remapTCCData(TCCData &tccData, MIntArray &vR, MIntArray &polyOrder, MIntArray &circShift);

    //////////////////
    // Private Data //
    //////////////////
    bool          fDoDgModifier;

    size_t        fnV;
    MIntArray     fnFVc, fFc;
    MIntArray     fPole;
    MIntArray     fCorner;
    MIntArray     fT;
    MIntArray     fEqc;
    MDoubleArray  fItv;
    MIntArray     fErr;
    
    MDGModifier   fDgModifier;
    MDagPath      fSrcDagPath;
    MObject       fTCCnode;
    MObject       fComponent;
};

#endif
