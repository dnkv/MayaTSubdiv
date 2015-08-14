import maya.OpenMaya as om
import maya.OpenMayaMPx as ompx
import maya.cmds as cmds
import mayaTools as mt

import sys;

kPluginCmdName = "assignUnitUVsCmd"

#####################################################################
## COMMAND ##########################################################
#####################################################################
#
# This command assigns "unit" UVs to a once subdivided Catmull-Clark 
# mesh, i.e. it assigns a quad corner [0..0.5]x[0..0.5] to each 
# quad face. The result is a UV map that contains all original 
# mesh edges either on U=0 or V=0.
#
# (C) 2013 Denis Kovacs / NYU

class assignUnitUVsClass(ompx.MPxCommand):

    def __init__(self):
        ompx.MPxCommand.__init__(self)

    def doIt(self, args):
        dagPath = mt.selected_dagPath();
        meshFn = om.MFnMesh(dagPath);

        (nFV, F) = mt.getTopology(meshFn);

        mU = om.MFloatArray(len(F),0.5);
        mV = om.MFloatArray(len(F),0.5);
        mI = om.MIntArray(len(F),0);

        TU = [0.5, 0.5, 0,   0];
        TV = [0.5,   0, 0, 0.5];

        idx=0;
        vp1 = om.MPoint(); vp2 = om.MPoint(); vp3 = om.MPoint();

        for kF in range(len(nFV)):
                V = mt.face_vertices(dagPath, kF);
                FC = 0; maxV = 0;
                for kV in range(4):
                    if V[kV]>maxV:
                        FC=kV; maxV=V[kV];
                        
                for kFV in range(4):
                    mU[idx] = TU[(kFV-FC)%4];
                    mV[idx] = TV[(kFV-FC)%4];
                    mI[idx] = idx;
                    idx = idx + 1;
                continue;
                
        meshFn.clearUVs();
        meshFn.setUVs(mU, mV);
        meshFn.assignUVs(nFV, mI);
                    
def cmdCreator():
    return ompx.asMPxPtr(assignUnitUVsClass())

def initializePlugin(mobject):
    mplugin = ompx.MFnPlugin(mobject, "Autodesk", "1.0", "Any")
    try:
        mplugin.registerCommand(kPluginCmdName, cmdCreator)
    except:
        sys.stderr.write( "Failed to register command: %s\n" % kPluginCmdName)
        raise

def uninitializePlugin(mobject):
    mplugin = ompx.MFnPlugin(mobject)
    try:
        mplugin.deregisterCommand(kPluginCmdName)
    except:
        sys.stderr.write("Failed to unregister command: %s\n" % kPluginCmdName)
        raise
        