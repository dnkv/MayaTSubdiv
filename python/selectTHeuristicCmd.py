import maya.OpenMaya as om
import maya.OpenMayaMPx as ompx
import maya.cmds as cmds
import maya.mel as mel

import mayaTools as mt

import sys;
import math;

kPluginCmdName = "selectTHeuristic"

#####################################################################
## COMMAND ##########################################################
#####################################################################

class selTHeuristic(ompx.MPxCommand):

    def __init__(self):
        ompx.MPxCommand.__init__(self)

    def doIt(self, args):
        dagPath = mt.selected_dagPath();
        meshFn = om.MFnMesh(dagPath);

        mel.eval("SelectVertexFaceMask");
        cmds.select(cl=True);

        vp1 = om.MPoint(); vp2 = om.MPoint(); vp3 = om.MPoint();

        (nFV, F) = mt.getTopology(meshFn);
        for kF in range(len(nFV)):
            if nFV[kF]==5:
                V = mt.face_vertices(dagPath, kF);
                maxAngle = 0.0;
                for kV in range(5):
                    meshFn.getPoint(V[(kV-1)%5], vp1)
                    meshFn.getPoint(V[kV],vp2);
                    meshFn.getPoint(V[(kV+1)%5], vp3);
                    a=vp1 - vp2; 
                    if a*a>0: a=a/a.length();
                    b=vp3 - vp2; 
                    if b*b>0: b=b/b.length();
                    angle = math.acos(a*b);
                    if angle>maxAngle:
                        T=V[kV]; maxAngle = angle;
                cmds.select(dagPath.fullPathName()+".vtxFace["+str(T)+"]["+str(kF)+"]", tgl=True);
                    
def cmdCreator():
    return ompx.asMPxPtr(selTHeuristic())

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
        