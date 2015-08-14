import maya.OpenMaya as om
import maya.OpenMayaMPx as ompx
import maya.cmds as cmds
import maya.mel as mel

import mayaTools as mt
from vec3 import *

import sys;

kPluginCmdName = "selectValence2Vertices"

#####################################################################
## COMMAND ##########################################################
#####################################################################

class tccKnotIntervals(ompx.MPxCommand):

    def __init__(self):
        ompx.MPxCommand.__init__(self)

    def doIt(self, args):
        dagPath = mt.selected_dagPath();
        meshFn = om.MFnMesh(dagPath);
        mel.eval("SelectVertexMask");
        cmds.select(cl=True);
        for k in range(meshFn.numVertices()):
            ve = mt.vertex_edges(dagPath, k);
            vf = mt.vertex_faces(dagPath, k);
            if len(ve)==2 and len(vf)>1:
                cmds.select(dagPath.fullPathName()+".vtx["+str(k)+"]", tgl=True);
        
                    
def cmdCreator():
    return ompx.asMPxPtr(tccKnotIntervals())

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
        