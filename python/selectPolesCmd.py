import maya.OpenMaya as om
import maya.OpenMayaMPx as ompx
import maya.cmds as cmds
import maya.mel as mel

import mayaTools as mt
from vec3 import *

import sys;

kPluginCmdName = "selectPoles"

#####################################################################
## COMMAND ##########################################################
#####################################################################

class selPoles(ompx.MPxCommand):

    def __init__(self):
        ompx.MPxCommand.__init__(self)

    def doIt(self, args):
        dagPath = mt.selected_dagPath();
        meshFn = om.MFnMesh(dagPath);

        mel.eval("SelectVertexMask");
        cmds.select(cl=True);

        tccNode = mt.get_TCC_DGnode(dagPath);
        T = om.MIntArray(meshFn.numVertices(), 0);
        if tccNode is not None:
            Tb = cmds.getAttr(tccNode.name()+".T");
            (nFV, F) = mt.getTopology(meshFn);
            for k in range(len(F)):
                T[F[k]] += Tb[k];

        polesArray = [];
        connectedEdges = om.MIntArray()
        vertexIter=om.MItMeshVertex(dagPath)
        while not vertexIter.isDone():
            vertexIter.getConnectedEdges(connectedEdges);
            valence = len(connectedEdges) + T[vertexIter.index()];
            if (not vertexIter.onBoundary()) and (valence!=4):
                polesArray.append(vertexIter.index());
            vertexIter.next()
            
        for k in polesArray:
            cmds.select(dagPath.fullPathName()+".vtx["+str(k)+"]", tgl=True);

        
                    
def cmdCreator():
    return ompx.asMPxPtr(selPoles())

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
        