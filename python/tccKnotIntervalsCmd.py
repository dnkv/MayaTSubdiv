import maya.OpenMaya as om
import maya.OpenMayaMPx as ompx
import maya.cmds as cmds

import mayaTools as mt
from vec3 import *

import sys;

kPluginCmdName = "tccKnotIntervals"

#####################################################################
## COMMAND ##########################################################
#####################################################################

class tccKnotIntervals(ompx.MPxCommand):

    def __init__(self):
        ompx.MPxCommand.__init__(self)

    def doIt(self, args):

        if args.length()<1:
            raise RuntimeError("have to specify 'insert'/'remove'!");
            
        if args.asString(0)=="double":
            multiplier = 2.0;
        elif args.asString(0)=="half":
            multiplier = 0.5;
        else:
            raise RuntimeError("have to specify 'double'/'half'!");
    
        (dagPath, comp) = mt.selected_dagPath_component();
        dagPath.extendToShape();
        
        meshFn = om.MFnMesh(dagPath);
        nFV = om.MIntArray(); F = om.MIntArray(); 
        meshFn.getVertices(nFV, F);
        
        edgeIt = om.MItMeshEdge(dagPath, comp);
        
        
        outMeshAttr = meshFn.attribute( "outMesh" );
        outMeshPlug = om.MPlug(dagPath.node(), outMeshAttr );
        plugArray = om.MPlugArray(); outMeshPlug.connectedTo(plugArray, True, True);
                    
        foundTCC = False;
        for k in range(0, plugArray.length()):
            TCCnode = plugArray[k].node();
            TCCDepNode = om.MFnDependencyNode( TCCnode );
            if TCCDepNode.typeName() == "TCC":
                foundTCC = True; break;
                
        if not foundTCC: 
            raise Exception('did not find TCC node');

        eqc = cmds.getAttr(TCCDepNode.name()+".eqc");
        itv = cmds.getAttr(TCCDepNode.name()+".itv");
        
        he = mt.halfedges(dagPath, edgeIt.index()); he = he[0];
        
        targetEqc = mt.get_halfedge_data(dagPath, eqc, he);

        for k in range(len(eqc)): 
            if eqc[k]==targetEqc: itv[k] = itv[k] * multiplier;
            
        cmds.setAttr(TCCDepNode.name()+".itv", itv, type="doubleArray");
        
                    
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
        