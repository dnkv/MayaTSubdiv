import maya.OpenMaya as om;
import maya.OpenMayaMPx as ompx;
import maya.cmds as cmds;

import mayaTools as mt;

import sys;

kPluginCmdName = "tccToggleCorner"

#####################################################################
## COMMAND ##########################################################
#####################################################################

class tccToggleCorner(ompx.MPxCommand):

    def __init__(self):
        ompx.MPxCommand.__init__(self)

    def doIt(self, args):
        v = mt.selected_vertices();
        dagPath = mt.selected_dagPath();
        TCC = mt.get_TCC_DGnode(dagPath);
        corner = cmds.getAttr(TCC.name()+".corner");
        
        # check if all elements in selection have same value
        val = corner[v[0]];
        allEqual = True;
        for cv in v:
            if corner[cv] != val:
                allEqual = False; break;
        if allEqual:
            val = (val+1) % 3;
        else:
            val = 0;
                
        for cv in v: corner[cv] = val;
        
        cornerType = ["off", "on", "auto"];
        print "All vertices set to corner="+cornerType[val];
    
        cmds.setAttr(TCC.name()+".corner", corner, type="Int32Array");
            
def cmdCreator():
    return ompx.asMPxPtr(tccToggleCorner())

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
    