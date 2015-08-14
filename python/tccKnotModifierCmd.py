import maya.OpenMaya as om
import maya.OpenMayaMPx as ompx
import maya.cmds as cmds

import mayaTools as mt
import tccKnotModifier as km
from vec3 import *

import sys;

kPluginCmdName = "tccKnotModifier"

#####################################################################
## COMMAND ##########################################################
#####################################################################

class tccKnotModifier(ompx.MPxCommand):

    def __init__(self):
        ompx.MPxCommand.__init__(self)

    def doIt(self, args):
        if args.length()<1:
            raise RuntimeError("have to specify 'insert'/'remove'!");
            
        if args.asString(0)=="insert":
            data = km.knot_modifier( km.KnotModifier.Insert );
        elif args.asString(0)=="remove":
            data = km.knot_modifier( km.KnotModifier.Remove );
        else:
            raise RuntimeError("have to specify 'insert'/'remove'!");
        
        dagPath = mt.selected_dagPath();
        dagPath.extendToShape();
        
        nameStr = dagPath.fullPathName();
        
        for vp in data:
            (v, p) = vp;
            cmds.move(p.x, p.y, p.z, nameStr+".vtx["+str(v)+"]", absolute=True, objectSpace=True);
            
def cmdCreator():
    return ompx.asMPxPtr(tccKnotModifier())

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
