import maya.OpenMaya as om;
import maya.cmds as cmds;

def selected_dagPath():
    selection = om.MSelectionList(); 
    om.MGlobal.getActiveSelectionList(selection);
    selIt = om.MItSelectionList(selection);
    if selIt.isDone():
        raise Exception("No object selected!");
    dagPath = om.MDagPath();
    selIt.getDagPath(dagPath);
    return dagPath;

def selected_dagPath_component():
    selection = om.MSelectionList(); 
    om.MGlobal.getActiveSelectionList(selection);
    selIt = om.MItSelectionList(selection);
    if selIt.isDone():
        raise Exception("No object selected!");
    dagPath = om.MDagPath(); component = om.MObject(); 
    selIt.getDagPath(dagPath, component);
    return (dagPath, component);

def selected_DGNode():
    selection = om.MSelectionList(); 
    om.MGlobal.getActiveSelectionList(selection);
    selIt = om.MItSelectionList(selection);
    if selIt.isDone():
        raise Exception("No object selected!");
    dagPath = om.MDagPath(); component = om.MObject(); 
    depNode = om.MObject();
    stat = selIt.getDependNode(depNode);
    DGNode = om.MFnDependencyNode( depNode );
    return DGNode;

def selected_edges():
    (dagPath, comp) = selected_dagPath_component();
    
    if comp.apiType() != om.MFn.kMeshEdgeComponent:
        raise Exception("No edge selection!");

    e = om.MIntArray();
    compFn = om.MFnSingleIndexedComponent(comp);
    compFn.getElements(e);
    return e;
    
def selected_vertices():
    (dagPath, comp) = selected_dagPath_component();
    
    if comp.apiType() != om.MFn.kMeshVertComponent:
        raise Exception("No vertex selection!");

    v = om.MIntArray();
    compFn = om.MFnSingleIndexedComponent(comp);
    compFn.getElements(v);
    return v;
    
def transfer_Selection():
    hilite = mel.eval("ls -hl");
    vlist = selected_vertices()
    cmds.select(cl=True);
    for v in vlist:
        cmds.select(hilite[1]+".vtx["+str(v)+"]", tgl=True);

def set_double_blinddata(meshFn, compType, bdId, attrName, valArr):
    if meshFn.isBlindDataTypeUsed(bdId):
        meshFn.clearBlindData(compType, bdId, attrName);
    else:
        meshFn.createBlindDataType(bdId, [attrName], [attrName], ["double"]);
    
    idx = om.MIntArray(len(valArr));
    val = om.MDoubleArray(len(valArr));
    for k in range(len(val)):
        idx[k] = k;
        val[k] = valArr[k];
    
    meshFn.setDoubleBlindData(idx, compType, bdId, attrName, val);

def edge_vertices(dagPath, eI):
    MSU = om.MScriptUtil(); a = MSU.asIntPtr();
    eIt = om.MItMeshEdge(dagPath);
    eIt.setIndex(eI, a); 
    v = om.MIntArray(2);
    v[0] = eIt.index(0);
    v[1] = eIt.index(1);
    return v;


def edge_faces(dagPath, eI):
    MSU = om.MScriptUtil(); a = MSU.asIntPtr();
    eIt = om.MItMeshEdge(dagPath);
    eIt.setIndex(eI, a); 
    f = om.MIntArray();  eIt.getConnectedFaces(f);
    return f;
    
def vertex_faces(dagPath, vI):
    MSU = om.MScriptUtil(); a = MSU.asIntPtr();
    vIt = om.MItMeshVertex(dagPath);
    vIt.setIndex(vI, a); 
    f = om.MIntArray(); vIt.getConnectedFaces(f);
    return f;
    
def face_edges(dagPath, fI):
    MSU = om.MScriptUtil(); a = MSU.asIntPtr();
    fIt = om.MItMeshPolygon(dagPath);
    fIt.setIndex(fI, a); 
    e = om.MIntArray(); fIt.getEdges(e);
    return e;

def face_faces(dagPath, fI):
    MSU = om.MScriptUtil(); a = MSU.asIntPtr();
    fIt = om.MItMeshPolygon(dagPath);
    fIt.setIndex(fI, a); 
    f = om.MIntArray(); fIt.getConnectedFaces(f);
    return f;

def face_vertices(dagPath, fI):
    MSU = om.MScriptUtil(); a = MSU.asIntPtr();
    fIt = om.MItMeshPolygon(dagPath);
    fIt.setIndex(fI, a); 
    v = om.MIntArray(); fIt.getVertices(v);
    return v;

def vertex_edges(dagPath, vI):
    MSU = om.MScriptUtil(); a = MSU.asIntPtr();
    vIt = om.MItMeshVertex(dagPath);
    vIt.setIndex(vI, a); 
    e = om.MIntArray(); vIt.getConnectedEdges(e);
    return e;

def face_edge_index(e, eI):
    for k in range(len(e)):
        if e[k] == eI:
            return k;
    return -1;
    

def getTopology(meshFn):
    nFV = om.MIntArray(); F = om.MIntArray();
    meshFn.getVertices(nFV, F);
    return (nFV, F);

# Halfedge data structure wrapper

def halfedge(dP, f, eI):
    e0 = face_edges(dP, f);
    k0 = face_edge_index(e0, eI);
    return (f, k0, e0, dP);

def halfedges(dP, eI):
    f = edge_faces(dP, eI);
    he = [];
    for cF in f: he.append( halfedge(dP, cF, eI) );
    return he;
    
def edge(he):
    return he[2][he[1]];

def face(he):
    return he[0];

def next(he, dir=1):
    kN = (he[1] + dir) % len(he[2]);
    return (he[0], kN, he[2], he[3]);

def prev(he, dir=1):
    kP = (he[1] - dir) % len(he[2]);
    return (he[0], kP, he[2], he[3]);

def twin(he):
    f = edge_faces(he[3], edge(he));
    if f[0]==he[0]:
        f = f[1];
    else:
        f = f[0];
    e = face_edges(he[3], f);
    k = face_edge_index(e, edge(he));
    return (f, k, e, he[3]);
    
def tip(he, dir=1):
    v = face_vertices(he[3], he[0]);
    if (dir==1): he = next(he);
    return v[ he[1] ];
# below: paranoid version: does not rely on the assumption  v0 -e0-> v1 -e1-> v2 ...
#
#    va = edge_vertices(he[3], edge(he));
#    vb = edge_vertices(he[3], edge(next(he, dir)));
#    if va[0] == vb[0] or va[0] == vb[1]:
#        v = va[0];
#    else:
#        v = va[1];
#    return v;
  
  
def get_halfedge_data(dagPath, data, he):
    meshFn = om.MFnMesh(dagPath);
    
    MSU = om.MScriptUtil(); idxPtr = MSU.asIntPtr();
    meshFn.getFaceVertexBlindDataIndex(he[0], tip(he,1), idxPtr);
    idx = om.MScriptUtil(idxPtr).asInt();
    
    if (len(data)<idx):
        raise Exception("data index out of range. Did you modify mesh topology and forget to update?");
    
#    print("idx: ",idx, "e", he[2][he[1]]);
#    print(data[idx]);

    return data[idx];    

def get_TCC_DGnode(dagPath):
    dagPath.extendToShape();
    meshShapeFn = om.MFnDependencyNode(dagPath.node());
            
    outMeshAttr = meshShapeFn.attribute( "outMesh" );
    outMeshPlug = om.MPlug(dagPath.node(), outMeshAttr );
            
    plugArray = om.MPlugArray();
    outMeshPlug.connectedTo(plugArray, True, True);

    for k in range(plugArray.length()):
        TCCnode = plugArray[k].node();
        TCCDepNode = om.MFnDependencyNode( TCCnode );
         
        if TCCDepNode.typeName() == "TCC":
            return TCCDepNode;
            
    Exception("No TCC node found");




