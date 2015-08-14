from vec3 import *;
import mayaTools as mt; reload(mt); 

import maya.OpenMaya as om;
import maya.cmds as cmds;
from collections import defaultdict;


class KnotModifier:
    Insert, Remove = range(2);

def get_CV_line(dagPath, v, eK1, eK2):
    meshFn = om.MFnMesh(dagPath);
    
    e = mt.vertex_edges(dagPath, v);
    
    if e[0]==eK1 or e[0]==eK2:
        e = [e[1], e[3]];
    else:
        e = [e[0], e[2]];

    tccNode = mt.get_TCC_DGnode(dagPath);
    T = cmds.getAttr(tccNode.name()+".T");
      
    he = []; hv = [];
    for eK in e:
        f = mt.edge_faces(dagPath, eK);
    
        he0 = mt.halfedge(dagPath, f[0], eK1);
        if he0[1]==-1:
            he0 = mt.halfedge(dagPath, f[0], eK2);
        eN = mt.edge(mt.next(he0));
        dir = 0;
        if eN == e[0] or eN == e[1]:
            dir = 1;
        else:
            dir = -1;
            
        he0 = mt.next(he0, dir); 
        T0 = mt.get_halfedge_data(dagPath, T, he0) if dir ==1 else mt.get_halfedge_data(dagPath, T, mt.prev(he0));
        he1 = mt.next(he0, dir) if T0 else mt.next(mt.twin(mt.next(he0, dir)),dir);
        T1 = mt.get_halfedge_data(dagPath, T, he1) if dir ==1 else mt.get_halfedge_data(dagPath, T, mt.prev(he1));
        he2 = mt.next(he1, dir) if T1 else mt.next(mt.twin(mt.next(he1, dir)),dir);
        
        he.append(he0);
        he.append(he1);
        he.append(he2);
        
        hv.append(mt.tip(he0, dir));
        hv.append(mt.tip(he1, dir));
        hv.append(mt.tip(he2, dir));
        
        
    he = [ he[2], he[1], he[0], he[3], he[4], he[5] ];
    hv = [ hv[2], hv[1], hv[0], v, hv[3], hv[4], hv[5] ];
    
    return (he, hv);
        
        
def compute_knot_removal(p, i):
#  --i[0]--|--i[1]--|--i[2]--|--i[3]--|--i[4]--|--i[5]--
#         p[0]     p[1]    (p[2])    p[3]     p[4]

    # relocate neighbors 
    p1 = ( p[1] - p[0]*(i[3] / (i[0]+i[1]+i[2]+i[3]) ) ) / ( (i[0]+i[1]+i[2]) / (i[0]+i[1]+i[2]+i[3]) );
    p3 = ( p[3] - p[4]*(i[2] / (i[2]+i[3]+i[4]+i[5]) ) ) / ( (i[3]+i[4]+i[5]) / (i[2]+i[3]+i[4]+i[5]) );

    # this is the position we would get if we did another knot insertion again.
    L = (p[1] * (i[0] + i[1] + i[2] + i[3]) - p[0] * (i[3])) / (i[0] + i[1] + i[2]);
    R = (p[3] * (i[5] + i[4] + i[3] + i[2]) - p[4] * (i[2])) / (i[5] + i[4] + i[3]);
    p2 = ( L * (i[3] + i[4]) + R * (i[1] + i[2])) / (i[1] + i[2] + i[3] + i[4]);
    
    return (p1, p2, p3);


def compute_knot_insertion(p, i):
#  --i[0]--|--i[1]--|--i[2]--|--i[3]--|--i[4]--|--i[5]--
#         p[0]     p[1]    (p[2])    p[3]     p[4]
    p1 = p[0] * i[3]/(i[0]+i[1]+i[2]+i[3]) + p[1] * (i[0]+i[1]+i[2])/(i[0]+i[1]+i[2]+i[3]);
    p2 = p[1] * (i[3]+i[4]) / (i[1]+i[2]+i[3]+i[4]) + p[3] * (i[1]+i[2])/(i[1]+i[2]+i[3]+i[4]);
    p3 = p[3] * (i[3]+i[4]+i[5])/(i[2]+i[3]+i[4]+i[5]) + p[4] * (i[2])/(i[2]+i[3]+i[4]+i[5]);
    
    return (p1, p2, p3);


def build_knot_removal_data(dagPath):
    e = mt.selected_edges();
    ve = defaultdict(list)
    
    for eC in e:
        vI = mt.edge_vertices(dagPath, eC);
        ve[vI[0]].append(eC);
        ve[vI[1]].append(eC);
        
    for (v, e) in ve.items():
        if len(e) != 2:
            del ve[v];
            
    return ve;


def knot_modifier(modType):    
    (dagPath, component) = mt.selected_dagPath_component();
    
    meshFn = om.MFnMesh(dagPath);
    
    tccNode = mt.get_TCC_DGnode(dagPath);
    
    itv = cmds.getAttr(tccNode.name()+".itv");
    T = cmds.getAttr(tccNode.name()+".T");
    
    ve = build_knot_removal_data(dagPath);

    R = [];

    for (v,e) in ve.items():
        (he, v) = get_CV_line(dagPath, v, e[0], e[1]);

        i = [];
        for k in range(6):
            i.append(mt.get_halfedge_data(dagPath, itv, he[k]));
            
        p = [];
        for k in range(1,6):
            pk = om.MPoint(); meshFn.getPoint(v[k], pk); 
            p.append(Vec3(pk.x, pk.y, pk.z));

        if modType == KnotModifier.Insert:
            (p2, p3, p4) = compute_knot_insertion(p, i);
        elif modType == KnotModifier.Remove:
            (p2, p3, p4) = compute_knot_removal(p, i);
        else:
            raise Exception("Unknown Modifier");

        R.append( (v[2], p2) )
        R.append( (v[3], p3) )
        R.append( (v[4], p4) )
    
    return R;


def insert_knot():
    dagPath = mt.selected_dagPath();
    meshFn = om.MFnMesh(dagPath);
        
    for d in knot_modifier(KnotModifier.Insert):
            (v, p) = d;
            meshFn.setPoint(v, om.MPoint(p.x, p.y, p.z));

def remove_knot():
    dagPath = mt.selected_dagPath();
    meshFn = om.MFnMesh(dagPath);
        
    for d in knot_modifier(KnotModifier.Remove):
            (v, p) = d;
            meshFn.setPoint(v, om.MPoint(p.x, p.y, p.z));
