//
//  hds.c
//  MSocketHDS
//
//  Created by Denis Kovacs on 3/31/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "hds.h"

#include <map>
#include <utility>
#include <iostream>
#include <algorithm>

using namespace std;

/*
void dumpArray(const char *name, std::vector<uint64_t> &v)
{
    std::cout<<name<<"["<<v.size()<<"]: ";
    for (int k=0; k<v.size(); k++) std::cout<<" "<<v[k]<<" ";
    std::cout<<std::endl;
}

void dumpRawArray(const char *name, RawArray<uint64_t> &v)
{
    std::cout<<v.p<<"["<<v.n<<"]: ";
    for (int k=0; k<v.n; k++) std::cout<<" "<<v.p[k]<<" ";
    std::cout<<std::endl;
}


void dumpHDS(HDS &hds)
{
    dumpArray("nFV", hds.nFV);
    dumpArray("F2H", hds.F2H);
    dumpArray("H2F", hds.H2F);
    dumpArray("next", hds.next);
    dumpArray("prev", hds.prev);
    dumpArray("twin", hds.twin);
    dumpArray("tip", hds.tip);
    dumpArray("V2H", hds.V2H);
}

void dumpHDS_R(HDS_R hds_r)
{
    dumpRawArray("nFV", hds_r.nFV);
    dumpRawArray("F2H", hds_r.F2H);
    dumpRawArray("H2F", hds_r.H2F);
    dumpRawArray("next", hds_r.next);
    dumpRawArray("prev", hds_r.prev);
    dumpRawArray("twin", hds_r.twin);
    dumpRawArray("tip", hds_r.tip);
    dumpRawArray("V2H", hds_r.V2H);
}
*/


typedef bool (*VV2Hcomp_ptr)(std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>);
bool VV2Hcomp (std::pair<uint64_t, uint64_t> a, std::pair<uint64_t, uint64_t> b) {
    if (a.first==b.first) return a.second<b.second;
    return a.first<b.first;
}

bool finalize_HDS(HDS &hds)
// build <next, prev, opp> HDS fields from <nFV, tip> (nFV is optional, default is 4)
{
    if (hds.tip.n==0) return false;
    bool has_nFV = hds.nFV.n > 0;
    size_t nH = hds.tip.n;
    size_t nF = has_nFV ? hds.nFV.n : nH / 4;
    size_t nV = hds.tip.v.max()+1;
    
    hds.next.resizeDim(1, nH, uint64_t(-1));
    hds.prev.resizeDim(1, nH, uint64_t(-1));
    hds.twin.resizeDim(1, nH, uint64_t(-1));
    
    hds.F2H.resizeDim(1, nF, uint64_t(-1));
    hds.H2F.resizeDim(1, nH, uint64_t(-1));
    hds.V2H.resizeDim(1, nV, uint64_t(-1));
    
    uint64_t cH = 0;
    for (size_t k=0; k<nF; k++)
    {
        hds.F2H[k] = cH;
        
        uint64_t nFV = has_nFV ? hds.nFV[k] : 4; // default nFV is 4
        for (size_t cFV=0; cFV<nFV; cFV++)
        {
            hds.next[cH+cFV]=(cFV==nFV-1) ? cH : cH+cFV+1 ;
            hds.prev[cH+cFV]=(cFV==0) ? cH+nFV-1 : cH+cFV-1;
            hds.H2F[cH+cFV]=k;
            hds.V2H[hds.tip[cH+cFV]]=cH+cFV;
        }
        cH+=nFV;
    }
    
    std::map< std::pair<uint64_t, uint64_t>, uint64_t, VV2Hcomp_ptr> VV2Hmap(VV2Hcomp);

    for (size_t k=0; k<nH; k++)
    {
        std::pair<uint64_t,uint64_t> vvTwin(hds.tip[k],hds.tip[hds.prev[k]]), vv(hds.tip[hds.prev[k]], hds.tip[k]);
                
        if ( VV2Hmap.find(vv)!=VV2Hmap.end())
        {
            std::cerr<<"Nonmanifold topology! Halfedge "<<hds.tip[hds.prev[k]]<<" - "<<hds.tip[k]<<" appears twice."<<std::endl;
            return false;
        }
        
        if ( VV2Hmap.find(vvTwin)==VV2Hmap.end() )
        {
            VV2Hmap[vv]=k;
        }
        else
        {
            VV2Hmap[vv]=k;
            uint64_t twin=VV2Hmap[vvTwin];
            hds.twin[k] = twin;
            hds.twin[twin] = k;
        }
    }
    
    size_t nHB = nH;

    // fix border halfedges
    
    for (size_t k=0; k<nH; k++)
    {
        if (hds.twin[k]==uint64_t(-1)) nHB++;
    }
    
    hds.tip.resizeDim(1, nHB, uint64_t(-1));
    hds.next.resizeDim(1, nHB, uint64_t(-1));
    hds.prev.resizeDim(1, nHB, uint64_t(-1));
    hds.twin.resizeDim(1, nHB);
    
    cH = nH;
    for (size_t k=0; k<nH; k++)
    {
        if (hds.twin[k]==uint64_t(-1))
        {
            hds.twin[k]=cH;
            hds.twin[cH]=k;
            hds.tip[cH] = hds.tip[ hds.prev[k] ];
            
            cH++;
        }
    }
    
    // fix prev/next loop at the tip vertex of each border halfedge
    for (size_t k=nH; k<nHB; k++)
    {
      
        uint64_t bd_next = hds.twin[ hds.prev[ hds.twin[k] ] ];
        while (hds.prev[bd_next]!=uint64_t(-1))
        {
            bd_next = hds.twin[ hds.prev[ bd_next ] ];
        }
        
        hds.next[k] = bd_next;
        hds.prev[bd_next] = k;
    }
    
    return true;
}



bool get_face_vertices(HDS_R &hds, uint64_t fI, KVCMutableIndexArray<uint64_t> &va)
{    
    if (fI>=hds.F2H.n) return false;
    uint64_t h=hds.F2H[fI], nFV = hds.nFV[fI];
    
    va.resizeDim(1, nFV);
    for (size_t k=0; k<nFV; k++)
    {
        va[k] = hds.tip[h];
        h = hds.next[h];
    }

    return true;
}

bool get_face_halfedges(HDS_R &hds, uint64_t fI, KVCMutableIndexArray<uint64_t> &he)
{    
    if (fI>=hds.F2H.n) return false;
    uint64_t h=hds.F2H[fI], nFV = hds.nFV[fI];
    
    he.resizeDim(1, nFV);
    for (size_t k=0; k<nFV; k++)
    {
        he[k] = h;
        h = hds.next[h];
    }
    
    return true;
}


bool get_halfedge_onering(HDS_R &hds, uint64_t h, KVCMutableIndexArray<uint64_t> &va)
{    
    if (h>=hds.tip.n) return false;
    uint64_t hS=h;
    
    std::vector<uint64_t> v;
    do 
    {
        v.push_back(hds.tip[h]);
        h = hds.next[h];
    }
    while(h!=hS);
    
    hS=hds.twin[hS];
    h=hds.next[hS];
    hS=hds.prev[hS];
    
    do {
        v.push_back(hds.tip[h]);
        h = hds.next[h];
    } while (h!=hS);
    
    va.resizeDim( 1, v.size() ); copy( v.begin(), v.end(), &va[0] );  
    return true;
}

bool get_vertex_onering(HDS_R &hds, uint64_t vI, KVCMutableIndexArray<uint64_t> &va)
{    
    if (vI>=hds.V2H.n) return false;
    uint64_t h = hds.twin[hds.V2H[vI]], hS=h;
    
    std::vector<uint64_t> v;
    while (1)
    {
        v.push_back(hds.tip[h]);
        h = hds.next[h];
        if (hds.tip[h]==vI)
        {
            h = hds.twin[h];
            if (h==hS) break; 
            h = hds.next[h];
        }
    }
    v.pop_back();
    
    va.resizeDim( 1, v.size() ); copy( v.begin(), v.end(), &va[0] );  
    return true;
}

bool get_vertex_halfedges(HDS_R &hds, uint64_t vI, KVCMutableIndexArray<uint64_t> &he)
{    
    if (vI>=hds.V2H.n) return false;
    uint64_t h = hds.V2H[vI], hS=h;

    std::vector<uint64_t> heV;
    while (1)
    {
        heV.push_back(h);
        h = hds.twin[hds.next[h]];
        if (h==hS) break;
    }
    
    he.resizeDim( 1, heV.size() ); copy( heV.begin(), heV.end(), &he[0] );  
    return true;
}

bool get_halfedge_halfedges(HDS_R &hds, uint64_t h, KVCMutableIndexArray<uint64_t> &he)
{    
    if (h>=hds.tip.n) return false;
    uint64_t hS=h;

    std::vector<uint64_t> heV;
    while (1)
    {
        heV.push_back(h);
        h = hds.twin[hds.next[h]];
        if (h==hS) break;
    }
    
    he.resizeDim( 1, heV.size() ); copy( heV.begin(), heV.end(), &he[0] );  
    return true;
}


bool subdivide_HDS(HDS &hds, HDS &hds_sd)
{
    bool has_nFV = hds.nFV.n > 0;
    
    size_t nF = hds.F2H.n;
    size_t nV = hds.V2H.n;
    size_t nHB = hds.twin.n;
    
    size_t nV_prev = nV;
    
    std::vector<uint64_t> v_nFV;
    std::vector<uint64_t> v_tip;
    std::vector<double>   v_itv;
    std::vector<bool>     v_T;
    
    std::vector<uint64_t> H2EV(nHB, -1); // halfedge-to-edge vertex map

    // nF new vertices for facet centers
    nV+=nF;
    
    // nE new vertices for edge centers
    for (size_t iH=0; iH<nHB; iH++)
    {
        if (H2EV[iH]==uint64_t(-1)) 
        {
            H2EV[iH]=H2EV[hds.twin[iH]]=nV; nV++;
        }
    }
    
    
    /* Picture for face creation:
     *    *-------T1------*           *---*---*-------*        *-------*-------*
     *    ^                           |       |                |       |
     * cH>|                           *       |                |       |
     *    |                           |       |                |       |
     *    T2                  ->      *-------*           or   *-------*
     *    |                           |                        |
     *    |                           |                        |
     *    |                           |                        |
     *    *                           *                        *
     */
    
    uint64_t cH=0, new_nFV, ncH, pcH;
    for (size_t iF=0; iF<nF; iF++)
    {
        size_t nFV = has_nFV ? hds.nFV[iF] : 4;
        for (size_t iFV=0; iFV<nFV; iFV++, cH++)
        {
            if (hds.T[cH]) continue;
            ncH = hds.next[cH];
            pcH = hds.prev[cH];
            bool isT1 = hds.T[ncH];
            bool isT2 = hds.T[pcH];

            new_nFV = 4;
            v_tip.push_back( hds.tip[cH] );
            v_tip.push_back( H2EV[ncH] );
            if (isT1)
            {
                v_tip.push_back( hds.tip[ncH] ); new_nFV++;
            }
            v_tip.push_back( nV_prev + iF );
            
            if (isT2)
            {
                v_tip.push_back( hds.tip[pcH] ); new_nFV++;
            }
            
            v_tip.push_back( H2EV[cH] );             
            

            // store T-joint flags 
            v_T.push_back(0); 
            v_T.push_back(hds.T[ncH]); 
            if (isT1) v_T.push_back(0); 
            v_T.push_back(0); 
            if (isT2) v_T.push_back(0); 
            v_T.push_back(hds.T[pcH]);
            
            
            // store knot intervals
            char b = isT1; if (isT2) b+=2;
            switch (b)
            {
                case 0:
                    /*
                     *    1-------2-------*
                     *    |       |      
                     *    |       |      
                     *    |       |      
                     *    4-------3      
                     *    |              
                     *    |              
                     *    |              
                     *    *              
                     */
                    v_itv.push_back(hds.itv[cH]/2.0);
                    v_itv.push_back(hds.itv[ncH]/2.0);
                    v_itv.push_back(hds.itv[cH]/2.0);
                    v_itv.push_back(hds.itv[ncH]/2.0);
                    break;
                case 1:
                    /*
                     *    1---2---3-------*
                     *    |       |      
                     *    |       |      
                     *    |       |      
                     *    5-------4      
                     *    |              
                     *    |              
                     *    |              
                     *    *              
                     */
                    v_itv.push_back(hds.itv[cH]/2.0);
                    v_itv.push_back(hds.itv[ncH]/2.0);
                    v_itv.push_back(hds.itv[ncH]/2.0);
                    v_itv.push_back(hds.itv[cH]/2.0);
                    v_itv.push_back(hds.itv[ncH]);
                    break;
                case 2:
                    /*
                     *    1-------2-------*
                     *    |       |      
                     *    5       |      
                     *    |       |      
                     *    4-------3      
                     *    |              
                     *    |              
                     *    |              
                     *    *              
                     */
                    v_itv.push_back(hds.itv[cH]/2.0);
                    v_itv.push_back(hds.itv[ncH]/2.0);
                    v_itv.push_back(hds.itv[cH]);
                    v_itv.push_back(hds.itv[ncH]/2.0);
                    v_itv.push_back(hds.itv[cH]/2.0);
                    break;
                case 3:
                    /*
                     *    1---2---3-------*
                     *    |       |      
                     *    6       |      
                     *    |       |      
                     *    5-------4      
                     *    |              
                     *    |              
                     *    |              
                     *    *              
                     */
                    v_itv.push_back(hds.itv[cH]/2.0);
                    v_itv.push_back(hds.itv[ncH]/2.0);
                    v_itv.push_back(hds.itv[ncH]/2.0);
                    v_itv.push_back(hds.itv[cH]);
                    v_itv.push_back(hds.itv[ncH]);
                    v_itv.push_back(hds.itv[cH]/2.0);
                    break;
            }
            
            v_nFV.push_back(new_nFV);
        }
    }

    hds_sd.nFV.resizeDim(1, v_nFV.size()); copy( v_nFV.begin(), v_nFV.end(), &hds_sd.nFV[0] );  
    hds_sd.tip.resizeDim(1, v_tip.size()); copy( v_tip.begin(), v_tip.end(), &hds_sd.tip[0] );  
    hds_sd.T.resizeDim(1, v_T.size());     copy( v_T.begin(), v_T.end(), &hds_sd.T[0] );  
    hds_sd.itv.resizeDim(1, v_itv.size()); copy( v_itv.begin(), v_itv.end(), &hds_sd.itv[0] );  

    finalize_HDS(hds_sd);
    
    return true;
}

