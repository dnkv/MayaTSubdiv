//
//  hds.c
//  MSocketHDS
//
//  Created by Denis Kovacs on 3/31/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include "hds.h"

#include <map>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <algorithm>
using namespace std;
// using namespace tr1;
#include <hds_hash.h>


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


void test_hash()
{
    std::unordered_map<uint64_t ,double> mymap;
    
    mymap[10] = 10.0;
    
    uint64_t input = 10;
    
    std::unordered_map<uint64_t,double>::const_iterator got = mymap.find (input);
    
    if ( got == mymap.end() )
        std::cout << "not found";
    else
        std::cout << got->first << " is " << got->second;
    
    std::cout << std::endl;
}





bool finalize_HDS(HDS &hds)
// build <next, prev, opp> HDS fields from <nFV, tip> (nFV is optional, default is 4)
{
    if (hds.tip.size()==0) return false;
    bool has_nFV = hds.nFV.size() > 0;
    size_t nH = hds.nHE();
    size_t nF = hds.nF();
    size_t nV = hds.tip.v.max()+1; // V might not exist
    
    hds.next.resize(1, nH, uint64_t(-1));
    hds.prev.resize(1, nH, uint64_t(-1));
    hds.twin.resize(1, nH, uint64_t(-1));
    
    hds.F2H.resize(1, nF, uint64_t(-1));
    hds.H2F.resize(1, nH, uint64_t(-1));
    hds.V2H.resize(1, nV, uint64_t(-1));
    
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
    
//    std::map< std::pair<uint64_t, uint64_t>, uint64_t, VV2Hcomp_ptr> VV2Hmap(VV2Hcomp);
    
    typedef std::pair<uint64_t, uint64_t> VVpair;
    std::unordered_map< VVpair, uint64_t> VV2Hmap;

    for (size_t k=0; k<nH; k++)
    {
        uint64_t v0=hds.tip[k], v1=hds.tip[hds.prev[k]];
        bool v01 = v0<v1;

        VVpair vv(v01 ? v0:v1, v01 ? v1:v0);
        
        std::unordered_map<VVpair, uint64_t>::const_iterator got = VV2Hmap.find(vv);

        if ( got==VV2Hmap.end() )
        {
            VV2Hmap[vv]=k;
        }
        else
        {
            uint64_t twin = got->second;
            if (hds.twin[twin]!=uint64_t(-1))
            {
                std::cerr<<"Nonmanifold topology! Halfedge "<<hds.tip[hds.prev[k]]<<" - "<<hds.tip[k]<<" appears twice."<<std::endl;
                return false;
            }

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
    
    hds.tip.resize(1, nHB, uint64_t(-1));
    hds.next.resize(1, nHB, uint64_t(-1));
    hds.prev.resize(1, nHB, uint64_t(-1));
    hds.twin.resize(1, nHB);
    
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



bool get_face_vertices(HDS_R &hds, uint64_t fI, KVCIndexArray<uint64_t> &va)
{    
    if (fI>=hds.nF()) return false;
    uint64_t h=hds.F2H[fI], nFV = hds.nFV[fI];
    
    va.resize(1, nFV);
    for (size_t k=0; k<nFV; k++)
    {
        va[k] = hds.tip[h];
        h = hds.next[h];
    }

    return true;
}

bool get_face_halfedges(HDS_R &hds, uint64_t fI, KVCIndexArray<uint64_t> &he)
{    
    if (fI>=hds.nF()) return false;
    uint64_t h=hds.F2H[fI], nFV = hds.nFV[fI];
    
    he.resize(1, nFV);
    for (size_t k=0; k<nFV; k++)
    {
        he[k] = h;
        h = hds.next[h];
    }
    
    return true;
}


bool get_halfedge_onering(HDS_R &hds, uint64_t h, KVCIndexArray<uint64_t> &va)
{    
    if (h>=hds.nHE()) return false;
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
    
    va.resize( 1, v.size() ); copy( v.begin(), v.end(), &va[0] );
    return true;
}

bool get_vertex_onering(HDS_R &hds, uint64_t vI, KVCIndexArray<uint64_t> &va)
{    
    if (vI>=hds.nV()) return false;
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
    
    va.resize(1, v.size() ); copy( v.begin(), v.end(), &va[0] );
    return true;
}

bool get_vertex_halfedges(HDS_R &hds, uint64_t vI, KVCIndexArray<uint64_t> &he)
{    
    if (vI>=hds.nV()) return false;
    uint64_t h = hds.V2H[vI], hS=h;

    std::vector<uint64_t> heV;
    while (1)
    {
        heV.push_back(h);
        h = hds.twin[hds.next[h]];
        if (h==hS) break;
    }
    
    he.resize( 1, heV.size() ); copy( heV.begin(), heV.end(), &he[0] );
    return true;
}

bool get_halfedge_halfedges(HDS_R &hds, uint64_t h, KVCIndexArray<uint64_t> &he)
{    
    if (h>=hds.nHE()) return false;
    uint64_t hS=h;

    std::vector<uint64_t> heV;
    while (1)
    {
        heV.push_back(h);
        h = hds.twin[hds.next[h]];
        if (h==hS) break;
    }
    
    he.resize( 1, heV.size() ); copy( heV.begin(), heV.end(), &he[0] );  
    return true;
}


bool subdivide_HDS(HDS &hds, HDS &hds_sd)
{
    size_t nF = hds.nF();
    size_t nV = hds.nV();
    size_t nHB = hds.nHE();
    
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
        size_t nFV = hds.nFV[iF];
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

    hds_sd.nFV.resize(1, v_nFV.size()); copy( v_nFV.begin(), v_nFV.end(), &hds_sd.nFV[0] );
    hds_sd.tip.resize(1, v_tip.size()); copy( v_tip.begin(), v_tip.end(), &hds_sd.tip[0] );
    hds_sd.T.resize(1, v_T.size());     copy( v_T.begin(), v_T.end(), &hds_sd.T[0] );
    hds_sd.itv.resize(1, v_itv.size()); copy( v_itv.begin(), v_itv.end(), &hds_sd.itv[0] );  

    finalize_HDS(hds_sd);
    
    return true;
}

