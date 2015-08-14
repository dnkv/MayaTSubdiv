//
//  updateTCC.cpp
//  MatlabKVCHDS
//
//  Created by Denis Kovacs on 9/15/12.
//  Copyright 2012 Denis Kovacs. All rights reserved.
//

#include <iostream>

#include "updateTCC.h"
#include "hds.h"
#include "varray.h"

using namespace std;

struct TCCTempData
{
    v_size_t valence;
    v_size_t TFvalence;

    v_bool   is_openTF; // bool array, true if T-face open
    v_size_t openTF;    // indices of open T-faces
    v_size_t openTtip;
    v_size_t openTH;


    v_bool mismatched;
    v_bool inconsistent;
    v_bool open_eqc;
    
    v_double eqm;
    v_bool eq_consistent;
};

void override_T(HDS &hds, TCCNodeData &nd);
v_size_t compute_valence(HDS &hds);
v_bool compute_open_Tfaces(HDS &hds, TCCNodeData &nd);
v_size_t compute_open_Tface_valence(HDS &hds, TCCNodeData &nd, TCCTempData &td);

bool heuristic1(HDS &hds, TCCNodeData &nd, TCCTempData &td);
bool heuristic2(HDS &hds, TCCNodeData &nd, TCCTempData &td);

void propagate_eqc(HDS &hds, TCCNodeData &nd, TCCTempData &td);
void recompute_equivalence_classes(HDS &hds, TCCNodeData &nd, TCCTempData &td);
void override_eqc(HDS &hds, TCCNodeData &nd);

void compute_eqc(HDS &hds, TCCNodeData &nd, TCCTempData &td);


void print_v_size_t(v_size_t &va)
{
    std::cout<<"["<<va.dim[0]<<" "<<va.dim[1]<<"]";
    for (size_t k=0; k<va.v.size(); k++)
    {
        std::cout<<" "<<va[k];
    }
    std::cout<<std::endl;
}

void print_v_uint64(v_uint64 &va)
{
    std::cout<<"["<<va.dim[0]<<" "<<va.dim[1]<<"]";
    for (size_t k=0; k<va.v.size(); k++)
    {
        std::cout<<" "<<va[k];
    }
    std::cout<<std::endl;
}

void print_v_bool(v_bool &va)
{
    std::cout<<"["<<va.dim[0]<<" "<<va.dim[1]<<"]";
    for (size_t k=0; k<va.v.size(); k++)
    {
        std::cout<<" "<<int(va[k]);
    }
    std::cout<<std::endl;
}

void print_v_double(v_double &va)
{
    std::cout<<"["<<va.dim[0]<<" "<<va.dim[1]<<"]";
    for (size_t k=0; k<va.v.size(); k++)
    {
        std::cout<<" "<<va[k];
    }
    std::cout<<std::endl;
}

void updateTCC(TCCNodeData &nd)
{
    HDS hds;
    TCCTempData td;
    
    hds.nFV  = nd.nFV;  
    hds.tip  = nd.tip;
    hds.V    = nd.V;
    hds.itv  = nd.itv;
    hds.T    = nd.T;
    hds.pole = nd.pole;
    
    finalize_HDS(hds);
    
    if (hds.nIHE() != hds.nHE())
    {
        v_size_t bI = v_range(hds.nIHE(), hds.nHE()-1);
        nd.eqc = v_cat<uint64_t>(v_col, nd.eqc, nd.eqc(hds.twin(bI)));
        hds.T = v_cat(v_col, hds.T, v_bool(1, length(bI), false));
        hds.itv = v_cat<double>(v_col, nd.itv, nd.itv(hds.twin(bI)));
    }
    
    
    // check first for bad topology
    nd.err.resizeDim(1, hds.nF());
    v_bool BT = (hds.nFV!=4) & (hds.nFV!=5);
    if (any(BT))
    {
        nd.err(BT) = 1;
        return;
    }

    //  remove T and pole data from faces with unassigned eqc/itv
    v_bool newHE = (nd.eqc == 0);
    hds.pole( hds.tip(newHE) ) = 2;
    hds.T(newHE)              = 0;

    //  override user-specified pole and T-tags
    override_T(hds, nd);
    //override_pole(Povr);

    //  compute valence
    td.valence = compute_valence(hds);
    
    //  mark all T-faces without tagged T-joints
    td.is_openTF  = compute_open_Tfaces(hds, nd);
    td.openTF = find(td.is_openTF);
    
    //  compute T-face valence (how many unassigned T-faces border a vertex?)
    td.TFvalence = compute_open_Tface_valence(hds, nd, td);
    
    v_size_t A = hds.F2H(td.openTF); v_transpose_1D(A);
    td.openTH = v_cat(v_col,v_cat(v_col,v_cat(v_col,v_cat(v_col, A,  A+1), A+2), A+3), A+4);
    td.openTtip = hds.tip(td.openTH);
    
    //  loop over open T-faces, find unique T-joints using two criteria
    while (1)
    {
        bool found = heuristic1(hds, nd, td) | heuristic2(hds, nd, td);
    
        if (!found || isempty(td.openTH)) break;
    }


    //  fill in T-joints through equivalence classes
    propagate_eqc(hds, nd, td);
    
    
    //  recompute eqc (and itv if prev. unassigned), set openTF/mismatched/inconsistent/open_eqc
    recompute_equivalence_classes(hds, nd, td);
    
    //  if edge given: update equivalence class relative to edge's itv value.
//    override_eqc(eqcOvr);
    
    //  fix obvious poles (valence > 4) or ((no nb open t-faces) and (valence != 4))
    hds.pole( (td.TFvalence==0 & td.valence!=4) | (td.valence > 4) ) = 1;
    hds.pole( td.valence ==4 ) = 0;

    //if !any(hds.open_eqc) //  all T-joints marked -> check for valid topology
    //     BT = is_bad_topology(hds);  //  BROKEN!
    //end

    td.open_eqc.resize(1, hds.nF());
    td.is_openTF.resize(1, hds.nF());
    td.mismatched.resize(1, hds.nF());
    td.inconsistent.resize(1, hds.nF());
    
    nd.err(td.open_eqc) = 4;
    nd.err(td.is_openTF)   = 3;
    nd.err(td.mismatched) = 2;
    nd.err(td.inconsistent) = 1;

    // all nd arrays include border halfedges
    nd.pole = hds.pole;
    nd.T = hds.T; nd.T.resize(1, hds.nIHE());
    nd.itv = hds.itv; nd.itv.resize(1, hds.nIHE());
    nd.eqc.resize(1, hds.nIHE());
}    


void override_T(HDS &hds, TCCNodeData &nd)
{
    if (isempty(nd.selHE)) return;
    for (size_t k=0; k<length(nd.selHE); k++) 
    {
        size_t kHE = nd.selHE[k];
        size_t cF = hds.H2F(kHE);
        if (hds.nFV(cF)!=5) continue;
        hds.T(v_range(0,4)+hds.F2H(cF)) = 0; //  wipe out T data on this face
        hds.T(kHE) = 1;                      //  set new T tag for this face
    }
}

void override_pole(HDS &hds, TCCNodeData &nd)
{
    if (isempty(nd.selV)) return;
    hds.pole(hds.pole == 2) = 0;
    hds.pole(nd.selV) = 1;
}

v_bool compute_open_Tfaces(HDS &hds, TCCNodeData &nd)
{
    v_bool oT = (hds.nFV==5);
    
    for (size_t kT = 0; kT < hds.nIHE(); kT++)
    {
        if (hds.T[kT])
        {
            oT[hds.H2F[kT]] = 0;
        }
    }
    
    return oT;
}


v_size_t compute_valence(HDS &hds)
{
    v_size_t val(1, size(hds.V, v_col), 0);
    for (size_t k=0; k<hds.nHE(); k++)
    {
        val[hds.tip[k]]++;
    }
    if (hds.nIHE()!=hds.nHE())
    {
        val(hds.tip(v_range(hds.nIHE(), hds.nHE()-1))) = 4;       //  border is all valence 4 for now
    }
    val(hds.tip(hds.T)) = val(hds.tip(hds.T)) + 1;   //  increment valence at T-joints
    
    return val;
}


v_size_t compute_open_Tface_valence(HDS &hds, TCCNodeData &nd, TCCTempData &td)
{
    v_size_t Tval(1, size(hds.V, v_col), 0);
    
    for (size_t k=0; k<length(td.openTF); k++)   // find(nd.openTF)
    {
        size_t kHE = hds.F2H[td.openTF[k]];
        for (size_t kFH = 0; kFH<5; kFH++) 
        {
            Tval[hds.tip(kHE + kFH)] += 1;
        }
    }
    if (hds.nIHE() != hds.nHE())
    {
        Tval(hds.tip(v_range(hds.nIHE(), hds.nHE()-1))) = 4;      //  borders are never a problem
    }
    return Tval;
}


bool heuristic1(HDS &hds, TCCNodeData &nd, TCCTempData &td)
{
    bool found = false;
    
    //  HEURISTIC 1: only one of the 5 T-face vertices is regular.
    // 
    v_bool oTreg  = (hds.pole(td.openTtip)==0) & (td.valence(td.openTtip)==4);
    v_bool uniqueT = bsum(oTreg, v_col) == 4; //  4 of the 5 vertices are regular?
    if (any(uniqueT))
    {
        found = true;
        
        v_size_t newTH   = td.openTH(uniqueT,_);
        v_size_t newTtip = td.openTtip(uniqueT,_);
        v_bool newTirreg = !oTreg(uniqueT,_);
        v_size_t newT = newTH(newTirreg);
        hds.T(newT) = 1;
        td.valence(hds.tip(newT)) = td.valence(hds.tip(newT)) + 1; //  increase valence at new T-joint
        td.TFvalence(newTtip) = td.TFvalence(newTtip) - 1; //  decrease Tface-valence for all vertices of new T-faces
        
        // remove these
        v_bool remainingTF = !uniqueT;
        td.openTH = td.openTH(remainingTF,_);
        td.openTtip = td.openTtip(remainingTF,_);
    }
    
    return found;
}

bool heuristic2(HDS &hds, TCCNodeData &nd, TCCTempData &td)
{
    bool found = false;
    //  HEURISTIC 2: non-pole vertices with valence 3 and one neighboring
    //               unassigned T-face
    v_bool uniqueV = (td.TFvalence==1) & (hds.pole==0) & (td.valence==3);
    if (any(uniqueV))
    {
        found = true;
        
        for (size_t kV=0; kV<length(uniqueV); kV++)
        {
            if (uniqueV[kV])
            {
                for (size_t rc = 0; rc<length(td.openTtip); rc++)
                {
                    if (td.openTtip[rc]==kV)
                    {
                        size_t r = rc / 5;
                        size_t c = rc % 5;
                        
                        hds.T(td.openTH(r,c)) = 1;
                        td.valence(kV) = td.valence(kV) + 1;
                        td.TFvalence(td.openTtip(r,_)) = td.TFvalence(td.openTtip(r,_)) - 1;
                        
                        
                        v_bool rmidx(size(td.openTH, v_row), 1, true);
                        rmidx(r) = false;
                        td.openTH = td.openTH(rmidx, _);
                        td.openTtip = td.openTtip(rmidx, _);
                    }
                }
            }
        }
    }
    
    return found;
}

void recompute_equivalence_classes(HDS &hds, TCCNodeData &nd, TCCTempData &td)
{    
    //  recompute equivalence classes (without changing itv)
    compute_eqc(hds, nd, td);

    td.mismatched = v_bool(1,hds.nF()+1, false);    // +1 = border "face"
    td.inconsistent = v_bool(1,hds.nF()+1, false); // +1 = border "face"
    td.open_eqc = v_bool(1,hds.nF()+1, false);     // +1 = border "face"
    
    
    v_size_t H2Fb(1, hds.nHE(), hds.nF());
    H2Fb(v_range(0, hds.nIHE()-1)) = hds.H2F;

    
    //  assign hds.eqm to equivalence classes that have been previously unassigned faces in same equivalence class as openTF:
    td.open_eqc(H2Fb(nd.eqc==0)) = true;
    
    //  assign itv to previously unassigned eqc, mark mismatched/inconsistent
    //  faces
    for (size_t k=0; k<length(td.eq_consistent); k++)
    {
        v_bool c_eqc = nd.eqc==(k+1);
        
        if (td.eq_consistent(k))
        {
            if (all( hds.itv(c_eqc)==0 ))
            {
                hds.itv(c_eqc) = td.eqm(c_eqc);
            } 
            else 
            {
                //  check if the currently assigned itv values match
                v_double ratios = hds.itv(c_eqc)/td.eqm(c_eqc);
                
                //  ignore unassigned itvs
                double maxRatio       = max(ratios);
                v_bool equal_ratios   = ratios == maxRatio;
                v_bool nonzero_ratios = ratios > 0;
                
                if (!all(equal_ratios(nonzero_ratios)))
                {
                    td.mismatched(H2Fb(c_eqc)) = true;
                } 
                else 
                {
                    //  fill in itv=0 eqc data
                    hds.itv(c_eqc) = maxRatio * td.eqm(c_eqc);
                }
            }
        } 
        else 
        {
            td.inconsistent(H2Fb(c_eqc)) = true;
        }
    }
    
}

void propagate_eqc(HDS &hds, TCCNodeData &nd, TCCTempData &td)
{
    v_size_t heQva = find( nd.eqc(v_range(0, hds.nIHE()-1))==0 & nd.eqc(hds.twin(v_range(0, hds.nIHE()-1))) > 0 );
    
    vector<size_t> heQ(length(heQva)), next_heQ;
    for (size_t kQ = 0; kQ<length(heQ); kQ++) heQ[kQ] = heQva[kQ];
    
    while (1)
    {
        bool do_next = false;
        for (size_t kQ = 0; kQ<length(heQ); kQ++)
        {
            size_t kHE = heQ[kQ];

            if (nd.eqc(kHE) != 0) continue;
            
            size_t kF = hds.H2F(kHE);
            size_t twin_eqc = nd.eqc(hds.twin(kHE));
            nd.eqc(kHE) = twin_eqc;
            
            if (hds.nFV(kF)==4) //  quad? know how to propagate
            {
                size_t oppHE = hds.next(hds.next(kHE));
                if (nd.eqc(oppHE) != 0)  //  if opposite edge is already assigned, use their eqc
                {
                    nd.eqc(kHE) = nd.eqc(oppHE); continue;
                }
                
                nd.eqc(oppHE) = twin_eqc;
                size_t twin_oppHE = hds.twin(oppHE);
                if (nd.eqc(twin_oppHE)==0 && twin_oppHE<=hds.nIHE())
                {
                    next_heQ.push_back(twin_oppHE);
                }
                do_next = true;
                continue;
            }
            
            v_size_t Tf_HEs = hds.F2H(kF)+v_range(0,4);
            v_size_t TposVA = find(hds.T(Tf_HEs));
            v_uint64 Tf_eqc = nd.eqc(Tf_HEs);
            
            int Tpos = 0;
            if (isempty(TposVA))  // T-face with no known T-joint?
            {
                v_uint64 uTf_eqc = unique<uint64_t>(Tf_eqc(Tf_eqc!=0));
                
                v_size_t counts(1, length(uTf_eqc), 0);
                for (size_t kEq=0; kEq<length(uTf_eqc); kEq++)
                {
                    counts[kEq] = bsum(Tf_eqc==uTf_eqc[kEq], v_col)[0];
                }
                    
                v_size_t kM = find(counts>1);
                if (isempty(kM) || length(counts)<2)  //  have to have at least two eqc, one with 2 edges to proceed
                {
                    continue; //  don't know how to propagate yet
                }
                
                td.is_openTF(kF) = 0;
                size_t itv_kM = uTf_eqc(kM[0]);
                size_t Tf_k = find(Tf_eqc==itv_kM)[0];
                
                // determine where the T-joint is (4 possibilities)
                if (Tf_eqc(Tf_k+1)==itv_kM)             //  [*1 1 0 0 0]
                {
                    Tpos = Tf_k;
                }
                else if (Tf_eqc(Tf_k+2)==itv_kM)        //  [1 0 1 *0 0]
                {
                    Tpos = (Tf_k + 3) % 5;
                }
                else if (Tf_eqc(Tf_k+3)==itv_kM)        //  [1 *0 0 1 0]
                {
                    Tpos = (Tf_k + 1) % 5;
                }
                else if (Tf_eqc(Tf_k+4)==itv_kM)        //  [1 0 0 0 *1]
                {
                    Tpos = (Tf_k - 1) % 5;
                }
                hds.T(hds.F2H(kF) + Tpos) = 1;
                
                size_t kV = hds.tip(hds.F2H(kF) + Tpos);
                td.valence(kV) = td.valence(kV) + 1;
                td.TFvalence(hds.tip(Tf_HEs)) = td.TFvalence(hds.tip(Tf_HEs)) - 1;
            }
            
            //  now T-face has known T-joint -> propagate all known edges
            
            v_bool nonTdir = circshift(cv_bool(false, false, true, false, true), Tpos, v_col);
            v_bool Tdir = !nonTdir;
            
            if (any(Tf_eqc(nonTdir)>0))
            {
                v_uint64 _Tf = Tf_eqc(nonTdir);
                Tf_eqc(nonTdir) = max(_Tf);
            }
            
            if (any(Tf_eqc(Tdir)>0))
            {
                v_uint64 _Tf = Tf_eqc(Tdir);
                Tf_eqc(Tdir) = max(_Tf);
            }
            
            nd.eqc(Tf_HEs) = Tf_eqc;
            v_size_t tHE = hds.twin(Tf_HEs(Tf_eqc>0));
            tHE = tHE(nd.eqc(tHE)==0 & tHE<=hds.nIHE());
            if (!isempty(tHE))
            {
                for (size_t ktHE=0; ktHE<length(tHE); ktHE++) 
                {
                    next_heQ.push_back(tHE[ktHE]);
                }
                do_next = true;
            }
        }
        
        if (!do_next) break;
        heQ = next_heQ;
    }
 
    td.openTF = find(td.is_openTF);
}

size_t propagate_eqc(HDS &hds, v_double &itv, size_t he0, v_bool &done, bool &consistent);

void compute_eqc(HDS &hds, TCCNodeData &nd, TCCTempData &td)
{
    
    td.eq_consistent.resize(1, 0);
    nd.eqc = v_uint64(1, hds.nHE());
    v_bool allDone(1, hds.nHE());
    
    v_double itv(1, hds.nHE(), 0);
    
    
    // fill out unassigned T-faces first and compute their equivalence class
    for (size_t k=0; k<length(td.openTF); k++)
    {
        allDone(hds.F2H(td.openTF[k])+v_range(0,4)) = true;
    }
    size_t nHE_done = 5*length(td.openTF);

    bool consistent;

    for (size_t k=0; k<length(td.openTF); k++)
    {
        v_size_t kHE = hds.F2H(td.openTF[k])+v_range(0,4);
        for (size_t kk=0; kk<5; kk++)
        {
            nHE_done += propagate_eqc(hds, itv, kHE[kk], allDone, consistent);
        }
    }
    
    
    std::vector<bool> eq_consistent;
    size_t eqcN = 0;
    size_t seed = 0;
    size_t nHE = hds.nHE();
    while (nHE_done<nHE)
    {
        eqcN = eqcN + 1;
        while (seed<nHE && allDone[seed]) seed++;
        v_bool done(1, hds.nHE(), false);
        nHE_done += propagate_eqc(hds, itv, seed, done, consistent);
        
        allDone(done) = true;
        if (consistent)
        {
            nd.eqc(done)=eqcN;
        }
        else
        {
            nd.eqc(done)=-eqcN;
        }
        eq_consistent.push_back(consistent);
    }
    
    td.eq_consistent.resize(1, eqcN);
    copy(eq_consistent.begin(), eq_consistent.end(), &td.eq_consistent[0]);
    
    
    td.eqm = itv;
    for (size_t k=0; k<eqcN; k++)
    {
        v_bool eq = nd.eqc==(k+1);
        if (eq_consistent[k])
        {
            td.eqm(eq) = td.eqm(eq) / min<double>(td.eqm(eq));
        }
        else
        {
            td.eqm(eq) = 0;
        }
    }

}

#define is_border_halfedge(h) (h>=hds.nIHE())

#define setItv(h, v) \
if (!done(h))\
{\
    done(h)=true; itv(h) = v;\
    nMarked = nMarked + 1;\
    if (!is_border_halfedge(h)) \
    {\
        Q.push_back(h);\
    }\
}\
else\
{\
    if (itv(h)!=v) \
    {\
        consistent = false;\
    }\
}

size_t propagate_eqc(HDS &hds, v_double &itv, size_t he0, v_bool &done, bool &consistent)
{
    size_t nMarked = 0;
    
    itv(he0) = 1;
    if (!done(he0))
    {
        done(he0) = 1;
        nMarked = nMarked + 1;
    }

    size_t iQ = 0;
    vector<size_t> Q; Q.push_back(he0);
    consistent = true;
    
    while (iQ<length(Q))
    {
        size_t he = Q[iQ]; 
        double cItv = itv(he);
        
        setItv(hds.twin(he), cItv);
        bool isHalf = false;
        
        // check neighbor edges on same side
        if (hds.T(hds.prev(he)))
        {
            setItv(hds.prev(he), cItv);
            isHalf = true;
        }
        
        if (hds.T(he))
        {
            he = hds.next(he);
            setItv(he, cItv);
            isHalf = true;
        }
        
        if (isHalf) cItv = 2*cItv;
        
        // go to opposite side of face
        he = hds.next(he); if (hds.T(he)) he = hds.next(he);
        he = hds.next(he);
        
        if (hds.T(he))
        {
            setItv(he, cItv / 2.0);
            setItv(hds.next(he), cItv / 2.0);
        }
        else
        {
            setItv(he, cItv);
        }
        
        iQ = iQ + 1;
    }

    return nMarked;
}

