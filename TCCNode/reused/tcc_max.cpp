//
//  hds.c
//  MSocketHDS
//
//  Created by Denis Kovacs on 3/31/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

/* MATLAB Conversion steps:
 *  - operator=: both sides have to have the same length!
 *  - replace % with //
 *  - replace ~ with !
 *  - indexing: 1-based -> 0-based (start with replacing 1->0, then 2->1, etc)
 *  - A(:,1) = X        ->  A.rcol(0) = X
 *  - X = A(:,1)        ->  X = A.col(0)
 *  - idx > v.size()    ->  idx >= v.size()
 *  - make sure using modulo that the argument is not negative with unsigned ints!
 */


#include "hds.h"

#include <map>
#include <utility>
#include <iostream>
#include <algorithm>
#include "varray.h"
#include "tcc_max.h"

using namespace std;

namespace TCC_MAX
{
    
    const size_t MAX_SIZE_T = size_t(-1);
    
    void p_vdouble(v_double &v)
    {
        cout<<"["<<v.dim[0]<<","<<v.dim[1]<<"] ";
        for (size_t k=0; k<v.v.size(); k++)
        {
            cout<<v.v[k]<<" ";
        }
        cout<<endl;
    }
    
    void p_vsizet(v_size_t &v)
    {
        cout<<"["<<v.dim[0]<<","<<v.dim[1]<<"] ";
        for (size_t k=0; k<v.v.size(); k++)
        {
            cout<<v.v[k]<<" ";
        }
        cout<<endl;
    }
    
    void p_vbool(v_bool &v)
    {
        cout<<"["<<v.dim[0]<<","<<v.dim[1]<<"] ";
        for (size_t k=0; k<v.v.size(); k++)
        {
            cout<<int(v.v[k])<<" ";
        }
        cout<<endl;
    }
    
    
    
    
    bool is_border_halfedge(HDS &hds, size_t he)
    {
        return he>=size(hds.H2F, v_col);
    };
    
    bool is_border_edge(HDS &hds, size_t he)
    {
        return (he>size(hds.H2F, v_col)) | (hds.twin(he)>size(hds.H2F, v_col));
    }
    
    v_bool is_border_halfedge(HDS &hds, v_size_t he)
    {
        return he>=size(hds.H2F, v_col);
    };
    
    v_bool is_border_edge(HDS &hds, v_size_t he)
    {
        return (he>size(hds.H2F, v_col)) | (hds.twin(he)>size(hds.H2F, v_col));
    }
    
    bool is_corner(HDS &hds, size_t he) // is vertex pointed to by he a corner?
    {
        size_t v = hds.tip(he);
        if ((hds.corner.size() <= v) || (hds.corner[v]==2)) return hds.twin(hds.next(he)) == hds.prev(hds.twin(he));
        return (hds.corner[v]==1);
    }
    
    
    void copyHDS(HDS &a, HDS &b) // copies from a to b
    {
        b.nFV = a.nFV;
        b.F2H = a.F2H;
        
        b.H2F = a.H2F;
        b.next = a.next;
        b.prev = a.prev;
        b.twin = a.twin;
        
        b.tip = a.tip;
        b.tipK = a.tipK;
        
        b.T = a.T;
        b.pole = a.pole;
        b.val  = a.val;
        
        b.itv = a.itv;
        
        b.V2H = a.V2H;
        
        b.V = a.V;
        b.VK = a.VK;
        b.VFK = a.VFK;
//        b.UV = a.UV;
    }
    
    void extend_homogeneous(v_double &V)
    {
        size_t nR = V.nR(), nC = V.nC();
        
        v_double oV(V); // store original
        V.v.resize((nR+1) * nC); V.dim[0]=nR+1; // don't use varray.resize - no need to copy
        V(v_range(0,nR-1), _) = oV;
        V(nR, _) = 1;
    }
    
    void normalize_homogeneous(v_double &V)
    {
        size_t nR = V.nR();
        
        v_double Vnor = V(nR-1,_);
        
        for (size_t k=0; k<nR-1; k++)
        {
            V(k,_) /= Vnor;
        }
        V(nR-1,_) = 1;
    }
    
    void reduce_homogeneous(v_double &V) // assumes already normalized!
    {
        size_t nR = V.nR(), nC = V.nC();
        
        v_double oV(V);
        V.v.resize((nR-1) * nC); V.dim[0]=nR-1;
        V = oV(v_range(0,nR-2), _);
    }
    
    
    void compute_valence(HDS &hds)
    {
        size_t nV   = hds.V.nC();
        size_t nHE  = hds.tip.nC();
        size_t nIHE = hds.H2F.nC();
        
        hds.val = v_double(1, nV);
        
        for (size_t iHE=0; iHE<nHE; iHE++)
        {
            size_t v = hds.tip(iHE);
            hds.val(v)++;
            if (hds.T(iHE)) hds.val(v)++;
        }
        
        // for now all boundary vertices are regular!
        for (size_t iHE=nIHE; iHE<nHE; iHE++)
        {
            size_t v = hds.tip(iHE);
            hds.val(v) = 4;
        }
    }
    
    void subdivide(HDS &hds, size_t nSubdivs)
    {
        compute_valence(hds);
        
        // add homogeneous coordinate
        extend_homogeneous(hds.V);
        
        for (size_t iSD = 0; iSD<nSubdivs; iSD++)
        {
            HDS hds2;
            
            hds.T(v_range(length(hds.H2F), length(hds.tip)-1)) = false;
//            cout<<"Level: "<<iSD<<endl;
            
            subdivide_HDS(hds, hds2);
            
            v_size_t hds2BorderRange = v_range(length(hds2.T), length(hds2.tip)-1);
            v_double border_itv = hds2.itv(hds2.twin(hds2BorderRange));
            hds2.itv.resize(1, length(hds2.tip)); hds2.itv(hds2BorderRange) = border_itv;
            hds2.T.resize(1, length(hds2.tip), false);
            hds2.V.resize(size(hds.V, v_row), length(hds2.V2H));
            hds2.pole = hds.pole; hds2.pole.resize(1, length(hds2.V2H), false);
            hds2.val  = hds.val;  hds2.val.resize(1, length(hds2.V2H), 4);
            
            size_t nF  = length(hds.nFV);
            size_t nV  = size(hds.V, v_col);
            size_t nH  = length(hds.tip);
//            size_t nV2 = size(hds2.V,v_col);
            
            v_double V(hds2.V.nR(), hds2.V.nC());
            
//            hds2.UV.resize(2, nV2);
//            hds2.UV(_, v_range(0, nV-1)) = hds.UV;
            
            // build halfedge-to-edge and border face data structures
            v_size_t H2E(1, nH);
            size_t nE=0;
            
            for (size_t iH=0; iH<nH; iH++)
            {
                if (hds.twin(iH)<iH) { H2E(iH)=H2E(hds.twin(iH)); continue; }
                
//                hds2.UV(_,nV + nF + nE) = 0.5 * (hds.UV(_,hds.tip(iH)) + hds.UV(_,hds.tip(hds.twin(iH))));
                
                H2E(iH)=nE; nE++;
            }
            
            compute_Tknot_insertion(hds, hds.VK, hds.VFK, hds.tipK);
            
            for (size_t fI = 0; fI<nF; fI++)
            {
                v_size_t fH = hds.F2H(fI)+v_range(0, hds.nFV(fI)-1); // fH = face_halfedges(hds, fI);
                
                v_double pV = V(_, hds.tip(fH));
                v_double pF = V(_, nV + fI);
                
                v_size_t hef = nV + nF + H2E(fH);
                v_double pE = V(_, nV + nF + H2E(fH));
                
                compute_TCC_max_stencils(hds, fI, pF, pE, pV);
                
                V(_,hds.tip(fH)) = pV;
                V(_, nV + fI) = pF;
                V(_, nV + nF + H2E(fH)) = pE;
                
//                v_double UVcorners = hds.UV(_,hds.tip(fH(!hds.T(fH))));
//                hds2.UV(_, nV + fI) = (1.0/4.0) * sum(UVcorners, v_col);
            }
            
            normalize_homogeneous(V);
            
            // V = (3/N) * V + (1 - 3/N) * hds.V
            for (size_t k=0; k<nV; k++)
            {
                double w = 3.0 / hds.val(k);
                V(_, k) = V(_, k) * w + hds.V(_, k) * (1.0 - w);
            }
            
            if (length(hds.H2F)<length(hds.tip))
            {
                //knot doubling along the boundary
                //compute border edge midpoints
                v_bool bheDone(1, length(hds.tip));
                for (size_t k=length(hds.H2F); k<length(hds.tip); k++)
                {
                    if (bheDone(k)) continue;
                    
                    size_t heS = k;
                    v_size_t he = cv_size_t(hds.prev(heS), heS, hds.next(heS));
                    while (1)
                    {
                        v_double s = hds.itv(he);
                        if (is_corner(hds, he(0))) { s(0)=s(1); }
                        if (is_corner(hds, he(1))) { s(2)=s(1); }
                        V(_, nV + nF + H2E(he(1))) = ((s(0) + s(1)/2.0) * hds.V(_,hds.tip(he(1))) + (s(1)/2.0 + s(2))*hds.V(_,hds.tip(he(0))) ) / sum(s);
                        bheDone(he(1)) = true;
                        he = cv_size_t(he(1), he(2), hds.next(he(2)));
                        if (he(1)==heS) break;
                    }
                    
                    // compute border vertices
                    heS = k; he=cv_size_t(heS, hds.next(heS));
                    while (1)
                    {
                        v_double s = hds.itv(he);
                        size_t iV = hds.tip(he(0));
                        if (is_corner(hds, he(0)))
                        {
                            V(_,iV) = hds.V(_,iV);
                        } else {
                            V(_,iV) = (s(0) * V(_,nV+nF+H2E(he(1))) + (s(0)+s(1)) * hds.V(_,iV) + s(1) * V(_,nV+nF+H2E(he(0)))) / (2.0*sum(s));
                        }
                        he = cv_size_t(he(1), hds.next(he(1)));
                        if (he(0)==heS) { break; }
                    }
                }
            }
            
            hds2.V = V;
            
            copyHDS(hds2, hds);
        }
        
        reduce_homogeneous(hds.V);
    }
    
    void linear_subdivide(HDS &hds, size_t nSubdivs)
    {
        for (size_t iSD = 0; iSD<nSubdivs; iSD++)
        {
            HDS hds2;
            
            hds.T(v_range(length(hds.H2F), length(hds.tip)-1)) = false;
//            cout<<"Level: "<<iSD<<endl;
            
            subdivide_HDS(hds, hds2);
            
            v_size_t hds2BorderRange = v_range(length(hds2.T), length(hds2.tip)-1);
            v_double border_itv = hds2.itv(hds2.twin(hds2BorderRange));
            hds2.itv.resize(1, length(hds2.tip)); hds2.itv(hds2BorderRange) = border_itv;
            hds2.T.resize(1, length(hds2.tip), false);
            hds2.V.resize(size(hds.V, v_row), length(hds2.V2H));
            
            size_t nF  = length(hds.nFV);
            size_t nV  = size(hds.V, v_col);
            size_t nH  = length(hds.tip);
            //            size_t nV2 = size(hds2.V,v_col);
            
            v_double V(hds2.V.nR(), hds2.V.nC());
            
            for (size_t vI = 0; vI<nV; vI++)
            {
                V(_, vI) = hds.V(_, vI);
            }
            
            size_t nE=0;
            
            for (size_t iH=0; iH<nH; iH++)
            {
                if (hds.twin(iH)<iH) { continue; }
                
                V(_,nV + nF + nE) = 0.5 * (hds.V(_,hds.tip(iH)) + hds.V(_,hds.tip(hds.twin(iH))));
                
                nE++;
            }
            

            for (size_t fI = 0; fI<nF; fI++)
            {
                v_size_t fH = hds.F2H(fI)+v_range(0, hds.nFV(fI)-1); // fH = face_halfedges(hds, fI);
                
                v_double UVcorners = V(_,hds.tip(fH(!hds.T(fH))));
                V(_, nV + fI) = (1.0/4.0) * sum(UVcorners, v_col);
            }
            
            hds2.V = V;
            
            copyHDS(hds2, hds);
        }
    }
    
    // compute_Tknot_insertion.m ------------------------------------------------
    void process_Tline(HDS &hds, size_t TF, bool isLoop, v_double &VV, size_t &nVV, v_double &EV, v_size_t &Ktip);
    
    bool compute_Tknot_insertion(HDS &hds, v_double &VV, v_double &EV, v_size_t &Ktip)
    {
        size_t nVV = 0;
        VV.resize(size(hds.V, v_row), size(hds.V, v_col));
        EV.resize(size(hds.V, v_row), size(hds.nFV, v_col));
        Ktip.resize(1, size(hds.tip, v_col), size_t(-1));
        
        // find / process T-lines
        v_bool Tdone(1,size(hds.T, v_col));
        size_t nIH = size(hds.H2F, v_col);
        
        v_size_t v_he=find(hds.T);
        for (size_t k=0; k<size(v_he, v_col); k++)
        {
            size_t he = v_he[k];
            
            if (he<nIH && !Tdone(he))
            {
                bool isLoop = false;
                Tdone(he)=1;
                size_t heN = he;
                size_t FN = hds.H2F(he);
                while (!isLoop)
                {
                    heN = hds.next(hds.twin(hds.next(hds.next(heN))));
                    Tdone(heN) = 1;
                    if (is_border_halfedge(hds,heN) || !hds.T(heN)) { break; }
                    FN = hds.H2F(heN);
                    if (heN==he){ isLoop = true; }
                }
                heN = he;
                while (!isLoop)
                {
                    heN = hds.prev(hds.prev(hds.twin(hds.prev(heN))));
                    if (is_border_halfedge(hds,heN) || !hds.T(heN)) { break; }
                    Tdone(heN) = 1;
                    if (heN==he) { isLoop = true; }
                }
                
                process_Tline(hds, FN, isLoop, VV, nVV, EV, Ktip);
            }
        }
        
        VV.resize(size(VV, v_row), nVV);
        
        return true;
    }
    
    void process_Tline(HDS &hds, size_t TF, bool isLoop, v_double &VV, size_t &nVV, v_double &EV, v_size_t &Ktip)
    {
        size_t F = TF;
        
        //initialize knot intervals
        v_size_t FH = hds.F2H(F) + v_range(0,hds.nFV(F)-1); // face_halfedges(hds, F);
        v_bool T = hds.T(FH);
        v_size_t vTI = find(T); size_t TI = vTI[0];
        
        size_t heM = hds.prev(hds.prev(FH(TI)));
        size_t heA = hds.prev(hds.twin(hds.prev(heM)));
        size_t heB = hds.next(hds.twin(hds.next(heM)));
        
        v_size_t HES = cv_size_t(heA, heM, heB);
        
        // compute edge vertices
        v_double s = hds.itv(HES);
        if (is_border_halfedge(hds, heA)) { s(0)=s(1); }
        if (is_border_halfedge(hds, heB)) { s(2)=s(1); }
        
        while (1)
        {
            EV(_,F) = ((s(0) + s(1)/2.0) * hds.V(_,hds.tip(heM)) + (s(1)/2.0 + s(2)) * hds.V(_,hds.tip(heA)) ) / sum(s);
            
            // slide window to next T-face
            heA = heM;
            heM = heB;
            heB = hds.next(hds.twin(hds.next(heB)));
            if (is_border_halfedge(hds, heM)) { break; }
            s = cv_double(s(1), s(2), hds.itv(heB));
            if (is_border_halfedge(hds, heB)) { s(2)=s(1); }
            F = hds.H2F(heM);
            if (hds.nFV(F)==4 || hds.T(heM)) { break; }
            if (F==TF) { break; } // skip if loop
        }
        
        
        //update T-line interior original vertex positions
        v_size_t he=HES;
        s = hds.itv(he);
        
        if (!isLoop)
        {
            //only one knot insertion for first vertex
            if (is_border_halfedge(hds, he(0)))
            {
                s(0)=s(1);
                // if at border, copy border vertex and link to it
                VV(_, nVV) = hds.V(_,hds.tip(hds.prev(he(1)))); Ktip( hds.prev(he(1)) ) = nVV; nVV++;
            }
            else
            {
                size_t he0;
                if (hds.T(hds.prev(he(0))))
                {
                    he0 = hds.prev(he(0));
                }
                else
                {
                    he0 = hds.prev(hds.twin(hds.prev(he(0))));
                    //he0 might be a missing edge!
                    if (hds.T(he0)) { he0 = hds.prev(he0); }
                }
                double s0 = hds.itv(he0);
                
                VV(_, nVV) = ( (s0 + s(0) + s(1)/2.0) * hds.V(_,hds.tip(he(0))) + (s(1)/2.0)*hds.V(_,hds.tip(he0)) ) / (s0 + s(0) + s(1));
                Ktip( cv_size_t(he(0), hds.prev(he(1))) ) = nVV; nVV++;
            }
        }
        
        
        // knot doubling for all interior vertices
        if (!is_border_halfedge(hds,he(2)))
        {
            size_t F1 = TF;
            size_t F2 = hds.H2F(he(2));
            while (hds.nFV(F2)==5 && !hds.T(he(2)))
            {
                VV(_, nVV) = (s(1)*EV(_,F2) + (s(1)+s(2))*hds.V(_,hds.tip(he(1))) + s(2)*EV(_,F1)) / (2*(s(1)+s(2)));
                Ktip( cv_size_t(he(1), hds.twin(hds.next(he(1)))) ) = nVV; nVV++;
                
                he = cv_size_t(he(1), he(2), hds.next(hds.twin(hds.next(he(2)))));
                if (is_border_halfedge(hds, he(2))) { break; }
                s = cv_double(s(1), s(2), hds.itv(he(2)));
                F1 = F2;
                F2 = hds.H2F(he(2));
                
                if (F1==TF) { break; }
            }
        }
        
        if (!isLoop)
        {
            if (is_border_halfedge(hds, he(2)))
            {
                VV(_, nVV) = hds.V(_,hds.tip(he(1))); Ktip( cv_size_t(he(1), hds.twin(hds.next(he(1)))) ) = nVV; nVV++;
            }
            else
            {
                size_t he4;
                if (hds.T(he(2)))
                {
                    he4 = hds.next(he(2));
                }
                else
                {
                    he4 = hds.twin(hds.next(he(2)));
                    //he4 might be missing edge!
                    if (hds.T(he4)) { he4=hds.next(he4); }
                    he4 = hds.next(he4);
                }
                double s4 = hds.itv(he4);
                
                //only one knot insertion for last vertex
                VV(_, nVV) = ((s(1)/2.0)*hds.V(_,hds.tip(he(2))) + (s(1)/2.0 + s(2) + s4)*hds.V(_,hds.tip(he(1)))) / (s(1) + s(2) + s4);
                Ktip( cv_size_t(he(1), hds.twin(hds.next(he(1)))) ) = nVV; nVV++;
            }
        }
    }
    
    
    
    
    
    
    // compute_TCC_max_stencils.m ------------------------------------------------
    
    struct M_Neighborhood
    {
        v_double n, p, h;
        double n4k, p2k;
        
        M_Neighborhood(size_t N): n(1,N), p(1,N), h(1,N) {};
    };
    
    struct I_Neighborhood
    {
        v_double o, n, p, m, q, s, t;
        
        I_Neighborhood(size_t N): o(1,N), n(1,N), p(1,N), m(1,N), q(1,N), s(1,2), t(1,N) {};
    };
    
    struct H_Neighborhood
    {
        v_size_t e, o, n, p, m, q;
        
        H_Neighborhood(size_t N): e(1,N), o(1,N), n(1,N), p(1,N), m(1,N), q(1,N) {};
    };
    
    struct T_Neighborhood
    {
        v_bool o, n, p, r;
        
        T_Neighborhood(size_t N): o(1,N), n(1,N), p(1,N), r(1,N) {};
    };
    
    void fix_borders_4(HDS &hds, H_Neighborhood &h, T_Neighborhood &T, I_Neighborhood &I);
    void fix_borders_5(HDS &hds, H_Neighborhood &h, T_Neighborhood &T, I_Neighborhood &I);
    
    bool compute_TCC_max_stencils_4(HDS &hds, uint64_t fI, std::v_double &Fs, std::v_double &Es, std::v_double &Vs);
    bool compute_TCC_max_stencils_5(HDS &hds, uint64_t fI, std::v_double &Fs, std::v_double &Es, std::v_double &Vs);
    
    bool compute_TCC_max_stencils(HDS &hds, uint64_t fI, std::v_double &Fs, std::v_double &Es, std::v_double &Vs)
    {
        if (hds.nFV[fI]==4)
        {
            return compute_TCC_max_stencils_4(hds, fI, Fs, Es, Vs);
        }
        else
        {
            return compute_TCC_max_stencils_5(hds, fI, Fs, Es, Vs);
        }
        return true;
    }
    
    void find(std::valarray<bool> v, std::valarray<size_t> &r)
    {
        std::vector<size_t> idx;
        for (size_t k=0; k<v.size(); k++)
        {
            if (v[k]) idx.push_back(k);
        }
        
        r.resize(idx.size());
        copy(idx.begin(), idx.end(), &r[0]);
    }
    
    struct TCirc // T- Circulator
    {
        HDS *hds;
        size_t c, n;
        TCirc(HDS *hhds, size_t cc, size_t nn): hds(hhds), c(cc), n(nn) {}
        bool operator==(TCirc tc2) { return (c==tc2.c && n==tc2.n); }
    };
    
    TCirc init_TC(HDS &hds, size_t he)
    {
        TCirc tc(&hds, he, MAX_SIZE_T);
        if (!hds.T(he)) tc.n = hds.twin[hds.next[he]];
        return tc;
    }
    
    TCirc next_TC(TCirc tc)
    {
        HDS &hds = *tc.hds;
        TCirc tcN(&hds, tc.n, MAX_SIZE_T);
        if (tc.n==MAX_SIZE_T)
        {
            tcN.n = hds.twin[hds.next[tc.c]];
        }
        else
        {
            if (!hds.T[tcN.c]) tcN.n = hds.twin[hds.next[tcN.c]];
        }
        return tcN;
    }
    
    TCirc prev_TC(TCirc tc)
    {
        HDS &hds = *tc.hds;
        TCirc tcP(&hds, MAX_SIZE_T, MAX_SIZE_T);
        if (tc.c==MAX_SIZE_T)
        {
            tcP.c = hds.prev[hds.twin[tc.n]];
        }
        else
        {
            tcP.c = hds.prev[hds.twin[tc.c]];
            if (hds.T[tcP.c]) tcP.c = MAX_SIZE_T;
            tcP.n = tc.c;
        }
        return tcP;
    }
    
    double get_itv(TCirc tc)
    {
        HDS &hds = *tc.hds;
        double itv = 0;
        if (tc.c==MAX_SIZE_T)
        {
            itv = hds.itv(hds.next(tc.n));
        }
        else
        {
            itv = hds.itv(tc.c);
        }
        return itv;
    }
    
    double sector_max(TCirc tc_start, TCirc tc_end) // start/end included in computation
    {
        TCirc tc = tc_start;
        double itv = 0;
        
        while (1)
        {
            itv = max(itv, get_itv(tc));
            if (tc == tc_end) break;
            tc = next_TC(tc);
        }
        return itv;
    }
    
    
    
    void update_0T_max_Mp_Mn(HDS &hds, M_Neighborhood &IM, I_Neighborhood &I, T_Neighborhood &T, H_Neighborhood &h)
    {
        if (all( hds.val(hds.tip(h.e))==4 )) return;
        
        for (size_t kHE1 = 0; kHE1<4; kHE1++)
        {
            size_t he = h.e(kHE1);
            size_t v  = hds.tip(he);
            if (hds.val(v)==4) continue;   // skip valence 4 vertices
            size_t kHE2 = (kHE1+1) % 4;
            
            TCirc tc1 = prev_TC( init_TC(hds, h.e(kHE1)));
            TCirc tc2 = next_TC( init_TC(hds, h.o(kHE2)));
            
            if (hds.val(v)==3) // valence 3 vertex ?
            {
                IM.p(kHE1) = max(IM.p(kHE1), I.o(kHE1));
                IM.n(kHE2) = max(IM.n(kHE2), I.o(kHE2));
                continue;
            }
            
            tc1 = prev_TC( tc1 );
            tc2 = next_TC( tc2 );
            double Msec = sector_max(tc2, tc1);
            
            IM.p(kHE1) = max(IM.p(kHE1), Msec);
            IM.n(kHE2) = max(IM.n(kHE2), Msec);
        }
    }
    
//    void print_nbh( I_Neighborhood &I, M_Neighborhood &IM, T_Neighborhood &T, H_Neighborhood &h)
//    {
//        cout<<"I.q "; print_v_double(I.q);
//        cout<<"I.p "; print_v_double(I.p);
//        cout<<"I.o "; print_v_double(I.o);
//        cout<<"I.n "; print_v_double(I.n);
//        cout<<"I.m "; print_v_double(I.m);
//        
//        cout<<"IM.h "; print_v_double(IM.h);
//        cout<<"IM.n "; print_v_double(IM.n);
//        cout<<"IM.n4k "<<IM.n4k<<endl;
//        cout<<"IM.p "; print_v_double(IM.p);
//        cout<<"IM.p2k "<<IM.p2k<<endl;
//        
//        cout<<"T.o "; print_v_bool(T.o);
//        cout<<"T.p "; print_v_bool(T.p);
//        cout<<"T.n "; print_v_bool(T.n);
//        cout<<"T.r "; print_v_bool(T.r);
//        
//    }
    
    
    
    bool compute_TCC_max_stencils_4(HDS &hds, uint64_t fI, std::v_double &Fs, std::v_double &Es, std::v_double &Vs)
    {
        I_Neighborhood I(4);
        M_Neighborhood IM(4);
        T_Neighborhood T(4);
        H_Neighborhood h(4);
        
        size_t nC = hds.V.nR();
        
        h.e[0] = hds.F2H[fI]; for (size_t k=1; k<4; k++) h.e[k] = h.e[0] + k;
        
        h.o = hds.twin(h.e);
        h.n = hds.next(h.o);
        h.p = hds.prev(h.o);
        
        h.m = hds.next(hds.twin(h.n));
        h.q = hds.prev(hds.twin(h.p));
        
        T.o = hds.T(h.o);
        T.p = hds.T(h.p);
        T.n = hds.T(h.n);
        T.r = hds.T(hds.prev(h.p));
        
        // replace missing edges with opposing edges, opposing edges of
        // cannot have a T-joint or extensions would intersect.
        h.m(T.o) = h.n(T.o); h.n(T.o) = h.p(T.o);
        h.q(T.p) = h.p(T.p); h.p(T.p) = h.n(T.p);
        
        // the ~Tp and ~To prevent access to "fixed" missing hn / hp edges
        std::v_bool Tq = hds.T(h.q)           & !T.p; h.q(Tq) = hds.prev(h.q(Tq));
        std::v_bool Tm = hds.T(hds.twin(h.n)) & !T.o; h.m(Tm) = hds.next(h.m(Tm));
        
        I.o = hds.itv(h.e);
        I.n = hds.itv(h.n);
        I.p = hds.itv(h.p);
        I.m = hds.itv(h.m);
        I.q = hds.itv(h.q);
        
        fix_borders_4(hds, h, T, I);
        
        I.t = I.n; I.t(T.n) = I.p(T.n);
        I.s[0] = I.o[0]; I.s[1] = I.o[1];
        
        //std::v_double P(nC, 4, 0); P.v = hds.V.cols( hds.tip.v[h.e] );
        
        std::v_double P = hds.V(_, hds.tip(h.e));
        
        // edge ownership: if halfedge index smaller than twin's and no T-joints
        // "on the other side" (missing edge or T-joint on edge), then it is ours.
        std::v_bool EO = h.e < h.o;
        if ( (!any(T.o | T.p | T.n | T.r)) && all( hds.val(hds.tip(h.e))==4 ) )
        {
            // optimization for uniform case without T-joints
            double aUitv[] = { IM.n[1], IM.n[3], I.o[0], I.o[2], IM.p[1], IM.p[3] }; v_double Uitv(aUitv, 6);
            double aVitv[] = { IM.n[0], IM.n[2], I.o[1], I.o[3], IM.p[0], IM.p[2] }; v_double Vitv(aVitv, 6);
            
            if ( all(Uitv==Uitv[0]) && all(Vitv==Vitv[0]) )
            {
                Fs = ( P(_, 0) + P(_, 1) + P(_, 2) + P(_, 3) )/4;
                v_double FsE = Fs / 4.0, FsV = Fs * Uitv[0] * Vitv[0];
                
                for (size_t k=0; k<4; k++)
                {
                    size_t kP = (k-1) & 3;
                    Es(_, k) += FsE;
                    Vs(_, k) += FsV;
                    
                    if (EO[k])
                    {
                        v_double M = (P(_, k) + P(_, kP)) / 2.0;
                        
                        Es(_, k)  += M / 2.0;
                        Vs(_, k)  += M * 2.0 * Uitv[0]*Vitv[0];
                        Vs(_, kP) += M * 2.0 * Uitv[0]*Vitv[0];
                    }
                }
                return true;
            }
        }
        else
        {
            EO = EO & !T.n & !T.r; // if Tn and Tr, edges have to be computed by the T-face.
        }
        
        
        // convenience variables
        double aU1[] = { I.n[1], I.o[0], I.p[3]}, aU2[] = {I.p[1], I.o[2], I.n[3]}, aU[] = { I.t[1], I.o[0], I.t[3] };
        double aV1[] = { I.p[0], I.o[1], I.n[2]}, aV2[] = {I.n[0], I.o[3], I.p[2]}, aV[] = { I.t[0], I.o[1], I.t[2] };
        
        v_double U1(aU1, 3), U2(aU2, 3), U(aU, 3), V1(aV1, 3), V2(aV2, 3), V(aV, 3);
        
        /*
         // semi-standard corner case
         v_bool Tsemi(T.r & cshift(T.n, 1, v_col));
         if (any(Tsemi))
         {
         for (size_t kSC1=0; kSC1<4; kSC1++)
         {
         if (Tsemi[kSC1])
         {
         size_t kSC2 = (kSC1 + 1) & 3, kSC3 = (kSC1 + 2) & 3, kSC4 = (kSC1 + 3) & 3;
         
         v_double pSC = (I.t[kSC1] * I.t[kSC2] / ( 4.0 * sum(U) * sum(V))) * P(_, kSC3);
         P(_, kSC1) = P(_, kSC1) - pSC;
         
         double WA = (I.o[kSC1] + I.n[kSC2] + I.q[kSC1]), WB = (I.o[kSC1] + I.p[kSC4] + I.m[kSC1]), WAB = 2.0*(WA+WB);
         v_double MSC = (WB / WAB) * pSC;
         
         Es(_, kSC1) -= MSC;
         Vs(_, kSC4) -= MSC * (I.o[kSC4] + I.n[kSC1]) * (I.p[kSC4] + I.m[kSC1]);
         
         WA = (I.o[kSC2] + I.n[kSC3] + I.q[kSC2]); WB = (I.o[kSC2] + I.p[kSC1] + I.m[kSC2]); WAB = 2.0*(WA+WB);
         MSC = (WA / WAB) * pSC;
         Es(_, kSC2) -= MSC;
         Vs(_, kSC2) -= MSC * (I.o[kSC3] + I.p[kSC2]) * (I.n[kSC3] + I.q[kSC2]);
         }
         }
         }
         */
        
        IM.n = I.n;
        IM.p = I.p;
        
        update_0T_max_Mp_Mn(hds, IM, I, T, h);
        
        // Face components of stencils
        double Ao[] = {2*(IM.n[1] + I.o[0] + IM.p[3]), 2*(IM.p[0] + I.o[1] + IM.n[2]), 2*(IM.p[1] + I.o[2] + IM.n[3]), 2*(IM.p[2] + I.o[3] + IM.n[0])};
        
        double W0 = (2 * IM.n[2] + I.o[1]) * (2 * IM.p[3] + I.o[0]) / (Ao[0]*Ao[1]);
        double W1 = (2 * IM.p[0] + I.o[1]) * (2 * IM.n[3] + I.o[2]) / (Ao[1]*Ao[2]);
        double W2 = (2 * IM.n[0] + I.o[3]) * (2 * IM.p[1] + I.o[2]) / (Ao[2]*Ao[3]);
        double W3 = (2 * IM.p[2] + I.o[3]) * (2 * IM.n[1] + I.o[0]) / (Ao[3]*Ao[0]);
        double sW = W0+W1+W2+W3; W0/=sW; W1/=sW; W2/=sW; W3/=sW;
        
        // face stencil
        Fs = W0 * P(_,0) + W1 * P(_,1) + W2 * P(_,2) + W3 * P(_,3);
        
        
        // edge stencil face component
        // see pics/EdgeStencilConditionals.png for details on conditionals
        // weights use mid-face knot vectors, hence U,V instead of U1,V1,U2,V2
        
        Es(_, 0) += ( (W0 + T.r[0]*W1) * P(_,0) + !T.r[0]*W1 * P(_,1) + !T.n[0]*W2 * P(_,2) + (T.n[0]*W2 + W3) * P(_,3) ) * I.t[0]/(2.0*(I.t[0]+I.s[1]));
        Es(_, 1) += ( (W0 + T.n[1]*W3) * P(_,0) + (W1 + T.r[1]*W2) * P(_,1) + !T.r[1]*W2 * P(_,2) + !T.n[1]*W3 * P(_,3) ) * I.t[1]/(2.0*(I.t[1]+I.s[0]));
        Es(_, 2) += ( !T.n[2]*W0 * P(_,0) + (T.n[2]*W0 + W1) * P(_,1) + (W2 + T.r[2]*W3) * P(_,2) + !T.r[2]*W3 * P(_,3) ) * I.t[2]/(2.0*(I.t[2]+I.s[1]));
        Es(_, 3) += ( !T.r[3]*W0 * P(_,0) + !T.n[3]*W1 * P(_,1) + (T.n[3]*W1 + W2) * P(_,2) + (T.r[3]*W0 + W3) * P(_,3) ) * I.t[3]/(2.0*(I.t[3]+I.s[0]));
        
        // vertex stencil face component
        // see pics/VertexStencilConditionals.png for details on conditionals
        Vs(_,0) += ( W0 * P(_,0) + (W1 + T.r[1]*W2) * P(_,1) + !T.r[1]*!T.n[0]*W2 * P(_,2) + (T.n[0]*W2 + W3) * P(_,3) ) * I.p[0]*I.n[1];
        Vs(_,1) += ( (W0 + T.n[1]*W3) * P(_,0) + W1 * P(_,1) + (W2 + T.r[2]*W3) * P(_,2) + !T.n[1]*!T.r[2]*W3 * P(_,3) ) * I.p[1]*I.n[2];
        Vs(_,2) += ( !T.r[3]*!T.n[2]*W0 * P(_,0) + (W1 + T.n[2]*W0) * P(_,1) + W2 * P(_,2) + (T.r[3]*W0 + W3) * P(_,3) ) * I.p[2]*I.n[3];
        Vs(_,3) += ( (W0 + T.r[0]*W1) * P(_,0) + !T.r[0]*!T.n[3]*W1 * P(_,1) + (W2 + T.n[3]*W1) * P(_,2) + W3 * P(_,3) ) * I.p[3]*I.n[0];
        
        // Edge components of stencils, if we have edge ownership
        double WA, WB, WAB;
        v_double M(0.0, nC);
        
        if (EO[0])
        {
            WA = (I.o[0] + 2 * IM.n[1]); WB = (I.o[0] + 2 * IM.p[3]); WAB = 2*(WA+WB);
            
            M = (WA/WAB)*P(_, 3) + (WB/WAB)*P(_, 0);
            Es(_, 0) += M;
            Vs(_, 3) += M * (I.o[3] + I.n[0]) * (I.p[3] + I.m[0]);
            Vs(_, 0) += M * (I.o[1] + I.p[0]) * (I.n[1] + I.q[0]);
        }
        
        if (EO[1])
        {
            WA = (I.o[1] + 2 * IM.n[2]); WB = (I.o[1] + 2 * IM.p[0]); WAB = 2*(WA+WB);
            
            M = (WA/WAB)*P(_, 0) + (WB/WAB)*P(_, 1);
            Es(_, 1) += M;
            Vs(_, 0) += M * (I.o[0] + I.n[1]) * (I.p[0] + I.m[1]);
            Vs(_, 1) += M * (I.o[2] + I.p[1]) * (I.n[2] + I.q[1]);
        }
        
        if (EO[2])
        {
            WA = (I.o[2] + 2 * IM.n[3]); WB = (I.o[2] + 2 * IM.p[1]); WAB = 2*(WA+WB);
            
            M = (WA/WAB)*P(_, 1) + (WB/WAB)*P(_, 2);
            Es(_, 2) += M;
            Vs(_, 1) += M * (I.o[1] + I.n[2]) * (I.p[1] + I.m[2]);
            Vs(_, 2) += M * (I.o[3] + I.p[2]) * (I.n[3] + I.q[2]);
        }
        
        if (EO[3])
        {
            WA = (I.o[3] + 2 * IM.n[0]); WB = (I.o[3] + 2 * IM.p[2]); WAB = 2*(WA+WB);
            
            M = (WA/WAB)*P(_, 2) + (WB/WAB)*P(_, 3);
            Es(_, 3) += M;
            Vs(_, 2) += M * (I.o[2] + I.n[3]) * (I.p[2] + I.m[3]);
            Vs(_, 3) += M * (I.o[0] + I.p[3]) * (I.n[0] + I.q[3]);
        }
        
        return true;
    }
    
    void fix_borders_4(HDS &hds, H_Neighborhood &h, T_Neighborhood &T, I_Neighborhood &I)
    {
        //deal with border edges
        
        //Important:
        // we assume that border edges have no T-joints set,
        // so previous code doesn't wreak havoc on supposedly "missing"
        // halfedges
        
        size_t nIH = hds.H2F.size();
        
        v_bool Bo = (h.o >= nIH);
        if (any(Bo))
        {
            v_size_t BI = find(Bo);
            
            v_size_t BIP = (BI - size_t(1)) & 3;
            v_size_t BIN = (BI + size_t(1)) & 3;
            
            I.n(Bo) = I.o(BIP); T.n(Bo) = false;
            I.p(Bo) = I.o(BIN); T.p(Bo) = false; T.r(Bo) = false;
            I.q(Bo) = I.n(BIN);
            I.m(Bo) = I.p(BIP);
        }
        
        v_bool Bn = !Bo & hds.twin(h.n)>=nIH;
        if (any(Bn))
        {
            v_size_t BI = find(Bn);
            v_size_t BIP = (BI - 1) & 3;
            I.m(Bn) = I.p(BIP);
        }
        
        v_bool Bp = !Bo & hds.twin(h.p)>=nIH;
        if (any(Bp))
        {
            v_size_t BI = find(Bp);
            v_size_t BIN = (BI + 1) & 3;
            I.q(Bp) = I.n(BIN);
        }
        
    }
    
    
    void update_1T_max_Mp_Mn(HDS &hds, M_Neighborhood &IM, I_Neighborhood &I, T_Neighborhood &T, H_Neighborhood &h)
    {
        if (all(hds.val(hds.tip(h.e))==4)) return;
        
        for (size_t kHE1 = 1; kHE1<5; kHE1++)
        {
            size_t he = h.e(kHE1);
            size_t v  = hds.tip(he);
            if (hds.val(v)==4) continue; // skip valence 4 vertices
            size_t kHE2 = (kHE1+1) % 5;
            
            TCirc tc1 = prev_TC( init_TC(hds, h.e(kHE1)));
            TCirc tc2 = next_TC( init_TC(hds, h.o(kHE2)));
            
            if (hds.val(v)==3)
            {
                IM.p(kHE1) = max(IM.p(kHE1), I.o(kHE1));
                IM.n(kHE2) = max(IM.n(kHE2), I.o(kHE2));
                
                if (kHE1==2) IM.p2k = max(IM.p2k, I.o(kHE1));
                if (kHE2==4) IM.n4k = max(IM.n4k, I.o(kHE2));
                continue;
            }
            
            tc1 = prev_TC( tc1 );
            tc2 = next_TC( tc2 );
            double Msec = sector_max(tc2, tc1);
            
            IM.p(kHE1) = max( IM.p(kHE1), Msec );
            IM.n(kHE2) = max( IM.n(kHE2), Msec );
            
            // update (possibly) half-knot intervals
            if (kHE1==2) IM.p2k = max(IM.p2k, Msec);
            if (kHE2==4) IM.n4k = max(IM.n4k, Msec);
        }
        
        // T-vertex
        size_t v = hds.tip(h.e(0));
        if (hds.val(v)!=4)  // extraordinary T-joint vertex ?
        {
            TCirc tc1 = init_TC(hds, h.e(0)), tc1p = prev_TC(tc1), tc1pp = prev_TC(tc1p);
            TCirc tc2 = init_TC(hds, h.o(1)), tc2n = next_TC(tc2), tc2nn = next_TC(tc2n);
            
            IM.h(0) = sector_max(tc2, tc1pp);
            IM.h(1) = sector_max(tc2nn, tc1);
            IM.h(2) = sector_max(tc2n, tc1p);
        }
        
    }
    
    
    bool compute_TCC_max_stencils_5(HDS &hds, uint64_t fI, std::v_double &Fs, std::v_double &Es, std::v_double &Vs)
    {
        I_Neighborhood I(5);
        M_Neighborhood IM(5);
        T_Neighborhood T(5);
        H_Neighborhood h(5);
        
        size_t nC = hds.V.nR();
        
        h.e[0] = hds.F2H[fI];
        size_t csh = 0;
        for (size_t k=1; k<5; k++)
        {
            h.e[k] = h.e[0] + k; if (hds.T.v[h.e[0] + k]) csh = k;
        }
        h.e = cshift(h.e, csh, v_col);
        Es  = cshift(Es, csh, v_col);
        Vs  = cshift(Vs, csh, v_col);
        
        h.o = hds.twin(h.e);
        h.n = hds.next(h.o);
        h.p = hds.prev(h.o);
        
        h.m = hds.next(hds.twin(h.n));
        h.q = hds.prev(hds.twin(h.p));
        
        T.o = hds.T(h.o);
        T.p = hds.T(h.p);
        T.n = hds.T(h.n);
        T.r = hds.T(hds.prev(h.p));
        
        // replace missing edges with opposing edges, opposing edges of
        // cannot have a T-joint or extensions would intersect.
        h.m(T.o) = h.n(T.o); h.n(T.o) = h.p(T.o);
        h.q(T.p) = h.p(T.p); h.p(T.p) = h.n(T.p);
        
        // the ~Tp and ~To prevent access to "fixed" missing hn / hp edges
        v_bool Tq = hds.T(h.q) & !T.p;           h.q(Tq) = hds.prev(h.q(Tq));
        v_bool Tm = hds.T(hds.twin(h.n)) & !T.o; h.m(Tm) = hds.next(h.m(Tm));
        
        // compute knot intervals
        I.o = hds.itv(h.e);
        I.n = hds.itv(h.n);
        I.p = hds.itv(h.p);
        I.m = hds.itv(h.m);
        I.q = hds.itv(h.q);
        
        fix_borders_5(hds, h, T, I);
        
        IM.p = I.p;
        IM.n = I.n;
        IM.h[0] = I.o[1]; IM.h[1] = I.o[0]; IM.h[2] = I.p[0];
        IM.n4k = T.n[4]*IM.n[4] + !T.n[4]*IM.p[4];
        IM.p2k = T.r[2]*IM.p[2] + !T.r[2]*IM.n[2];
        
        update_1T_max_Mp_Mn(hds, IM, I, T, h);
        
        I.t = I.n; I.t(T.n) = I.p(T.n);
        I.s[0] = I.o[3]; I.s[1] = I.o[2];
        
        
        // convenience variables
        double aU1[]  = { I.n[2], I.o[3], I.p[4] };
        double aU2[]  = { I.p[2], I.o[3], I.n[4] };
        double aV[]   = { I.p[1], I.o[2], I.n[3] }; // ToDO!!! Ip(1) is not a good value at extraordinary T-vertices
        v_double U1(aU1,3), U2(aU2,3), V(aV,3);
        
        v_double P( nC, 5, 0 ), PK(nC, 3, 0);
        P = hds.V(_, hds.tip(h.e));
        
        PK(_, 0) = hds.VK(_, hds.tipK( h.e[2] ) );
        PK(_, 1) = hds.VFK(_, fI );
        PK(_, 2) = hds.VK(_, hds.tipK( h.e[3] ) );
        
        // for edges 2 and 4 ownership is ours if no neighboring T-face or if our edge index wins,
        // if neighbor face nonexistent h.o[x]>fH[x] by construction
        std::valarray<bool> EO(false, 5);
        EO[0] = h.e[0]<h.o[0];
        EO[1] = h.e[1]<h.o[1];
        EO[2] = (h.e[2]<h.o[2]) || (hds.nFV.v[hds.H2F.v[h.o[2]]]==4);
        EO[3] = h.e[3]<h.o[3];
        EO[4] = (h.e[4]<h.o[4]) || (hds.nFV.v[hds.H2F.v[h.o[4]]]==4);
        
        // Face components of stencils
        double Ao[] = { 2.0*(IM.p[4]+I.o[0]+I.o[1]+IM.n[2]), 2.0*(IM.p[1]+I.o[2]+IM.n[3]), 2.0*(IM.p[2]+I.o[3]+IM.n[4]), 2.0*(IM.p[3]+I.o[4]+IM.n[0])};
        double At = 2*(IM.h[2] + I.s[1] + I.t[3]);
        
        double W0 = (I.s[1] + 2*I.t[3]) / At;
        double W2 = (I.o[2] + 2*IM.p[1]) * (I.o[3] + 2*IM.n[4]) / (Ao[1]*Ao[2]);
        double W3 = (I.o[4] + 2*IM.n[0]) * (I.o[3] + 2*IM.p[2]) / (Ao[2]*Ao[3]);
        double sW = W0+W2+W3; W0/=sW; W2/=sW; W3/=sW;
        
        // face stencil
        Fs = W0 * P(_,0) + W2 * P(_,2) + W3 * P(_,3);
        
        // knot-inserted if no t-joints on upper side, else take the upper knot intervals
        
        
        double AoR[] = { 2.0*(IM.p[4]+I.o[0]+IM.h[0]), 2.0*(IM.h[2]+I.s[1]+I.t[3]), 2.0*(I.o[3]+IM.n4k), 2.0*(IM.p[3]+I.o[4]+IM.n[0])};
        double WR0 = (I.o[0] + 2*IM.p[4]) * (I.s[1] + 2*I.t[3] ) / (AoR[0]*AoR[1]);
        double WR1 = (I.o[0] + 2*IM.n4k ) * (I.s[1] + 2*IM.h[2]) / (AoR[1]*AoR[2]);
        double WR2 =           (3*I.o[0]) * (I.o[4] + 2*IM.n[0]) / (AoR[2]*AoR[3]);
        double WR3 = (I.o[0] + 2*IM.h[0]) * (I.o[4] + 2*IM.p[3]) / (AoR[3]*AoR[0]);
        double sWR = WR0+WR1+WR2+WR3; WR0/=sWR; WR1/=sWR; WR2/=sWR; WR3/=sWR;
        
        double AoL[] = { 2.0*(IM.h[1]+I.o[1]+IM.n[2]), 2.0*(IM.p[1]+I.o[2]+IM.n[3]), 2.0*(IM.p2k+I.o[3]), 2.0*(I.t[3]+I.s[1]+IM.h[2])};
        double WL0 = (I.o[1] + 2*IM.h[1]) * (I.o[2] + 2*IM.n[3]) / (AoL[0]*AoL[1]);
        double WL1 =           (3*I.o[1]) * (I.o[2] + 2*IM.p[1]) / (AoL[1]*AoL[2]);
        double WL2 = (I.o[1] + 2*IM.p2k ) * (I.s[1] + 2*IM.h[2]) / (AoL[2]*AoL[3]);
        double WL3 = (I.o[1] + 2*IM.n[2]) * (I.s[1] + 2*I.t[3] ) / (AoL[3]*AoL[0]);
        double sWL = WL0+WL1+WL2+WL3; WL0/=sWL; WL1/=sWL; WL2/=sWL; WL3/=sWL;
        
        // for split faces, VR2U -> VR1U, VR3U -> VR4U, VL2U -> VL1U, VL3U -> VL4U
        //v_double FH1 = WR0 * P(_,0) + WR1 * PK(_,1) + WR2 * PK(_,2) + WR3 * P(_,4);
        //v_double FH2 = WL0 * P(_,1) + WL1 * PK(_,0) + WL2 * PK(_,1) + WL3 * P(_,0);
        
        
        // edge stencil face componentp
        // see pics/EdgeStencilConditImals.png for details on conditImals
        // weights use mid-face knot vectors, hence U,V instead of U1,V1,U2,V2
        // for split faces, VR2U -> VR1U, VR3U -> VR4U
        Es(_,0) += ((WR0 + T.r[0]*WR1) * P(_,0) + !T.r[0]*WR1 * PK(_,1) + !T.n[0]*WR2 * PK(_,2) + (WR3 + T.n[0]*WR2) * P(_,4)) * I.t[0] / (2*(I.o[4]+I.t[0]));
        Es(_,1) += ((WL0 + T.r[1]*WL1) * P(_,1) + !T.r[1]*WL1 * PK(_,0) + !T.n[1]*WL2 * PK(_,1) + (WL3 + T.n[1]*WL2) * P(_,0)) * I.t[1] / (2*(I.o[2]+I.t[1]));
        Es(_,2) +=            (!T.n[2] * W0 * P(_,0) + T.n[2] * W0 * P(_,1) + (W2 + T.r[2]*W3) * P(_,2) + !T.r[2]*W3 * P(_,3)) * I.t[2] / (2*(I.o[3]+I.t[2]));
        Es(_,3) +=                                                                   (W0 * P(_,0) + W2 * P(_,2) + W3 * P(_,3)) * I.t[3] / (2*(I.o[2]+I.t[3]));
        Es(_,4) +=            (!T.r[4] * W0 * P(_,0) + !T.n[4]*W2 * P(_,2) + (W3 + T.n[4]*W2) * P(_,3) + T.r[4] * W0 * P(_,4)) * I.t[4] / (2*(I.o[3]+I.t[4]));
        
        // vertex stencil face component
        // see pics/VertexStencilConditImals.png for details on conditioals
        Vs(_,0) +=                        (WR0 * P(_,0) + WR1 * PK(_,1) + !T.n[0]*WR2 * PK(_,2) + (WR3 + T.n[0]*WR2) * P(_,4)) * (I.o[1]*I.p[0]) +
        ((WL0 + T.r[1]*WL1) * P(_,1) + !T.r[1]*WL1 * PK(_,0) + WL2 * PK(_,1) + WL3 * P(_,0)) * (I.o[0]*I.n[1]);
        Vs(_,1) += (WL0 * P(_,1) + (WL1 + T.r[2]*WL2) * PK(_,0) + !T.r[2]*!T.n[1]*WL2 * PK(_,1) + (WL3 + T.n[1]*WL2) * P(_,0)) * (I.p[1]*I.n[2]);
        Vs(_,2) +=                                  (!T.n[2] * W0 * P(_,0) + T.n[2] * W0 * P(_,1) + W2 * P(_,2) + W3 * P(_,3)) * (I.p[2]*I.n[3]);
        Vs(_,3) +=                                  (!T.r[4] * W0 * P(_,0) + W2 * P(_,2) + W3 * P(_,3) + T.r[4] * W0 * P(_,4)) * (I.p[3]*I.n[4]);
        Vs(_,4) += ((WR0 + T.r[0]*WR1) * P(_,0) + !T.r[0]*!T.n[4]*WR1 * PK(_,1) + (WR2 + T.n[4]*WR1) * PK(_,2) + WR3 * P(_,4)) * (I.p[4]*I.n[0]);
        
        // Edge components of stencils, if we have edge ownership
        // (T-faces have automatic edge ownership for all edges except edge 4)
        
        // this is the missing edge for the T-joint
        v_double M(0.0, nC);
        
        double WA = (2*I.t[3] + I.s[1]); double WB = (I.s[1] + 2*IM.h[2]); double WAB = 2*(WA+WB);
        M = (WA/WAB)*P(_, 0) + (WB/WAB)*PK(_, 1);
        Vs(_, 0) += M * I.o[3] * (I.p[0] + I.n[1]);
        
        if (EO[0])
        {
            WA = (I.o[0] + 2*IM.h[0]); WB = (I.o[0] + 2*IM.p[4]); WAB = 2*(WA+WB);
            M = (WA/WAB)*P(_, 4) + (WB/WAB)*P(_, 0);
            Es(_, 0) += M;
            Vs(_, 4) += M * (I.o[4] + I.n[0]) * (I.p[4] + I.m[0]);
            Vs(_, 0) += M * (I.o[4] + I.p[0]) * (I.o[1] + I.q[0]);
        }
        
        if (EO[1])
        {
            WA = (I.o[1] + 2*IM.n[2]); WB = (I.o[1] + 2*IM.h[1]); WAB = 2*(WA+WB);
            M = (WA/WAB)*P(_, 0) + (WB/WAB)*P(_, 1);
            Es(_, 1) += M;
            Vs(_, 0) += M * (I.o[2] + I.n[1]) * (I.o[0] + I.m[1]);
            Vs(_, 1) += M * (I.o[2] + I.p[1]) * (I.n[2] + I.q[1]);
        }
        
        if (EO[2] || !T.n[2])
        {
            WA = (I.o[2] + 2*IM.n[3]); WB = (I.o[2] + 2*IM.p[1]); WAB = 2*(WA+WB);
            M = (WA/WAB)*P(_, 1) + (WB/WAB)*PK(_, 0);  // PK instead of P
            Vs(_, 1) += M * (I.o[1] + I.n[2]) * (I.p[1] + I.m[2]);
            if (EO[2])
            {
                M = (WA/WAB)*P(_, 1) + (WB/WAB)*P(_, 2);
                Es(_, 2) = Es(_, 2) + M;
                if (!T.r[2])
                {
                    Vs(_, 2) += M * (I.o[3] + I.p[2]) * (I.n[3] + I.q[2]);
                }
            }
        }
        
        if (EO[3])
        {
            WA = (I.o[3] + 2*IM.n[4]); WB = (I.o[3] + 2*IM.p[2]); WAB = 2*(WA+WB);
            
            M = (WA/WAB)*P(_, 2) + (WB/WAB)*P(_, 3);
            Es(_, 3) += M;
            Vs(_, 2) += M * (I.o[2] + I.n[3]) * (I.p[2] + I.m[3]);
            Vs(_, 3) += M * (I.o[4] + I.p[3]) * (I.n[4] + I.q[3]);
        }
        
        
        if (EO[4] || !T.r[4])
        {
            WA = (I.o[4] + 2*IM.n[0]); WB = (I.o[4] + 2*IM.p[3]); WAB = 2*(WA+WB);
            M = (WA/WAB)*PK(_, 2) + (WB/WAB)*P(_, 4);  //PK instead of P
            Vs(_, 4) += M * (I.o[0] + I.p[4]) * (I.n[0] + I.q[4]);
            if (EO[4])
            {
                M = (WA/WAB)*P(_, 3) + (WB/WAB)*P(_, 4);
                Es(_, 4) += M;
                if (!T.n[4])
                {
                    Vs(_, 3) += M * (I.o[3] + I.n[4]) * (I.p[3] + I.m[4]);
                }
            }
        }
        
        Es = cshift(Es, -csh, v_col);
        Vs = cshift(Vs, -csh, v_col);
        
        return true;
    }
    
    void fix_borders_5(HDS &hds, H_Neighborhood &h, T_Neighborhood &T, I_Neighborhood &I)
    {
        //deal with border edges
        
        //Important:
        // we assume that border edges have no T-joints set,
        // so previous code doesn't wreak havoc on supposedly "missing"
        // halfedges
        
        size_t nIH = hds.H2F.size();
        
        v_bool Bo = (h.o >= nIH);
        if (any(Bo))
        {
            v_size_t BI = find(Bo);
            v_size_t BIP = (BI + size_t(4)) % 5;
            v_size_t BIN = (BI + size_t(1)) % 5;
            
            I.n(Bo) = I.o(BIP); T.n(Bo) = BI==2;
            I.p(Bo) = I.o(BIN); T.r(Bo) = BI==4; T.p(Bo) = false;
            I.q(Bo) = I.n(BIN);
            I.m(Bo) = I.p(BIP);
        }
        
        v_bool Bn = !Bo & hds.twin(h.n)>=nIH;
        if (any(Bn))
        {
            v_size_t BI = find(Bn);
            v_size_t BIP = (BI + size_t(4)) % 5;
            I.m(Bn) = I.p(BIP);
        }
        
        v_bool Bp = !Bo & hds.twin(h.p)>=nIH;
        if (any(Bp))
        {
            v_size_t BI = find(Bp);
            v_size_t BIN = (BI + size_t(1)) % 5;
            I.q(Bp) = I.n(BIN);
        }
    }
    
}
