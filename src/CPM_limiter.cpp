/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include<vector>
#include<algorithm>

/*--------------------------------------------------------------------
Grid-limited step

CPM is not DEM: parcels never see each other, they see the background mesh.
The packing limit is a kinematic constraint on the mesh instead of a stiff stress:

  - a parcel moves at most a fraction of a cell per sub-step (CFL S 14)
  - it enters a new cell only if that cell has free volume:
        A = V_max - V_occupied + V_leaving,   V_max = max(theta_max V, theta_0 V + V_parcel)
    i.e. at least one parcel more than the packing of the bed, so that a sheared bed can move
    (occupancy = parcel volume per cell, nearest cell)
  - first come, first served: the requests for a cell are ordered by their arrival time
    (fraction of the step at which the parcel reaches the cell), the cell accepts until
    its free volume is used up; the result is a cut-off arrival time per cell
  - a rejected parcel stops at the face, its velocity into the full cell is removed
  - solid bodies (solid level set < 0) have no free volume

Passes: an accepted move is final; its volume leaves the source cell and fills the
target cell in the following passes, the pending requests compete for the remaining
free volume (monotone, the occupancy never exceeds theta_max).

Parallel: for cells with requests from other ranks, the requests are binned by arrival
time into NB bins and summed into the owner of the cell with start4a_sum, the owner
accepts whole bins that fit and returns the cut-off by a ghost cell exchange. Cells with
local requests only are served exactly in the order of arrival. No parcel lists cross
the ranks, the result does not depend on the order of the parcels.
--------------------------------------------------------------------*/

void CPM::limiter(lexer *p, fdm *a, ghostcell *pgc, double *X0, double *Y0, double *Z0,
                  double *X1, double *Y1, double *Z1, double *PU, double *PV, double *PW)
{
    const int NB = 4;
    field *hist[4] = {&Lh0,&Lh1,&Lh2,&Lh3};
    
    const double vpar = P.ParcelFactor*Vp;
    
    int i0,j0,k0,i1,j1,k1;
    double ta,t;
    
    std::vector<int> ci(P.index),cj(P.index),ck(P.index);
    std::vector<int> si(P.index),sj(P.index),sk(P.index);
    std::vector<double> tarr(P.index);
    
    // 0 stays in its cell, 1 pending request (rejected if still pending at the end), 2 accepted, 4 accepted without check
    std::vector<char> state(P.index,0);
    
    auto local = [&](int ii, int jj, int kk) {return ii>=0 && ii<p->knox && jj>=0 && jj<p->knoy && kk>=0 && kk<p->knoz;};
    auto bin = [&](double tt) {return MIN(NB-1, int(tt*NB));};
    
    // occupancy at the start of the step
    for(i=-1;i<p->knox+1;++i)
    for(j=-1;j<p->knoy+1;++j)
    for(k=-1;k<p->knoz+1;++k)
    Locc(i,j,k) = 0.0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        i0 = p->posc_i(X0[n]);
        j0 = p->j_dir==1 ? p->posc_j(Y0[n]) : 0;
        k0 = p->posc_k(Z0[n]);
        
        if(local(i0,j0,k0))
        Locc(i0,j0,k0) += vpar;
    }
    
    // requests
    int nreq=0;
    int nclip=0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        i0 = p->posc_i(X0[n]);
        j0 = p->j_dir==1 ? p->posc_j(Y0[n]) : 0;
        k0 = p->posc_k(Z0[n]);
        
        // a move of more than one cell is shortened to stay within the neighbour cell
        {
            double hx = p->DXN[i0+marge], hz = p->DZN[k0+marge];
            double hy = p->j_dir==1 ? p->DYN[j0+marge] : 1.0e20;
            double sc = 1.0;
            
            if(fabs(X1[n]-X0[n])>0.9*hx) sc = MIN(sc, 0.9*hx/fabs(X1[n]-X0[n]));
            if(fabs(Y1[n]-Y0[n])>0.9*hy) sc = MIN(sc, 0.9*hy/fabs(Y1[n]-Y0[n]));
            if(fabs(Z1[n]-Z0[n])>0.9*hz) sc = MIN(sc, 0.9*hz/fabs(Z1[n]-Z0[n]));
            
            if(sc<1.0)
            {
                X1[n] = X0[n] + sc*(X1[n]-X0[n]);
                Y1[n] = Y0[n] + sc*(Y1[n]-Y0[n]);
                Z1[n] = Z0[n] + sc*(Z1[n]-Z0[n]);
                ++nclip;
            }
        }
        
        i1 = p->posc_i(X1[n]);
        j1 = p->j_dir==1 ? p->posc_j(Y1[n]) : 0;
        k1 = p->posc_k(Z1[n]);
        
        si[n]=i0; sj[n]=j0; sk[n]=k0;
        ci[n]=i1; cj[n]=j1; ck[n]=k1;
        
        if(i0==i1 && j0==j1 && k0==k1)
        continue;
        
        // arrival time: the last face crossed
        ta = 0.0;
        
        if(i1!=i0 && fabs(X1[n]-X0[n])>1.0e-20)
        {
            t = i1>i0 ? (p->XN[i0+1+marge]-X0[n])/(X1[n]-X0[n]) : (p->XN[i0+marge]-X0[n])/(X1[n]-X0[n]);
            ta = MAX(ta,t);
        }
        
        if(j1!=j0 && fabs(Y1[n]-Y0[n])>1.0e-20)
        {
            t = j1>j0 ? (p->YN[j0+1+marge]-Y0[n])/(Y1[n]-Y0[n]) : (p->YN[j0+marge]-Y0[n])/(Y1[n]-Y0[n]);
            ta = MAX(ta,t);
        }
        
        if(k1!=k0 && fabs(Z1[n]-Z0[n])>1.0e-20)
        {
            t = k1>k0 ? (p->ZN[k0+1+marge]-Z0[n])/(Z1[n]-Z0[n]) : (p->ZN[k0+marge]-Z0[n])/(Z1[n]-Z0[n]);
            ta = MAX(ta,t);
        }
        
        tarr[n] = MAX(0.0,MIN(ta,1.0));
        
        // a move of more than one cell or across a processor corner is not limited
        int out = (i1<0||i1>=p->knox) + (j1<0||j1>=p->knoy) + (k1<0||k1>=p->knoz);
        
        if(std::abs(i1-i0)>1 || std::abs(j1-j0)>1 || std::abs(k1-k0)>1 || out>1)
        state[n]=4;
        
        else
        {
            state[n]=1;
            ++nreq;
        }
    }
    
    nreq = pgc->globalisum(nreq);
    nclip_step += nclip;
    
    if(nreq==0)
    return;
    
    std::vector<int> order;
    order.reserve(P.index);
    
    // monotone passes: an accepted move is final, its volume leaves the source cell
    // and fills the target cell in the following passes; the pending requests compete
    // for the remaining free volume
    const int NPASS = 4;
    
    for(int pass=0; pass<NPASS; ++pass)
    {
        for(int b=0;b<NB;++b)
        for(i=-1;i<p->knox+1;++i)
        for(j=-1;j<p->knoy+1;++j)
        for(k=-1;k<p->knoz+1;++k)
        (*hist[b])(i,j,k) = 0.0;
        
        for(i=-1;i<p->knox+1;++i)
        for(j=-1;j<p->knoy+1;++j)
        for(k=-1;k<p->knoz+1;++k)
        {
        Lout(i,j,k) = 0.0;
        Lin(i,j,k) = 0.0;
        }
        
        for(n=0;n<P.index;++n)
        if(P.Flag[n]==ACTIVE)
        {
            if(state[n]==2 || state[n]==4)
            {
                if(local(si[n],sj[n],sk[n]))
                Lout(si[n],sj[n],sk[n]) += vpar;
                
                Lin(ci[n],cj[n],ck[n]) += vpar;
            }
            
            if(state[n]==1)
            (*hist[bin(tarr[n])])(ci[n],cj[n],ck[n]) += 1.0;
        }
        
        BASELOOP
        Lloc(i,j,k) = (*hist[0])(i,j,k) + (*hist[1])(i,j,k) + (*hist[2])(i,j,k) + (*hist[3])(i,j,k);
        
        for(int b=0;b<NB;++b)
        pgc->start4a_sum(p,*hist[b],1);
        
        pgc->start4a_sum(p,Lin,1);
        
        // free volume and cut-off per cell
        //   only local requests: exact first come first served below (Ltc > NB)
        //   requests from other ranks: whole arrival time bins that fit (Ltc = last bin, -1 none)
        BASELOOP
        {
            double V = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
            // capacity: theta_max, at least the bed packing plus one parcel (shear of the bed needs room)
            double cap = MAX(theta_max*V, theta_0*V + vpar*(1.0+1.0e-6));
            double A = cap - Locc(i,j,k) + Lout(i,j,k) - Lin(i,j,k);
            
            if(a->solid(i,j,k)<0.0)
            A = -1.0;
            
            LA(i,j,k) = A;
            
            double cnt = (*hist[0])(i,j,k) + (*hist[1])(i,j,k) + (*hist[2])(i,j,k) + (*hist[3])(i,j,k);
            
            if(cnt>Lloc(i,j,k)+0.5)
            {
                double slots = floor(A/vpar + 1.0e-9);
                double cum = 0.0;
                double tc = -1.0;
                
                for(int b=0;b<NB;++b)
                {
                    double hb = (*hist[b])(i,j,k);
                    
                    if(cum+hb<=slots+1.0e-9)
                    {
                        cum+=hb;
                        tc = double(b);
                    }
                    else
                    break;
                }
                
                Ltc(i,j,k) = tc;
            }
            else
            Ltc(i,j,k) = 2.0*NB;
        }
        
        pgc->start4a(p,Ltc,1);
        
        order.clear();
        int nacc=0;
        
        for(n=0;n<P.index;++n)
        if(P.Flag[n]==ACTIVE && state[n]==1)
        {
            double tc = Ltc(ci[n],cj[n],ck[n]);
            
            if(tc>double(NB) && local(ci[n],cj[n],ck[n]))
            order.push_back(n);
            
            else if(double(bin(tarr[n]))<=tc+1.0e-9)
            {
                state[n]=2;
                ++nacc;
            }
        }
        
        // exact first come first served for cells with local requests only, ties by index;
        // a worklist resolves chains: a parcel leaving a cell frees volume for the requests into it
        std::sort(order.begin(), order.end(), [&](int n1, int n2)
        {
            if(ci[n1]!=ci[n2]) return ci[n1]<ci[n2];
            if(cj[n1]!=cj[n2]) return cj[n1]<cj[n2];
            if(ck[n1]!=ck[n2]) return ck[n1]<ck[n2];
            if(tarr[n1]!=tarr[n2]) return tarr[n1]<tarr[n2];
            return n1<n2;
        });
        
        auto cid = [&](int ii, int jj, int kk) {return (long(ii)*p->knoy + jj)*p->knoz + kk;};
        
        std::vector<long> ckey;
        std::vector<size_t> cbeg, cend, cnext;
        
        for(size_t q=0; q<order.size(); ++q)
        {
            int nn = order[q];
            long key = cid(ci[nn],cj[nn],ck[nn]);
            
            if(ckey.empty() || ckey.back()!=key)
            {
                ckey.push_back(key);
                cbeg.push_back(q);
                cend.push_back(q+1);
            }
            else
            cend.back() = q+1;
        }
        
        cnext = cbeg;
        
        std::vector<size_t> work(ckey.size());
        std::vector<char> inwork(ckey.size(),1);
        
        for(size_t c=0;c<ckey.size();++c)
        work.at(c) = c;
        
        while(!work.empty())
        {
            size_t c = work.back();
            work.pop_back();
            inwork.at(c)=0;
            
            int n1 = order[cbeg.at(c)];
            double &A = LA(ci[n1],cj[n1],ck[n1]);
            
            while(cnext.at(c)<cend.at(c) && A + 1.0e-9*vpar >= vpar)
            {
                int nn = order[cnext.at(c)];
                
                state[nn]=2;
                ++nacc;
                A -= vpar;
                ++cnext.at(c);
                
                // the source cell gets the volume back
                if(local(si[nn],sj[nn],sk[nn]))
                {
                    LA(si[nn],sj[nn],sk[nn]) += vpar;
                    
                    long skey = cid(si[nn],sj[nn],sk[nn]);
                    auto it = std::lower_bound(ckey.begin(), ckey.end(), skey);
                    
                    if(it!=ckey.end() && *it==skey)
                    {
                        size_t cs = it - ckey.begin();
                        
                        if(inwork.at(cs)==0 && cnext.at(cs)<cend.at(cs))
                        {
                            work.push_back(cs);
                            inwork.at(cs)=1;
                        }
                    }
                }
            }
        }
        
        nacc = pgc->globalisum(nacc);
        
        if(nacc==0)
        break;
    }
    
    // correction of the rejected parcels: no velocity into the full cell;
    // a move across one face stops at the face, a diagonal move stays at the start
    int nrej=0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE && state[n]==1)
    {
        int nface = (ci[n]!=si[n]) + (cj[n]!=sj[n]) + (ck[n]!=sk[n]);
        double f = nface==1 ? MAX(0.0, tarr[n] - 1.0e-6) : 0.0;
        
        X1[n] = X0[n] + f*(X1[n]-X0[n]);
        Y1[n] = Y0[n] + f*(Y1[n]-Y0[n]);
        Z1[n] = Z0[n] + f*(Z1[n]-Z0[n]);
        
        // safeguard at faces (round-off of the cell search): back to the start
        if(p->posc_i(X1[n])!=si[n] || (p->j_dir==1 && p->posc_j(Y1[n])!=sj[n]) || p->posc_k(Z1[n])!=sk[n])
        {
            X1[n] = X0[n];
            Y1[n] = Y0[n];
            Z1[n] = Z0[n];
        }
        
        if(ci[n]!=si[n])
        PU[n] = 0.0;
        
        if(cj[n]!=sj[n])
        PV[n] = 0.0;
        
        if(ck[n]!=sk[n])
        PW[n] = 0.0;
        
        ++nrej;
    }
    
    nrej_step += nrej;
}

// largest occupancy relative to theta_max, diagnostic
double CPM::occupancy_max(lexer *p, ghostcell *pgc)
{
    const double vpar = P.ParcelFactor*Vp;
    
    BASELOOP
    Locc(i,j,k) = 0.0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        i = p->posc_i(P.X[n]);
        j = p->posc_j(P.Y[n]);
        k = p->posc_k(P.Z[n]);
        
        if(i>=0 && i<p->knox && j>=0 && j<p->knoy && k>=0 && k<p->knoz)
        Locc(i,j,k) += vpar;
    }
    
    double omax=0.0;
    
    BASELOOP
    omax = MAX(omax, Locc(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]));
    
    return pgc->globalmax(omax);
}
