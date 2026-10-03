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

Exchanges: all requests are accepted first; a cell that ends up over capacity rejects
its latest arrivals until it fits, the rejected parcels stay in their source cells, and
this is repeated until no cell is over capacity (the accepted moves only decrease, the
start state fits, so it terminates). A parcel may thus enter a full cell while another
leaves it: a full row of cells can shear as a conveyor (no deadlock of the sheared bed).

Parallel: for cells with requests from other ranks, the requests are binned by arrival
time into NB bins and summed into the owner of the cell with start4a_sum, the owner
rejects whole bins (latest first) and returns the cut-off by a ghost cell exchange. Cells
with local requests only reject exactly in the order of arrival. No parcel lists cross
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
    
    // cell index of a position; beyond the local domain posc_i/j/k return knox+1 (knoy+1, knoz+1),
    // a move of less than one cell ends in the first ghost layer: knox (the neighbour's first cell,
    // or for a serial periodic side the copy of cell 0, see pfold)
    auto ci_of = [&](double x) {int q = p->posc_i(x); return q>p->knox ? p->knox : (q<-1 ? -1 : q);};
    auto cj_of = [&](double y) {int q = p->posc_j(y); return q>p->knoy ? p->knoy : (q<-1 ? -1 : q);};
    auto ck_of = [&](double z) {int q = p->posc_k(z); return q>p->knoz ? p->knoz : (q<-1 ? -1 : q);};
    
    // occupancy at the start of the step
    for(i=-1;i<p->knox+1;++i)
    for(j=-1;j<p->knoy+1;++j)
    for(k=-1;k<p->knoz+1;++k)
    Locc(i,j,k) = 0.0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]>=ACTIVE)
    {
        i0 = ci_of(X0[n]);
        j0 = p->j_dir==1 ? cj_of(Y0[n]) : 0;
        k0 = ck_of(Z0[n]);
        
        if(local(i0,j0,k0))
        Locc(i0,j0,k0) += vpar;
    }
    
    // requests
    int nreq=0;
    int nclip=0;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        i0 = ci_of(X0[n]);
        j0 = p->j_dir==1 ? cj_of(Y0[n]) : 0;
        k0 = ck_of(Z0[n]);
        
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
        
        i1 = ci_of(X1[n]);
        j1 = p->j_dir==1 ? cj_of(Y1[n]) : 0;
        k1 = ck_of(Z1[n]);
        
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
    
    // optimistic acceptance with fix-up (exchanges allowed):
    // all requests are accepted first, then cells that end up over capacity reject their latest
    // arrivals until they fit; a rejection keeps the parcel in its source cell, which may then
    // be over capacity in turn, so the fix-up is repeated until no cell is over capacity.
    // The accepted moves only decrease, the start state is within capacity, so it terminates.
    // A full row of cells can shear as a conveyor: a parcel enters a full cell while another leaves.
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE && state[n]==1)
    state[n]=2;
    
    const int MAXIT = 50;
    int it;
    
    for(it=0; it<MAXIT; ++it)
    {
        for(i=-1;i<p->knox+1;++i)
        for(j=-1;j<p->knoy+1;++j)
        for(k=-1;k<p->knoz+1;++k)
        {
        Lout(i,j,k) = 0.0;
        Lin(i,j,k) = 0.0;
        for(int b=0;b<NB;++b)
        (*hist[b])(i,j,k) = 0.0;
        }
        
        // volume leaving (local source) and arriving (any target, summed to the owner);
        // arrivals that can be rejected are counted per arrival time bin
        for(n=0;n<P.index;++n)
        if(P.Flag[n]==ACTIVE && (state[n]==2 || state[n]==4))
        {
            if(local(si[n],sj[n],sk[n]))
            Lout(si[n],sj[n],sk[n]) += vpar;
            
            Lin(ci[n],cj[n],ck[n]) += vpar;
            
            if(state[n]==2)
            (*hist[bin(tarr[n])])(ci[n],cj[n],ck[n]) += 1.0;
        }
        
        BASELOOP
        Lloc(i,j,k) = (*hist[0])(i,j,k) + (*hist[1])(i,j,k) + (*hist[2])(i,j,k) + (*hist[3])(i,j,k);
        
        for(int b=0;b<NB;++b)
        {
        pfold(p,*hist[b]);
        pgc->start4a_sum(p,*hist[b],1);
        }
        
        pfold(p,Lin);
        pgc->start4a_sum(p,Lin,1);
        
        // excess per cell, cut-off: arrivals in bins above Ltc are rejected
        //   Ltc = 2 NB: fits, nothing rejected
        //   Ltc > NB (= 1.5 NB): local arrivals only, the excess is rejected exactly (latest first)
        //   Ltc in [-1,NB-1]: arrivals from other ranks, whole bins are rejected
        int nover=0;
        
        BASELOOP
        {
            double V = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
            // capacity: packing of the contact network theta_0(I) (theta_0 at rest, lower in sheared layers)
            double t0 = p->Q12==2 ? T0e(i,j,k) : theta_0;
            double cap = MAX((t0 + theta_max - theta_0)*V, t0*V + vpar*(1.0+1.0e-6));
            
            if(a->solid(i,j,k)<0.0)
            cap = 0.0;
            
            double excess = Locc(i,j,k) - Lout(i,j,k) + Lin(i,j,k) - cap;
            
            LA(i,j,k) = excess;
            Ltc(i,j,k) = 2.0*NB;
            
            if(excess > 1.0e-9*vpar)
            {
                double cnt = (*hist[0])(i,j,k) + (*hist[1])(i,j,k) + (*hist[2])(i,j,k) + (*hist[3])(i,j,k);
                
                if(cnt<=0.0)
                continue;
                
                ++nover;
                
                if(cnt>Lloc(i,j,k)+0.5)
                {
                    // reject whole bins, latest first, until the excess is covered
                    double need = ceil(excess/vpar - 1.0e-9);
                    double rej = 0.0;
                    double tc = double(NB-1);
                    
                    for(int b=NB-1;b>=0;--b)
                    {
                        if(rej>=need)
                        break;
                        
                        rej += (*hist[b])(i,j,k);
                        tc = double(b-1);
                    }
                    
                    Ltc(i,j,k) = tc;
                }
                else
                Ltc(i,j,k) = 1.5*NB;
            }
        }
        
        nover = pgc->globalisum(nover);
        
        if(nover==0)
        break;
        
        pgc->start4a(p,Ltc,1);
        
        // reject: bins above the cut-off, exact for cells with local arrivals only
        order.clear();
        
        for(n=0;n<P.index;++n)
        if(P.Flag[n]==ACTIVE && state[n]==2)
        {
            double tc = Ltc(ci[n],cj[n],ck[n]);
            
            if(tc>double(NB) && tc<2.0*NB-0.5 && local(ci[n],cj[n],ck[n]))
            order.push_back(n);
            
            else if(tc<double(NB) && double(bin(tarr[n]))>tc+1.0e-9)
            state[n]=1;
        }
        
        // exact: latest arrivals first, ties by index
        std::sort(order.begin(), order.end(), [&](int n1, int n2)
        {
            if(ci[n1]!=ci[n2]) return ci[n1]<ci[n2];
            if(cj[n1]!=cj[n2]) return cj[n1]<cj[n2];
            if(ck[n1]!=ck[n2]) return ck[n1]<ck[n2];
            if(tarr[n1]!=tarr[n2]) return tarr[n1]>tarr[n2];
            return n1>n2;
        });
        
        for(size_t q=0; q<order.size(); ++q)
        {
            int nn = order[q];
            double &E = LA(ci[nn],cj[nn],ck[nn]);
            
            if(E > 1.0e-9*vpar)
            {
                state[nn]=1;
                E -= vpar;
            }
        }
    }
    
    nit_step = MAX(nit_step,it);
    
    // still over capacity after MAXIT: back to the start for all limited moves
    if(it==MAXIT)
    {
        if(p->mpirank==0)
        cout<<"CPM limiter: no consistent set of moves after "<<MAXIT<<" passes, limited moves rejected"<<endl;
        
        for(n=0;n<P.index;++n)
        if(P.Flag[n]==ACTIVE && state[n]==2)
        state[n]=1;
    }
    
    // correction of the rejected parcels: no velocity into the full cell;
    // a move across one face stops at the face, a diagonal move stays at the start
    int nrej=0;
    double gmag = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE && state[n]==1)
    {
        int nface = (ci[n]!=si[n]) + (cj[n]!=sj[n]) + (ck[n]!=sk[n]);
        double f = nface==1 ? MAX(0.0, tarr[n] - 1.0e-6) : 0.0;
        
        X1[n] = X0[n] + f*(X1[n]-X0[n]);
        Y1[n] = Y0[n] + f*(Y1[n]-Y0[n]);
        Z1[n] = Z0[n] + f*(Z1[n]-Z0[n]);
        
        // safeguard at faces (round-off of the cell search): back to the start
        if(ci_of(X1[n])!=si[n] || (p->j_dir==1 && cj_of(Y1[n])!=sj[n]) || ck_of(Z1[n])!=sk[n])
        {
            X1[n] = X0[n];
            Y1[n] = Y0[n];
            Z1[n] = Z0[n];
        }
        
        // ride-over (Q 56): a grain pushed against a full cell across the gravity direction
        // climbs over the grains ahead, part of the blocked velocity turns upward
        //   u_up = Q56 |u_blocked|   (Q56 = tan of the pivot angle)
        double ublock = 0.0;
        
        if(ci[n]!=si[n] && fabs(p->W20)<0.5*gmag)
        ublock = MAX(ublock,fabs(PU[n]));
        
        if(cj[n]!=sj[n] && fabs(p->W21)<0.5*gmag)
        ublock = MAX(ublock,fabs(PV[n]));
        
        if(ck[n]!=sk[n] && fabs(p->W22)<0.5*gmag)
        ublock = MAX(ublock,fabs(PW[n]));
        
        if(ci[n]!=si[n])
        PU[n] = 0.0;
        
        if(cj[n]!=sj[n])
        PV[n] = 0.0;
        
        if(ck[n]!=sk[n])
        PW[n] = 0.0;
        
        if(p->Q56>0.0 && ublock>0.0 && gmag>1.0e-10)
        {
            double up = p->Q56*ublock;
            double ex = -p->W20/gmag, ey = -p->W21/gmag, ez = -p->W22/gmag;
            double un = PU[n]*ex + PV[n]*ey + PW[n]*ez;
            
            if(un<up)
            {
                PU[n] += (up-un)*ex;
                PV[n] += (up-un)*ey;
                PW[n] += (up-un)*ez;
            }
        }
        
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
