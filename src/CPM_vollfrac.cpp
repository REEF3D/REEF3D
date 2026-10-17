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
Authors: Alexander Hanke, Hans Bihs
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

// cell centred tri-linear kernel
// particle contributions to ghost cells at physical boundaries are folded back into the domain,
// contributions to ghost cells at processor boundaries are summed up by start4a_sum
void CPM::kernel(lexer *p, double xs, double ys, double zs, int q)
{
    double wx,wy,wz;
    
    // x, posf_i returns knox in the last half cell at the outer boundary
    i = pfi(p,q,xs);
    
    if(i>=p->knox && p->nb4<0)
    i = p->knox-1;
    
    wx = (p->XP[IP1]-xs)/p->DXP[IP];
    wx = MAX(0.0,MIN(1.0,wx));
    ki[0]=i;
    ki[1]=i+1;
    
    if(ki[0]<0 && p->nb1<0 && perx!=1)
    ki[0]=ki[1];
    
    if(ki[1]>=p->knox && p->nb4<0 && perx!=1)
    ki[1]=ki[0];
    
    // y
    if(p->j_dir==0)
    {
        kj[0]=kj[1]=0;
        wy=1.0;
    }
    
    if(p->j_dir==1)
    {
        j = pfj(p,q,ys);
        
        if(j>=p->knoy && p->nb2<0)
        j = p->knoy-1;
        
        wy = (p->YP[JP1]-ys)/p->DYP[JP];
        wy = MAX(0.0,MIN(1.0,wy));
        kj[0]=j;
        kj[1]=j+1;
        
        if(kj[0]<0 && p->nb3<0 && pery!=1)
        kj[0]=kj[1];
        
        if(kj[1]>=p->knoy && p->nb2<0 && pery!=1)
        kj[1]=kj[0];
    }
    
    // z
    k = pfk(p,q,zs);
    
    if(k>=p->knoz && p->nb6<0)
    k = p->knoz-1;
    
    wz = (p->ZP[KP1]-zs)/p->DZP[KP];
    wz = MAX(0.0,MIN(1.0,wz));
    kk[0]=k;
    kk[1]=k+1;
    
    if(kk[0]<0 && p->nb5<0)
    kk[0]=kk[1];
    
    if(kk[1]>=p->knoz && p->nb6<0)
    kk[1]=kk[0];
    
    kw[0][0][0] = wx*wy*wz;
    kw[1][0][0] = (1.0-wx)*wy*wz;
    kw[0][1][0] = wx*(1.0-wy)*wz;
    kw[1][1][0] = (1.0-wx)*(1.0-wy)*wz;
    kw[0][0][1] = wx*wy*(1.0-wz);
    kw[1][0][1] = (1.0-wx)*wy*(1.0-wz);
    kw[0][1][1] = wx*(1.0-wy)*(1.0-wz);
    kw[1][1][1] = (1.0-wx)*(1.0-wy)*(1.0-wz);
}

bool CPM::wallcell(lexer *p, int ii, int jj, int kkk)
{
    if(perx!=1)
    if((ii<0 && p->nb1<0) || (ii>=p->knox && p->nb4<0))
    return true;
    
    if(p->j_dir==1 && pery!=1)
    if((jj<0 && p->nb3<0) || (jj>=p->knoy && p->nb2<0))
    return true;
    
    if((kkk<0 && p->nb5<0) || (kkk>=p->knoz && p->nb6<0))
    return true;
    
    return false;
}

void CPM::volfrac_update(lexer *p, ghostcell *pgc, sediment_fdm *s, double *PX, double *PY, double *PZ, double *PU, double *PV, double *PW)
{
    int qi,qj,qk;
    double w;
    
    for(i=-1;i<p->knox+1;++i)
    for(j=-1;j<p->knoy+1;++j)
    for(k=-1;k<p->knoz+1;++k)
    {
        cellSum(i,j,k) = 0.0;
        Us(i,j,k) = 0.0;
        Vs(i,j,k) = 0.0;
        Ws(i,j,k) = 0.0;
    }

    // deposit parcel volume and momentum
    for(size_t n=0;n<P.index;n++)
    if(P.Flag[n]>=ACTIVE)
    {
        kernel(p,PX[n],PY[n],PZ[n],int(n));
        
        for(qi=0;qi<2;++qi)
        for(qj=0;qj<2;++qj)
        for(qk=0;qk<2;++qk)
        {
            w = P.ParcelFactor*Vp*kw[qi][qj][qk];
            
            if(w>0.0)
            {
            // the parcels of the sub-grid bedload layer (Q 58) do not move the pore water of the bed
            double wu = P.Hop[n]>0.0 ? 0.0 : w;
            cellSum(ki[qi],kj[qj],kk[qk]) += w;
            Us(ki[qi],kj[qj],kk[qk]) += wu*PU[n];
            Vs(ki[qi],kj[qj],kk[qk]) += wu*PV[n];
            Ws(ki[qi],kj[qj],kk[qk]) += wu*PW[n];
            }
        }
    }
    
    pfold(p,cellSum);
    
    pgc->start4a_sum(p,cellSum,1);
    pfold(p,Us);
    pgc->start4a_sum(p,Us,1);
    pfold(p,Vs);
    pgc->start4a_sum(p,Vs,1);
    pfold(p,Ws);
    pgc->start4a_sum(p,Ws,1);
    
    BASELOOP
    {
        double vol = kvol(p,i,j,k);
        
        Ts(i,j,k) = cellSum(i,j,k)/vol;
        
        if(cellSum(i,j,k)>1.0e-10*vol)
        {
            Us(i,j,k) /= cellSum(i,j,k);
            Vs(i,j,k) /= cellSum(i,j,k);
            Ws(i,j,k) /= cellSum(i,j,k);
        }
        else
        {
            Us(i,j,k) = 0.0;
            Vs(i,j,k) = 0.0;
            Ws(i,j,k) = 0.0;
        }
    }
    
    pgc->start4a(p,Ts,1);
    
    smooth(p,pgc,Ts,p->Q27);
    
    pgc->start4a(p,Us,1);
    pgc->start4a(p,Vs,1);
    pgc->start4a(p,Ws,1);
}

// volume of the kernel of a node: the integral of the tri-linear weights of the node over the domain.
// On a uniform grid it is the cell volume; where the cell size changes, the weights of a node reach
// half way to the neighbouring centres on each side, 0.5 (DXP[i-1] + DXP[i]), and at a wall (weights
// folded back) from the wall to half way to the next centre. Dividing the deposited volume by the cell
// volume instead over-fills the last fine cell at a change 1:2 by 25 % and under-fills the first coarse
// cell by 12.5 %; the contact pressure then pushes the parcels of the bed out of the fine cell.
double CPM::kvol(lexer *p, int ii, int jj, int kk)
{
    auto dir = [&](const double *DN, const double *DP, int q, int nq, bool lo, bool hi) -> double
    {
        double h = DN[q+marge];
        double v = 0.5*(DP[q-1+marge] + DP[q+marge]);
        
        if(q==0 && lo)
        v = 0.5*(DN[q+marge] + DP[q+marge]);
        
        if(q==nq-1 && hi)
        v = 0.5*(DN[q+marge] + DP[q-1+marge]);
        
        if(nq==1 && lo && hi)
        v = h;
        
        return fabs(v-h) <= 1.0e-9*h ? h : v;
    };
    
    double vx = dir(p->DXN, p->DXP, ii, p->knox, p->nb1<0 && perx!=1, p->nb4<0 && perx!=1);
    double vy = p->j_dir==1 ? dir(p->DYN, p->DYP, jj, p->knoy, p->nb3<0 && pery!=1, p->nb2<0 && pery!=1) : p->DYN[jj+marge];
    double vz = dir(p->DZN, p->DZP, kk, p->knoz, p->nb5<0, p->nb6<0);
    
    return vx*vy*vz;
}

// conservative explicit filter, zero flux at physical walls; the exchange between two nodes is
// limited by the smaller kernel volume, so that the solid volume (Ts kvol) is conserved where the
// cell size changes (on a uniform grid the weights are 1)
void CPM::smooth(lexer *p, ghostcell *pgc, field &f, int passes)
{
    double c = p->j_dir==1 ? 1.0/12.0 : 1.0/8.0;
    double val,fc,vc;
    
    auto wgt = [&](int ii, int jj, int kk) -> double
    {
        double vn = kvol(p,ii,jj,kk);
        return vn<vc ? vn/vc : 1.0;
    };
    
    for(int qn=0;qn<passes;++qn)
    {
        BASELOOP
        {
            fc = f(i,j,k);
            vc = kvol(p,i,j,k);
            val = 0.0;
            
            if(!wallcell(p,i-1,j,k))
            val += (f(i-1,j,k) - fc)*wgt(i-1,j,k);
            
            if(!wallcell(p,i+1,j,k))
            val += (f(i+1,j,k) - fc)*wgt(i+1,j,k);
            
            if(p->j_dir==1)
            {
            if(!wallcell(p,i,j-1,k))
            val += (f(i,j-1,k) - fc)*wgt(i,j-1,k);
            
            if(!wallcell(p,i,j+1,k))
            val += (f(i,j+1,k) - fc)*wgt(i,j+1,k);
            }
            
            if(!wallcell(p,i,j,k-1))
            val += (f(i,j,k-1) - fc)*wgt(i,j,k-1);
            
            if(!wallcell(p,i,j,k+1))
            val += (f(i,j,k+1) - fc)*wgt(i,j,k+1);
            
            cellSum(i,j,k) = fc + c*val;
        }
        
        BASELOOP
        f(i,j,k) = cellSum(i,j,k);
        
        pgc->start4a(p,f,1);
    }
}
