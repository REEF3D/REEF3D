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
void CPM::kernel(lexer *p, double xs, double ys, double zs)
{
    double wx,wy,wz;
    
    // x, posf_i returns knox in the last half cell at the outer boundary
    i = p->posf_i(xs);
    
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
        j = p->posf_j(ys);
        
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
    k = p->posf_k(zs);
    
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
        kernel(p,PX[n],PY[n],PZ[n]);
        
        for(qi=0;qi<2;++qi)
        for(qj=0;qj<2;++qj)
        for(qk=0;qk<2;++qk)
        {
            w = P.ParcelFactor*Vp*kw[qi][qj][qk];
            
            if(w>0.0)
            {
            cellSum(ki[qi],kj[qj],kk[qk]) += w;
            Us(ki[qi],kj[qj],kk[qk]) += w*PU[n];
            Vs(ki[qi],kj[qj],kk[qk]) += w*PV[n];
            Ws(ki[qi],kj[qj],kk[qk]) += w*PW[n];
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
        double vol = p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        
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

// conservative explicit filter, zero flux at physical walls
void CPM::smooth(lexer *p, ghostcell *pgc, field &f, int passes)
{
    double c = p->j_dir==1 ? 1.0/12.0 : 1.0/8.0;
    double val,fc;
    
    for(int qn=0;qn<passes;++qn)
    {
        BASELOOP
        {
            fc = f(i,j,k);
            val = 0.0;
            
            if(!wallcell(p,i-1,j,k))
            val += f(i-1,j,k) - fc;
            
            if(!wallcell(p,i+1,j,k))
            val += f(i+1,j,k) - fc;
            
            if(p->j_dir==1)
            {
            if(!wallcell(p,i,j-1,k))
            val += f(i,j-1,k) - fc;
            
            if(!wallcell(p,i,j+1,k))
            val += f(i,j+1,k) - fc;
            }
            
            if(!wallcell(p,i,j,k-1))
            val += f(i,j,k-1) - fc;
            
            if(!wallcell(p,i,j,k+1))
            val += f(i,j,k+1) - fc;
            
            cellSum(i,j,k) = fc + c*val;
        }
        
        BASELOOP
        f(i,j,k) = cellSum(i,j,k);
        
        pgc->start4a(p,f,1);
    }
}
