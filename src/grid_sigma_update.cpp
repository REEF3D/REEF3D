/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"grid_sigma.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"fnpf_ddx_cds2.h"
#include"fnpf_ddx_cds4.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"
#include"grid_sigma_data.h"

#define WLVL (fabs(d->WL(i,j))>1.0e-20?d->WL(i,j):1.0e20)
#define HX (fabs(d->hx(i,j))>1.0e-20?d->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(d->WL(i,j)+d->WL(i+1,j)))>1.0e-20?0.5*(d->WL(i,j)+d->WL(i+1,j)):1.0e20)
#define HY (fabs(d->hy(i,j))>1.0e-20?d->hy(i,j):1.0e20)

void grid_sigma::sigma_update(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta, slice &eta_n, double alpha)
{
    double wl,sigval;
    double bx,by,ex,ey;
    
    // calculate: Ex,Ey,Exx,Eyy
    // 3D
    if(p->i_dir==1 && p->j_dir==1)
    SLICELOOP4
    {
    pd->Ex(i,j) = pdx->sx(p,eta,1.0);
    pd->Ey(i,j) = pdx->sy(p,eta,1.0);
    
    pd->Exx(i,j) = pddx->sxx(p,eta);
    pd->Eyy(i,j) = pddx->syy(p,eta);
    }
    
    // 2D
    if(p->j_dir==0 && p->A312!=1)
    SLICELOOP4
    {
    pd->Ex(i,j) = pdx->sx(p,eta,1.0);    
    pd->Exx(i,j) = pddx->sxx(p,eta);
    }
    
    double Pval,Qval;
    // 2D
    if(p->j_dir==0 && p->A312==1)
    SLICELOOP4
    {
    k=p->knoz-1;
    Pval = 0.5*(d->U[IJK]+d->U[Im1JK]);
    Qval = 0.5*(d->V[IJK]+d->V[IJm1K]);
    
    if(Pval>=0.0)
    pd->Ex(i,j) = (eta(i,j)-eta(i-1,j))/(p->DXP[IP]);
    
    if(Pval<0.0)
    pd->Ex(i,j) = (eta(i+1,j)-eta(i,j))/(p->DXP[IP]);

    pd->Exx(i,j) = pddx->sxx(p,eta);
    }
    
    pgc->gcsl_start4(p,pd->Ex,1);
    pgc->gcsl_start4(p,pd->Ey,1);
    
    // calculate: Bx,By,Bxx,Byy
    // 3D
    if(p->j_dir==1)
    SLICELOOP4
    {
    pd->Bx(i,j) = pdx->sx(p,d->depth,1.0);
    pd->By(i,j) = pdx->sy(p,d->depth,1.0);
    
    pd->Bxx(i,j) = pddx->sxx(p,d->depth);
    pd->Byy(i,j) = pddx->syy(p,d->depth);
    }

    // 2D
    if(p->j_dir==0 && p->A312!=1)
    SLICELOOP4
    {
    pd->Bx(i,j) = pdx->sx(p,d->depth,1.0);    
    pd->Bxx(i,j) = pddx->sxx(p,d->depth);
    }
    
    if(p->j_dir==0 && p->A312==1)
    SLICELOOP4
    {
    k=0;
    Pval = 0.5*(d->U[IJK]+d->U[Im1JK]);
    
    if(Pval>=0.0)
    pd->Bx(i,j) = (d->depth(i,j)-d->depth(i-1,j))/(p->DXP[IP]);
    
    if(Pval<0.0)
    pd->Bx(i,j) = (d->depth(i+1,j)-d->depth(i,j))/(p->DXP[IP]);
        
        
    pd->Bxx(i,j) = pddx->sxx(p,d->depth);
    }
    
    pgc->gcsl_start4(p,pd->Bx,1);
    pgc->gcsl_start4(p,pd->By,1);
    
    // sigx
    FLOOP
    p->sigx[FIJK] = (1.0 - p->sig[FIJK])*(pd->Bx(i,j)/WLVL) - p->sig[FIJK]*(pd->Ex(i,j)/WLVL);
    
    ULOOP
    {
    sigval = 0.25*(p->sig[FIJK]+p->sig[FIJKp1]+p->sig[FIp1JK]+p->sig[FIp1JKp1]);
    
    bx = (d->depth(i+1,j)-d->depth(i,j))/p->DXP[IP];
    ex = (eta(i+1,j)-eta(i,j))/p->DXP[IP]; 
    
    p->sigx1[IJK] = (1.0 - sigval)*(bx/HX) - sigval*(ex/HX);
    }
    
    LOOP
    {
    sigval = 0.5*(p->sig[FIJK]+p->sig[FIJKp1]);
    
    p->sigx4[IJK] = (1.0 - sigval)*(pd->Bx(i,j)/HXP) - sigval*(pd->Ex(i,j)/HXP);
    }
    
    // sigy
    FLOOP
    p->sigy[FIJK] = (1.0 - p->sig[FIJK])*(pd->By(i,j)/WLVL) - p->sig[FIJK]*(pd->Ey(i,j)/WLVL);
    
    // sigz
    SLICELOOP4
    {
    wl = MAX(0.0, eta(i,j) + p->wd - d->bed(i,j));
    wl = (fabs(wl)>1.0e-20?wl:1.0e20);
    
    p->sigz[IJ] = 1.0/wl;
    }

    // sigt
    FLOOP
    p->sigt[FIJK] = -(p->sig[FIJK]/WLVL)*(d->WL(i,j)-d->WL_n(i,j))/(p->dt);

    // sigxx
    FLOOP
    {
        p->sigxx[FIJK] = ((1.0 - p->sig[FIJK])/WLVL)*(pd->Bxx(i,j) - pow(pd->Bx(i,j),2.0)/WLVL) // xx
        
                      - (p->sig[FIJK]/WLVL)*(pd->Exx(i,j) - pow(pd->Ex(i,j),2.0)/WLVL)
                      
                      - (p->sigx[FIJK]/WLVL)*(pd->Bx(i,j) + pd->Ex(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(pd->Bx(i,j)*pd->Ex(i,j))
                      
                      
                      + ((1.0 - p->sig[FIJK])/WLVL)*(pd->Byy(i,j) - pow(pd->By(i,j),2.0)/WLVL) // yy
        
                      - (p->sig[FIJK]/WLVL)*(pd->Eyy(i,j) - pow(pd->Ey(i,j),2.0)/WLVL)
                      
                      - (p->sigy[FIJK]/WLVL)*(pd->By(i,j) + pd->Ey(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(pd->By(i,j)*pd->Ey(i,j));
    }
    
    
    // sig BC
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigx[FIJKm1] = p->sigx[FIJK];
            p->sigx[FIJKm2] = p->sigx[FIJK];
            p->sigx[FIJKm3] = p->sigx[FIJK];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigx[FIJKp1] = p->sigx[FIJK];
            p->sigx[FIJKp2] = p->sigx[FIJK];
            p->sigx[FIJKp3] = p->sigx[FIJK];
        } 
    }
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigy[FIJKm1] = p->sigy[FIJK];
            p->sigy[FIJKm2] = p->sigy[FIJK];
            p->sigy[FIJKm3] = p->sigy[FIJK];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigy[FIJKp1] = p->sigy[FIJK];
            p->sigy[FIJKp2] = p->sigy[FIJK];
            p->sigy[FIJKp3] = p->sigy[FIJK];
        } 
    }
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigxx[FIJKm1] = p->sigxx[FIJK];
            p->sigxx[FIJKm2] = p->sigxx[FIJK];
            p->sigxx[FIJKm3] = p->sigxx[FIJK];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigxx[FIJKp1] = p->sigxx[FIJK];
            p->sigxx[FIJKp2] = p->sigxx[FIJK];
            p->sigxx[FIJKp3] = p->sigxx[FIJK];
        } 
    }
    
    
    FLOOP
    p->ZSN[FIJK] = p->ZN[KP]*d->WL(i,j) + d->bed(i,j);
    
    
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*d->WL(i,j) + d->bed(i,j);

    
    pgc->start7S(p,p->sigx,1);
    pgc->start7S(p,p->sigy,1);
    pgc->start7S(p,p->sigxx,1);
    pgc->start7S(p,p->sigt,1);
    pgc->start7S(p,p->ZSN,1);
    pgc->gcslparaxijk(p, p->sigz, 1);
    
}

void grid_sigma::omega_update(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *W, slice &eta, slice &eta_n, double alpha)
{ 
    double wval,Pval,Qval,Rval;
    
    WLOOP
    {
    if(0.5*(U[Im1JKp1] + U[IJKp1])>=0.0)
    Pval=0.5*(U[Im1JK] + U[Im1JKp1]);
        
    if(0.5*(U[Im1JKp1] + U[IJKp1])<0.0)
    Pval=0.5*(U[IJK] + U[IJK+1]);
    
    
    if(0.5*(V[IJm1Kp1] + V[IJKp1])>=0.0)
    Qval=0.5*(V[IJm1K] + V[IJm1Kp1]);
        
    if(0.5*(V[IJm1Kp1] + V[IJKp1])<0.0)
    Qval=0.5*(V[IJK] + V[IJKp1]);
    
    
    
    if(W[IJK]>=0.0)
    Rval=0.5*(W[IJK] + W[IJKm1]);
        
    if(W[IJK]<0.0)
    Rval=0.5*(W[IJK] + W[IJKp1]);
    
    
        
    d->omega[IJK] =  p->sigt[FIJKp1]
                    
                    +  Pval*p->sigx[FIJKp1]
                    
                    +  Qval*p->sigy[FIJKp1]
                    
                    +  Rval*p->sigz[IJ];
                    
    }
    
    GC3LOOP
    if(p->gcb3[n][3]==6 && p->gcb3[n][4]==3)
    {
    i=p->gcb3[n][0];
    j=p->gcb3[n][1];
    k=p->gcb3[n][2];
    
        
        if(p->A516==1)
        for(int q=0;q<3;++q)
        d->omega[IJKp1+q] =  0.0;
        
        if(p->A516==2)
        for(int q=0;q<3;++q)
        d->omega[IJKp1+q] =   p->sigt[FIJKp2]
                    
                    +  0.5*(U[Im1JKp1] + U[IJKp1])*p->sigx[FIJKp2]
                    
                    +  0.5*(V[IJm1Kp1] + V[IJKp1])*p->sigy[FIJKp2]
                    
                    +  W[IJKp1]*p->sigz[IJ];
                    
        if(p->A516==3)
        for(int q=0;q<3;++q)
        {
        if(0.5*(U[Im1JKp1] + U[IJKp1])>=0.0)
        Pval=U[Im1JKp1];
        
        if(0.5*(U[Im1JKp1] + U[IJKp1])<0.0)
        Pval=U[IJKp1];
        
        d->omega[IJKp1+q] =   p->sigt[FIJKp2]
                    
                    +  Pval*p->sigx[FIJKp2]
                    
                    +  0.5*(V[IJm1Kp1] + V[IJKp1])*p->sigy[FIJKp2]
                    
                    +  W[IJKp1]*p->sigz[IJ];
        }
                    
        if(p->A516==4)
        for(int q=0;q<3;++q)
        d->omega[IJKp1+q] =  d->omega[IJK];
    }
    
    GC3LOOP
    if(p->gcb3[n][3]==5 && p->gcb3[n][4]==21)
    {
    i=p->gcb3[n][0];
    j=p->gcb3[n][1];
    k=p->gcb3[n][2];
    
        for(int q=0;q<3;++q)
        d->omega[IJKm1-q] =  0.0;
    }
    
    //pgc->start3(p,d->omega,17);
    //pgc->start3(p,d->omega,17);
}



