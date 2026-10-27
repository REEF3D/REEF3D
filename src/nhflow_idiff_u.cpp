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

#include"nhflow_idiff.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_thinbody.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"solver.h"
#include"slice.h"
#include"vrans.h"
#include"nhflow_curv.h"

nhflow_idiff::nhflow_idiff(lexer* p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=14;
	gcval_vh=15;
	gcval_wh=16;
}

nhflow_idiff::~nhflow_idiff()
{
}

void nhflow_idiff::diff_u(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, solver *psolv, double *UHdiff, double *UHin, double *UH, double *VH, double *WH, slice &WL, double alpha)
{
	starttime=pgc->timer();
    
    double visc_IP,visc_IM1,visc_IP1,visc_JM1,visc_JP1,visc_KM1,visc_KP1;
    double visc_IPh,visc_KM1h,visc_KP1h,cvt,cvb;
    
    // HLES split (A 561): horizontal stresses (xx, yy, xy) with nu_h = EVH, stresses with a z index with nu_v = EV
    double *EVh = (p->A561>0) ? d->EVH : d->EV;

    LOOP
    UHdiff[IJK] = UHin[IJK];
    
    pgc->start4V(p,UHdiff,gcval_uh);
    
    pflow->rkinflow_nhflow(p,d,pgc,UHdiff,UHin);
    
    if(p->geo_curv==1)
    nhflow_curv::diff_div(p,UH,VH);

    n=0;
    LOOP
	{
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {  
            visc_IP = d->VISC[IJK] + d->EV[IJK];
            visc_IPh = d->VISC[IJK] + EVh[IJK];
            
            visc_IM1 = 0.5*(visc_IPh + d->VISC[Im1JK] + EVh[Im1JK]);
            visc_IP1 = 0.5*(visc_IPh + d->VISC[Ip1JK] + EVh[Ip1JK]);
            
            visc_JM1 = 0.5*(visc_IPh + d->VISC[IJm1K] + EVh[IJm1K]);
            visc_JP1 = 0.5*(visc_IPh + d->VISC[IJp1K] + EVh[IJp1K]);
            
            visc_KM1 = 0.5*(visc_IP + d->VISC[IJKm1] + d->EV[IJKm1]);
            visc_KP1 = 0.5*(visc_IP + d->VISC[IJKp1] + d->EV[IJKp1]);
            
            visc_KM1h = 0.5*(visc_IPh + d->VISC[IJKm1] + EVh[IJKm1]);
            visc_KP1h = 0.5*(visc_IPh + d->VISC[IJKp1] + EVh[IJKp1]);
            
            // breaking viscosity vb, face values the mean of the two cells (both cells of a face use the same coefficient):
            //  A 550 1: in the top layer only; the face to the layer below takes vb/2 from both sides (k = knoz-1, knoz-2)
            //  A 550 3: in all layers on all faces, the eddy viscosity of Kennedy et al. (2000) as a turbulent viscosity of
            //           the column (with the horizontal faces only, the surf zone depended on the number of layers)
            if(p->A550==3 || k==p->knoz-1)
            {
            visc_IP +=  d->vb(i,j);
            visc_IPh += d->vb(i,j);
            visc_IM1 += 0.5*(d->vb(i,j) + d->vb(i-1,j));
            visc_IP1 += 0.5*(d->vb(i,j) + d->vb(i+1,j));
            visc_JM1 += 0.5*(d->vb(i,j) + d->vb(i,j-1));
            visc_JP1 += 0.5*(d->vb(i,j) + d->vb(i,j+1));
            }

            if(p->A550!=3 && k==p->knoz-1)
            {
            visc_KM1 += 0.5*d->vb(i,j);
            visc_KP1 += d->vb(i,j);
            visc_KM1h += 0.5*d->vb(i,j);
            visc_KP1h += d->vb(i,j);
            }

            if(p->A550!=3 && k==p->knoz-2)
            {
            visc_KP1 += 0.5*d->vb(i,j);
            visc_KP1h += 0.5*d->vb(i,j);
            }

            if(p->A550==3)
            {
            visc_KM1 += d->vb(i,j);
            visc_KP1 += d->vb(i,j);
            visc_KM1h += d->vb(i,j);
            visc_KP1h += d->vb(i,j);
            }

            
            // vertical coefficient at the faces k-1/2 (sigma metrics at FIJK) and k+1/2 (FIJKp1)
            sigxyz2_b = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            sigxyz2_t = pow(p->sigx[FIJKp1],2.0) + pow(p->sigy[FIJKp1],2.0) + pow(p->sigz[IJ],2.0);
            
            // vertical faces: the sigx, sigy parts belong to the horizontal stresses (nu_h), the sigz part to the
            // vertical one (nu_v); identical to visc*sigxyz2 without the split
            if(p->A561>0)
            {
            cvt = visc_KP1h*(pow(p->sigx[FIJKp1],2.0) + pow(p->sigy[FIJKp1],2.0)) + visc_KP1*pow(p->sigz[IJ],2.0);
            cvb = visc_KM1h*(pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0)) + visc_KM1*pow(p->sigz[IJ],2.0);
            }
            else
            {
            cvt = visc_KP1*sigxyz2_t;
            cvb = visc_KM1*sigxyz2_b;
            }
            
            // curvilinear grid (G 50): compact horizontal rows, the rest explicit (nhflow_curv.h)
            if(p->geo_curv==1)
            {
            double cIp,cIm,cJp,cJm;
            nhflow_curv::lapcoef(p,i,j,cIp,cIm,cJp,cJm);
            
            d->M.p[n]  =  visc_IP1*cIp + visc_IM1*cIm + visc_JP1*cJp + visc_JM1*cJm
                        
                        + cvt/(p->DZP[KP]*p->DZN[KP])
                        + cvb/(p->DZP[KM1]*p->DZN[KP])
                        
                        + CPORNH/(alpha*p->dt);
            
            d->M.n[n] = -visc_IP1*cIp;
            d->M.s[n] = -visc_IM1*cIm;
            d->M.w[n] = -visc_JP1*cJp;
            d->M.e[n] = -visc_JM1*cJm;
            
            d->M.t[n] = -cvt/(p->DZP[KP]*p->DZN[KP])
                        - visc_KP1h*p->sigxx[FIJK]/((p->DZP[KP]+p->DZP[KM1]));
                        
            d->M.b[n] = -cvb/(p->DZP[KM1]*p->DZN[KP])
                        + visc_KM1h*p->sigxx[FIJK]/((p->DZP[KP]+p->DZP[KM1]));
            
            d->rhsvec.V[n] = (CPORNH*UHin[IJK])/(alpha*p->dt) + (p->A561>0 ? nhflow_curv::diff_explicit_hv(p,1,i,j,k,UH,VH,WH,visc_IPh,visc_IP) : visc_IP*nhflow_curv::diff_explicit(p,1,i,j,k,UH,VH,WH));
            }
            
            if(p->geo_curv==0)
            {
            d->M.p[n]  =  2.0*visc_IP1/(p->DXP[IP]*p->DXN[IP])
                        + 2.0*visc_IM1/(p->DXP[IM1]*p->DXN[IP])
                        
                        + visc_JP1/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                        + visc_JM1/(p->DYP[JM1]*p->DYN[JP])*p->y_dir
                        
                        + cvt/(p->DZP[KP]*p->DZN[KP])
                        + cvb/(p->DZP[KM1]*p->DZN[KP])
                        
                        + CPORNH/(alpha*p->dt);
            

            d->M.n[n] = -2.0*visc_IP1/(p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -2.0*visc_IM1/(p->DXP[IM1]*p->DXN[IP]);

            d->M.w[n] = -visc_JP1/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
            d->M.e[n] = -visc_JM1/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

            d->M.t[n] = -cvt/(p->DZP[KP]*p->DZN[KP])
                        - visc_KP1h*p->sigxx[FIJK]/((p->DZP[KP]+p->DZP[KM1]));
                        
            d->M.b[n] = -cvb/(p->DZP[KM1]*p->DZN[KP])
                        + visc_KM1h*p->sigxx[FIJK]/((p->DZP[KP]+p->DZP[KM1]));
            
            
            d->rhsvec.V[n] = visc_IPh*((VH[Ip1Jp1K]-VH[Im1Jp1K]) - (VH[Ip1Jm1K]-VH[Im1Jm1K]))/((p->DXP[IP]+p->DXP[IM1])*(p->DYP[JP]+p->DYP[JM1]))*p->y_dir
                           + visc_IP*((WH[Ip1JKp1]-WH[Im1JKp1]) - (WH[Ip1JKm1]-WH[Im1JKm1]))/((p->DXP[IP]+p->DXP[IM1])*(p->DZP[KP]+p->DZP[KM1]))*p->sigz[IJ]


						 + (CPORNH*UHin[IJK])/(alpha*p->dt)
                            
                            
                            + visc_IPh*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(UH[Ip1JKp1] - UH[Im1JKp1] - UH[Ip1JKm1] + UH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZP[KP]+p->DZP[KM1]))
                        
                            + visc_IPh*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(UH[IJp1Kp1] - UH[IJm1Kp1] - UH[IJp1Km1] + UH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZP[KP]+p->DZP[KM1]))*p->y_dir;
            }
        }
        
        if(p->wet[IJ]==0 || p->flag4[IJK]<0 || p->DF[IJK]<0)
        {
        d->M.p[n]  =  1.0;

        d->M.n[n] = 0.0;
        d->M.s[n] = 0.0;

        d->M.w[n] = 0.0;
        d->M.e[n] = 0.0;

        d->M.t[n] = 0.0;
        d->M.b[n] = 0.0;
        
        d->rhsvec.V[n] =  0.0;
        }
        
	++n;
	}
    
    if(p->A513==1)
    {
    n=0;
	LOOP
	{
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
            if(p->flag4[Im1JK]<0 || p->DF[Im1JK]<0)
            {
            d->rhsvec.V[n] -= d->M.s[n]*UH[IJK];
            d->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0)
            {
            d->rhsvec.V[n] -= d->M.n[n]*UH[IJK];
            d->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0)
            {
            d->rhsvec.V[n] -= d->M.e[n]*UH[IJK]*p->y_dir;
            d->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0)
            {
            d->rhsvec.V[n] -= d->M.w[n]*UH[IJK]*p->y_dir;
            d->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
            {
            d->rhsvec.V[n] -= d->M.b[n]*UH[IJK];
            d->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
            {
            d->rhsvec.V[n] -= d->M.t[n]*UH[IJK];
            d->M.t[n] = 0.0;
            }
        }
	++n;
	}
    }
    
    // A 513 2: neighbour values; at a direct forcing solid (DF < 0 inside the domain) the mirrored value of the
    // cell, as in the ghost cells of the domain walls (start4V): tangential symmetric, normal antisymmetric.
    // The forced velocity (0) of the solid cell would make the solid a no-slip wall, the domain walls are slip
    // walls (wall friction only with A 519 2, nhflow_wall.h).
    if(p->A513==2)
    {
    n=0;
	LOOP
	{
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
            if(p->flag4[Im1JK]<0 || p->DF[Im1JK]<0)
            {
            d->rhsvec.V[n] -= d->M.s[n]*(p->DF[Im1JK]<0 && p->flag4[Im1JK]>0 ? -UH[IJK] : UH[Im1JK]);
            d->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0 || p->DF[Ip1JK]<0)
            {
            d->rhsvec.V[n] -= d->M.n[n]*(p->DF[Ip1JK]<0 && p->flag4[Ip1JK]>0 ? -UH[IJK] : UH[Ip1JK]);
            d->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0 || p->DF[IJm1K]<0)
            {
            d->rhsvec.V[n] -= d->M.e[n]*(p->DF[IJm1K]<0 && p->flag4[IJm1K]>0 ? UH[IJK] : UH[IJm1K])*p->y_dir;
            d->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0 || p->DF[IJp1K]<0)
            {
            d->rhsvec.V[n] -= d->M.w[n]*(p->DF[IJp1K]<0 && p->flag4[IJp1K]>0 ? UH[IJK] : UH[IJp1K])*p->y_dir;
            d->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0 || p->DF[IJKm1]<0)
            {
            d->rhsvec.V[n] -= d->M.b[n]*(p->DF[IJKm1]<0 && p->flag4[IJKm1]>0 ? UH[IJK] : UH[IJKm1]);
            d->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0 || p->DF[IJKp1]<0)
            {
            d->rhsvec.V[n] -= d->M.t[n]*UH[IJK];
            d->M.t[n] = 0.0;
            }
        }
	++n;
	}
    }
	
    // sharp thin bodies (X 330 'mobility sharp'): no diffusive exchange through the body
    if(d->thinbody!=nullptr)
    d->thinbody->matrix_walls(p,d,UH);
    
    psolv->startV(p,pgc,UHdiff,d->rhsvec,d->M,4);
    
    pgc->start4V(p,UHdiff,gcval_uh);
    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(4)<<time<<endl;
}
