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

#include"nhflow_idiff_2D.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_breaking_profile.h"
#include"nhflow_thinbody.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"solver.h"
#include"slice.h"
#include"vrans.h"

nhflow_idiff_2D::nhflow_idiff_2D(lexer* p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=14;
	gcval_vh=15;
	gcval_wh=16;
}

nhflow_idiff_2D::~nhflow_idiff_2D()
{
}

void nhflow_idiff_2D::diff_u(lexer *p, fdm_nhf *d, ghostcell *pgc, ioflow *pflow, solver *psolv, double *UHdiff, double *UHin, double *UH, double *VH, double *WH, slice &WL, double alpha)
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

    n=0;
    LOOP
    {
        if(p->wet[IJ]==1 && p->DF[IJK]>0)
        {
            visc_IP = d->VISC[IJK] + d->EV[IJK];
            visc_IPh = d->VISC[IJK] + EVh[IJK];
            
            visc_IM1 = 0.5*(visc_IPh + d->VISC[Im1JK] + EVh[Im1JK]);
            visc_IP1 = 0.5*(visc_IPh + d->VISC[Ip1JK] + EVh[Ip1JK]);
            
            visc_KM1 = 0.5*(visc_IP + d->VISC[IJKm1] + d->EV[IJKm1]);
            visc_KP1 = 0.5*(visc_IP + d->VISC[IJKp1] + d->EV[IJKp1]);
            
            visc_KM1h = 0.5*(visc_IPh + d->VISC[IJKm1] + EVh[IJKm1]);
            visc_KP1h = 0.5*(visc_IPh + d->VISC[IJKp1] + EVh[IJKp1]);
            
            // breaking viscosity vb, face values the mean of the two cells (both cells of a face use the same coefficient):
            //  A 550 1: in the top layer only; the face to the layer below takes vb/2 from both sides (k = knoz-1, knoz-2)
            //  A 550 3: in all layers on all faces, the eddy viscosity of Kennedy et al. (2000) as a turbulent viscosity of
            //           the column (with the horizontal faces only, the surf zone depended on the number of layers)
            //  A 550 4: as A 550 3 with the depth profile of the orbital velocity (nhflow_breaking_profile.h)
            if(p->A550>=3 || k==p->knoz-1)
            {
            visc_IP +=  nhflow_vbz(p,d,i,j,p->ZP[KP]);
            visc_IPh += nhflow_vbz(p,d,i,j,p->ZP[KP]);
            visc_IM1 += 0.5*(nhflow_vbz(p,d,i,j,p->ZP[KP]) + nhflow_vbz(p,d,i-1,j,p->ZP[KP]));
            visc_IP1 += 0.5*(nhflow_vbz(p,d,i,j,p->ZP[KP]) + nhflow_vbz(p,d,i+1,j,p->ZP[KP]));
            }

            if(p->A550<3 && k==p->knoz-1)
            {
            visc_KM1 += 0.5*d->vb(i,j);
            visc_KP1 += d->vb(i,j);
            visc_KM1h += 0.5*d->vb(i,j);
            visc_KP1h += d->vb(i,j);
            }

            if(p->A550<3 && k==p->knoz-2)
            {
            visc_KP1 += 0.5*d->vb(i,j);
            visc_KP1h += 0.5*d->vb(i,j);
            }

            if(p->A550>=3)
            {
            visc_KM1 += nhflow_vbz(p,d,i,j,p->ZN[KP]);
            visc_KP1 += nhflow_vbz(p,d,i,j,p->ZN[KP1]);
            visc_KM1h += nhflow_vbz(p,d,i,j,p->ZN[KP]);
            visc_KP1h += nhflow_vbz(p,d,i,j,p->ZN[KP1]);
            }

            
            // vertical coefficient at the faces k-1/2 (FIJK) and k+1/2 (FIJKp1)
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
            
            
            d->M.p[n]  =   2.0*visc_IP1/(p->DXP[IP]*p->DXN[IP])
                        + 2.0*visc_IM1/(p->DXP[IM1]*p->DXN[IP])
                        
                        + cvt/(p->DZP[KP]*p->DZN[KP])
                        + cvb/(p->DZP[KM1]*p->DZN[KP])
                        
                        + CPORNH/(alpha*p->dt);
            

            d->M.n[n] = -2.0*visc_IP1/(p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -2.0*visc_IM1/(p->DXP[IM1]*p->DXN[IP]);

            d->M.t[n] = -cvt/(p->DZP[KP]*p->DZN[KP])
                        - visc_KP1h*p->sigxx[FIJK]/((p->DZP[KP]+p->DZP[KM1]));
                        
            d->M.b[n] = -cvb/(p->DZP[KM1]*p->DZN[KP])
                        + visc_KM1h*p->sigxx[FIJK]/((p->DZP[KP]+p->DZP[KM1]));
            
            
            d->rhsvec.V[n] = visc_IP*((WH[Ip1JKp1]-WH[Im1JKp1]) - (WH[Ip1JKm1]-WH[Im1JKm1]))/((p->DXP[IP]+p->DXP[IM1])*(p->DZP[KP]+p->DZP[KM1]))*p->sigz[IJ]

						 + (CPORNH*UHin[IJK])/(alpha*p->dt)
                            
                            
                            + visc_IPh*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(UH[Ip1JKp1] - UH[Im1JKp1] - UH[Ip1JKm1] + UH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZP[KP]+p->DZP[KM1]));
        }
        
        if(p->wet[IJ]==0 || p->flag4[IJK]<0 || p->DF[IJK]<0)
        {
        d->M.p[n]  =  1.0;

        d->M.n[n] = 0.0;
        d->M.s[n] = 0.0;

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
