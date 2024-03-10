/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"nhflow_ediff.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"solver.h"

nhflow_ediff::nhflow_ediff(lexer* p)
{
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_uh=20;
	gcval_vh=21;
	gcval_wh=22;
}

nhflow_ediff::~nhflow_ediff()
{
}

void nhflow_ediff::diff_u(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *UHdiff, double *UHin, double *UH, double *VH, double *WH, slice &WL, double alpha)
{
    
	starttime=pgc->timer();
    
    LOOP
    UHdiff[IJK] = UHin[IJK];
    
    pgc->start4V(p,UHdiff,gcval_uh);
    
    LOOP
	{
    visc = d->VISC[IJK];
    
    sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
    

    d->F[IJK] += 2.0*visc*((UH[Ip1JK]-UH[IJK])/p->DXP[IP] - (UH[IJK]-UH[Im1JK])/p->DXP[IP])/p->DXN[IP]
    
                   + visc*((UH[IJp1K]-UH[IJK])/p->DYP[JP] - (UH[IJK]-UH[IJm1K])/p->DYP[JP])/p->DYN[JP]
                   
                   + visc*sigxyz2*((UH[IJKp1]-UH[IJK])/p->DZP[KP] - (UH[IJK]-UH[IJKm1])/p->DZP[KP])/p->DZN[KP]
    
    
        
         + visc*((VH[Ip1Jp1K]-VH[Im1Jp1K]) - (VH[Ip1Jm1K]-VH[Im1Jm1K]))/((p->DXP[IP]+p->DXP[IM1])*(p->DYN[JP]+p->DYN[JM1]))
         
         + visc*((WH[Ip1JKp1]-WH[Im1JKp1]) - (WH[Ip1JKm1]-WH[Im1JKm1]))/((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
        
        + visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(UH[Ip1JKp1] - UH[Im1JKp1] - UH[Ip1JKm1] + UH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
        + visc*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(UH[IJp1Kp1] - UH[IJm1Kp1] - UH[IJp1Km1] + UH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
		
	}
}

void nhflow_ediff::diff_v(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *VHdiff, double *VHin, double *UH, double *VH, double *WH, slice &WL, double alpha)
{
    LOOP
    VHdiff[IJK] = VHin[IJK];
    
    pgc->start4V(p,VHdiff,gcval_vh);
    
    
    LOOP
	{
    visc = d->VISC[IJK];
    
    sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
    

    d->G[IJK] +=    visc*((VH[Ip1JK]-VH[IJK])/p->DXP[IP] - (VH[IJK]-VH[Im1JK])/p->DXP[IP])/p->DXN[IP]
    
                   + 2.0*visc*((VH[IJp1K]-VH[IJK])/p->DYP[JP] - (VH[IJK]-VH[IJm1K])/p->DYP[JP])/p->DYN[JP]
                   
                   + visc*sigxyz2*((VH[IJKp1]-VH[IJK])/p->DZP[KP] - (VH[IJK]-VH[IJKm1])/p->DZP[KP])/p->DZN[KP]
    
    
        
         + visc*((UH[Ip1Jp1K]-UH[Ip1Jm1K]) - (UH[Im1Jp1K]-UH[Im1Jm1K]))/((p->DYN[JP]+p->DYN[JM1])*(p->DXP[IP]+p->DXP[IM1]))
         +  visc*((WH[IJp1Kp1]-WH[IJm1Kp1]) - (WH[IJp1Km1]-WH[IJm1Km1]))/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))
        
        + visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(UH[Ip1JKp1] - UH[Im1JKp1] - UH[Ip1JKm1] + UH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
        + visc*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(UH[IJp1Kp1] - UH[IJm1Kp1] - UH[IJp1Km1] + UH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
		
	}
}

void nhflow_ediff::diff_w(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *WHdiff, double *WHin, double *UH, double *VH, double *WH, slice &WL, double alpha)
{
    LOOP
    WHdiff[IJK] = WHin[IJK];
    
    pgc->start4V(p,WHdiff,gcval_wh);
    
    
    LOOP
	{
    visc = d->VISC[IJK];
    
    sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
    

    d->H[IJK] +=    visc*((WH[Ip1JK]-WH[IJK])/p->DXP[IP] - (WH[IJK]-WH[Im1JK])/p->DXP[IP])/p->DXN[IP]
    
                   + visc*((WH[IJp1K]-WH[IJK])/p->DYP[JP] - (WH[IJK]-WH[IJm1K])/p->DYP[JP])/p->DYN[JP]
                   
                   + 2.0*visc*sigxyz2*((WH[IJKp1]-WH[IJK])/p->DZP[KP] - (WH[IJK]-WH[IJKm1])/p->DZP[KP])/p->DZN[KP]
    
    
        
         + visc*((UH[Ip1JKp1]-UH[Ip1JKm1]) - (UH[Im1JKp1]-UH[Im1JKm1]))/((p->DZN[KP]+p->DZN[KM1])*(p->DXP[IP]+p->DXP[IM1]))
         +  visc*((VH[IJp1Kp1]-VH[IJp1Km1]) - (VH[IJp1Kp1]-VH[IJm1Km1]))/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))
        
        + visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(WH[Ip1JKp1] - WH[Im1JKp1] - WH[Ip1JKm1] + WH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
        + visc*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(WH[IJp1Kp1] - WH[IJm1Kp1] - WH[IJp1Km1] + WH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
		
	}
}

void nhflow_ediff::diff_scalar(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *UHdiff, double *UH_in, double *UH, double *VH, double *WH, double alpha)
{
}