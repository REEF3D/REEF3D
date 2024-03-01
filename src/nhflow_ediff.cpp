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

void nhflow_ediff::diff_u(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *UHdiff, double *UHin, double *UH, double *VH, double *WH, double alpha)
{
    /*
	starttime=pgc->timer();
    
    LOOP
    UHdiff[IJK] = UHin[IJK];
    
    pgc->start4V(p,UHdiff,gcval_uh);
    
    LOOP
	{
    visc = d->VISC[IJK];
    
    d->F(i,j,k) += 2.0*visc*((UH[Im1JK] - 2.0*UH[IJK] + UH[Ip1JK])/(p->XP[IM1] - 2.0*p->XP[IP] + p->XP[IP1]))
                       
        +   visc*((VH[IJm1K] - 2.0*VH[IJK] + VH[IJp1K])/(p->YP[JM1] - 2.0*p->YP[JP] + p->YP[JP1]))

        +   visc*((WH[IJKm1] - 2.0*WH[IJK] + ¨WH[IJKp1])/(p->ZP[KM1] - 2.0*p->ZP[KP] + p->ZP[KP1]))


        + visc*(VH[Ip1Jp1K]-VH[Im1Jp1K] - VH[Ip1Jm1K]+VH[Im1Jm1K])/(2.0*p->DXP[IP]*p->DYN[JP])
        
        + visc*(WH[Ip1JKp1]-WH[Im1JKp1] - WH[Ip1JKm1]+WH[Im1JKm1])/(2.0*p->DXP[IP]*p->DZN[JP]);
		
	}

    

    

    n=0;
    LOOP
	{
        if(p->wet[IJ]==1  && d->breaking(i,j)==0)
        {
            visc = d->VISC[IJK];
            
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            
            d->M.p[n]  =  2.0*visc/(p->DXP[IP]*p->DXN[IP])
                        + 2.0*visc/(p->DXP[IM1]*p->DXN[IP])
                        
                        + visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                        + visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir
                        
                        + (visc*sigxyz2)/(p->DZP[KM1]*p->DZN[KP])
                        + (visc*sigxyz2)/(p->DZP[KM1]*p->DZN[KM1]);


            d->M.n[n] = -2.0*visc/(p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -2.0*visc/(p->DXP[IM1]*p->DXN[IP]);

            d->M.w[n] = -visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
            d->M.e[n] = -visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

            d->M.t[n] = -(visc*sigxyz2)/(p->DZP[KP]*p->DZN[KP])     
                        - p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
                        
            d->M.b[n] = -(visc*sigxyz2)/(p->DZP[KM1]*p->DZN[KP]) 
                        + p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
            
            
            d->rhsvec.V[n] = visc*((VH[Ip1Jp1K]-VH[Im1Jp1K]) - (VH[Ip1Jm1K]-VH[Im1Jm1K]))/((p->DXP[IP]+p->DXP[IM1])*(p->DYN[JP]+p->DYN[JM1]))
						 +  visc*((WH[Ip1JKp1]-WH[Im1JKp1]) - (WH[Ip1JKm1]-WH[Im1JKm1]))/((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))

						 + (CPORNH*UHin[IJK])/(alpha*p->dt)
                            
                            
                            + visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(UH[Ip1JKp1] - UH[Im1JKp1] - UH[Ip1JKm1] + UH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
                            + visc*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(UH[IJp1Kp1] - UH[IJm1Kp1] - UH[IJp1Km1] + UH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        }
        
*/
}

void nhflow_ediff::diff_v(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *VHdiff, double *VHin, double *UH, double *VH, double *WH, double alpha)
{
}

void nhflow_ediff::diff_w(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *WHdiff, double *WHin, double *UH, double *VH, double *WH, double alpha)
{
}

void nhflow_ediff::diff_scalar(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *UHdiff, double *UH_in, double *UH, double *VH, double *WH, double alpha)
{
}