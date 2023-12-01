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

#include"nhflow_idiff.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"solver.h"

void nhflow_idiff::diff_w(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, double *WHdiff, double *WHin, double *UH, double *VH, double *WH, double alpha)
{
    /*
	starttime=pgc->timer();
    
    LOOP
    WHdiff[IJK] = WHin[IJK];
    
    pgc->start4V(p,WHdiff,gcval_wh);

    n=0;
    LOOP
	{
        if(p->wet[IJ]==1 && p->deep[IJ]==1 && d->breaking(i,j)==0)
        {
            visc = d->VISC[IJK];
            
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            
            d->M.p[n]  =  visc/(p->DXP[IP]*p->DXN[IP])
                        + visc/(p->DXP[IM1]*p->DXN[IP])
                        
                        + visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                        + visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir
                        
                        + 2.0*(visc*sigxyz2)/(p->DZP[KM1]*p->DZN[KP])
                        + 2.0*(visc*sigxyz2)/(p->DZP[KM1]*p->DZN[KM1]);


            d->M.n[n] = -visc/(p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -visc/(p->DXP[IM1]*p->DXN[IP]);

            d->M.w[n] = -visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
            d->M.e[n] = -visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

            d->M.t[n] = -2.0*(visc*sigxyz2)/(p->DZP[KM1]*p->DZN[KP])     
                        - 2.0*p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
                        
            d->M.b[n] = -2.0*(visc*sigxyz2)/(p->DZP[KM1]*p->DZN[KM1]) 
                        + 2.0*p->sigxx[FIJK]/((p->DZN[KP]+p->DZN[KM1]));
            
            
            d->rhsvec.V[n] = visc*((UH[Ip1JKp1]-UH[Ip1JKm1]) - (UH[Im1JKp1]-UH[Im1JKm1]))/((p->DZN[KP]+p->DZN[KM1])*(p->DXP[IP]+p->DXP[IM1]))
						 +  visc*((VH[IJp1Kp1]-VH[IJp1Km1]) - (VH[IJp1Kp1]-VH[IJm1Km1]))/((p->DYN[JP]+p->DYN[JM1])*(p->DZN[KP]+p->DZN[KM1]))

						 + (CPORNH*WHin[IJK])/(alpha*p->dt)
                            
                            
                            + visc*2.0*0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(WH[Ip1JKp1] - WH[Im1JKp1] - WH[Ip1JKm1] + WH[Im1JKm1])
                            /((p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
                            + visc*2.0*0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(WH[IJp1Kp1] - WH[IJm1Kp1] - WH[IJp1Km1] + WH[IJm1Km1])
                            /((p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        }
        
        if(p->wet[IJ]==0 || p->deep[IJ]==0 || p->flag4[IJK]<0 || d->breaking(i,j)==1)
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
    
    
    
    n=0;
	LOOP
	{
        if(p->wet[IJ]==1 && d->breaking(i,j)==0)
        {
            if(p->flag4[Im1JK]<0)
            {
            d->rhsvec.V[n] -= d->M.s[n]*WH[Im1JK];
            d->M.s[n] = 0.0;
            }
            
            if(p->flag4[Ip1JK]<0)
            {
            d->rhsvec.V[n] -= d->M.n[n]*WH[Ip1JK];
            d->M.n[n] = 0.0;
            }
            
            if(p->flag4[IJm1K]<0)
            {
            d->rhsvec.V[n] -= d->M.e[n]*WH[IJm1K]*p->y_dir;
            d->M.e[n] = 0.0;
            }
            
            if(p->flag4[IJp1K]<0)
            {
            d->rhsvec.V[n] -= d->M.w[n]*WH[IJp1K]*p->y_dir;
            d->M.w[n] = 0.0;
            }
            
            if(p->flag4[IJKm1]<0)
            {
            d->rhsvec.V[n] -= d->M.b[n]*WH[IJKm1];
            d->M.b[n] = 0.0;
            }
            
            if(p->flag4[IJKp1]<0 && p->flag4[IJKp1]>0)
            {
            d->rhsvec.V[n] -= d->M.t[n]*WH[IJKp1];
            d->M.t[n] = 0.0;
            }
  
        }
	++n;
	}
	
    psolv->startV(p,pgc,WHdiff,d->rhsvec,d->M,4);
    
	
    pgc->start4V(p,WHdiff,gcval_wh);
    
    
	time=pgc->timer()-starttime;
	p->witer=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && (p->count%p->P12==0))
	cout<<"wdiffiter: "<<p->witer<<"  wdifftime: "<<setprecision(3)<<time<<endl;*/
}
