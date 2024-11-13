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

#include"nhflow_poisson_phs.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"

nhflow_poisson_phs::nhflow_poisson_phs(lexer *p) 
{
}

nhflow_poisson_phs::~nhflow_poisson_phs()
{
}

void nhflow_poisson_phs::start(lexer* p, fdm_nhf *d, double *P)
{	
    double sigxyz2;
   
	n=0;
    LOOP
	{
        WETDRYDEEP
        //if(0.5*(d->FB[IJK]+d->FB[IJKm1])>=0.0 && 0.5*(d->SOLID[IJK]+d->SOLID[IJKm1]>=0.0))
        {
            sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            
            d->M.p[n]  =  (CPORNH*PORVALNH)/(p->W1*p->DXP[IP]*p->DXN[IP])
                        + (CPORNHm*PORVALNHm)/(p->W1*p->DXP[IM1]*p->DXN[IP])
                        
                        + (CPORNH*PORVALNH)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir
                        + (CPORNHm*PORVALNHm)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir
                        
                        + (sigxyz2*CPORNH*PORVALNH)/(p->W1*p->DZP[KM1]*p->DZN[KP])
                        + (sigxyz2*CPORNHm*PORVALNHm)/(p->W1*p->DZP[KM1]*p->DZN[KM1]);


            d->M.n[n] = -(CPORNH*PORVALNH)/(p->W1*p->DXP[IP]*p->DXN[IP]);
            d->M.s[n] = -(CPORNHm*PORVALNHm)/(p->W1*p->DXP[IM1]*p->DXN[IP]);

            d->M.w[n] = -(CPORNH*PORVALNH)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir;
            d->M.e[n] = -(CPORNHm*PORVALNHm)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir;

            d->M.t[n] = -(sigxyz2*CPORNH*PORVALNH)/(p->W1*p->DZP[KM1]*p->DZN[KP])     
                        - CPORNH*PORVALNH*p->sigxx[FIJK]/(p->W1*(p->DZN[KP]+p->DZN[KM1]));
                        
            d->M.b[n] = -(sigxyz2*CPORNHm*PORVALNHm)/(p->W1*p->DZP[KM1]*p->DZN[KM1]) 
                        + CPORNH*PORVALNH*p->sigxx[FIJK]/(p->W1*(p->DZN[KP]+p->DZN[KM1]));
            
            
            d->rhsvec.V[n] +=  CPORNH*PORVALNH*2.0*p->sigx[FIJK]*(P[FIp1JKp1] - P[FIm1JKp1] - P[FIp1JKm1] + P[FIm1JKm1])
                            /(p->W1*(p->DXP[IP]+p->DXP[IM1])*(p->DZN[KP]+p->DZN[KM1]))
                        
                            + CPORNH*PORVALNH*2.0*p->sigy[FIJK]*(P[FIJp1Kp1] - P[FIJm1Kp1] - P[FIJp1Km1] + P[FIJm1Km1])
                            /(p->W1*(p->DYP[JP]+p->DYP[JM1])*(p->DZN[KP]+p->DZN[KM1]))*p->y_dir;
        }
        
        if(p->wet[IJ]==0 || p->deep[IJ]==0 || p->flag7[FIJK]<0)//  || 0.5*(d->FB[IJK]+d->FB[IJKm1])<0.0 || 0.5*(d->SOLID[IJK]+d->SOLID[IJKm1]<0.0))
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
        WETDRYDEEP
        //if(0.5*(d->FB[IJK]+d->FB[IJKm1])>=0.0 && 0.5*(d->SOLID[IJK]+d->SOLID[IJKm1]>=0.0))
        {
            if(p->flag7[FIm1JK]<0 || p->wet[Im1J]==0 || p->deep[Im1J]==0)//  || 0.5*(d->FB[Im1JK]+d->FB[Im1JKm1])<0.0 || 0.5*(d->SOLID[Im1JK]+d->SOLID[Im1JKm1]<0.0))
            {
            d->rhsvec.V[n] -= d->M.s[n]*P[FIJK];
            d->M.s[n] = 0.0;
            }
            
            if(p->flag7[FIp1JK]<0 || p->wet[Ip1J]==0 || p->deep[Ip1J]==0)//  || 0.5*(d->FB[Ip1JK]+d->FB[Ip1JKm1])<0.0 || 0.5*(d->SOLID[Ip1JK]+d->SOLID[Ip1JKm1]<0.0))
            {
            d->rhsvec.V[n] -= d->M.n[n]*P[FIJK];
            d->M.n[n] = 0.0;
            }
            
            if(p->flag7[FIJm1K]<0 || p->wet[IJm1]==0 || p->deep[IJm1]==0)//  || 0.5*(d->FB[IJm1K]+d->FB[IJm1Km1])<0.0 || 0.5*(d->SOLID[IJm1K]+d->SOLID[IJm1Km1]<0.0))
            {
            d->rhsvec.V[n] -= d->M.e[n]*P[FIJK]*p->y_dir;
            d->M.e[n] = 0.0;
            }
            
            if(p->flag7[FIJp1K]<0 || p->wet[IJp1]==0 || p->deep[IJp1]==0)//  || 0.5*(d->FB[IJp1K]+d->FB[IJp1Km1])<0.0 || 0.5*(d->SOLID[IJp1K]+d->SOLID[IJp1Km1]<0.0))
            {
            d->rhsvec.V[n] -= d->M.w[n]*P[FIJK]*p->y_dir;
            d->M.w[n] = 0.0;
            }
            /*
            // TOP
            if(0.5*(d->FB[IJK]+d->FB[IJKp1])<0.0)// || 0.5*(d->SOLID[IJK]+d->SOLID[IJKp1]<0.0))
            {
            d->rhsvec.V[n] -= d->M.t[n]*P[FIJK]*p->y_dir;
            d->M.t[n] = 0.0;
            }
            
            // BOTTOM
            if(0.5*(d->FB[IJKm1]+d->FB[IJKm2])<0.0)// || 0.5*(d->SOLID[IJKm1]+d->SOLID[IJKm2]<0.0))
            {
            d->rhsvec.V[n] -= d->M.b[n]*P[FIJK]*p->y_dir;
            d->M.b[n] = 0.0;
            }*/
            
            
            // BED
            if(p->flag7[FIJKm1]<0)
            {
            d->rhsvec.V[n] -= d->M.b[n]*P[FIJK];
            d->M.b[n] = 0.0;
            }
            
            // FSFBC
            if(p->flag7[FIJKp2]<0 && p->flag7[FIJKp1]>0)
            {
            d->rhsvec.V[n] -= 0.0; // fsf: p=0
            d->M.t[n] = 0.0;
            }
            
            
            
            

  
        }
	++n;
	}
}


