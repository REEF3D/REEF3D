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

#include"vrans_nhflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void vrans_nhflow_f::force_calc(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, int val)
{
    Fx=Fy=Fz=0.0;
    
    FLOOP
    P[FIJK] = (p->wd + d->eta(i,j) - p->ZSN[FIJK])*p->W1*fabs(p->W22) + d->P[FIJK];
    
    // VRANS force
    if(p->B200>0)
    LOOP
	{
        H = Hporface(p,d,0,0,0);  
        
        porval   = d->POR[IJK];
        partval  = d->PORPART[IJK];
        alphaval = APOR[IJK];
        betaval  = BPOR[IJK];
        viscval  = d->VISC[IJK];
		
        
        Aporval = Apor(porval,partval,alphaval,viscval);
        Bporval = Bpor(porval,partval,betaval);
        
        dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*p->WL[IJ];
            

        Fx += dV*H*(Aporval*d->U[IJK] + Bporval*d->U[IJK]*fabs(d->U[IJK]) 
        
                + d->RO[IJK]*p->B260*(d->U[IJK]-UN[IJK])/(alpha*p->dt*PORVALNH)*cmfac
                
                - (1.0-PORVALNH)*(((0.5*(P[FIp1JKp1]+P[FIp1JK])-0.5*(P[FIm1JKp1]+P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1]))
                
                + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((P[FIJKp1]-P[FIJK])/p->DZN[KP])));
        
        
        Fy += dV*H*(Aporval*d->V[IJK] + Bporval*d->V[IJK]*fabs(d->V[IJK])
        
                //+ d->RO[IJK]*p->B260*(d->V[IJK]-VN[IJK])/(alpha*p->dt*PORVALNH)*cmfac
                
                + (1.0-PORVALNH)*(((0.5*(P[FIJp1Kp1]+P[FIJp1K])-0.5*(P[FIJm1Kp1]+P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1]))
                
                + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((P[FIJKp1]-P[FIJK])/p->DZN[KP])));
        
        
        Fz += dV*H*(Aporval*d->W[IJK] + Bporval*d->W[IJK]*fabs(d->W[IJK])
        
                //+ d->RO[IJK]*p->B260*(d->W[IJK]-WN[IJK])/(alpha*p->dt*PORVALNH)*cmfac
                
                - (1.0-PORVALNH)*((P[FIJKp1]-P[FIJK])/(p->DZN[KP])));
	}
    
    Fx = pgc->globalsum(Fx);
    Fy = pgc->globalsum(Fy);
    Fz = pgc->globalsum(Fz);
    
    LOOP
    {
    UN[IJK] = d->U[IJK];
    VN[IJK] = d->V[IJK];
    WN[IJK] = d->W[IJK];
    cmfac = 1.0;
    }
    
    // print
    if(val==1)
    print_force(p,d,pgc);

}



