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
    double rho,cmval,nfac,pfac;
    double Du,Dv,Dw,dpdx,dpdy,dpdz;
    
    Fx=Fy=Fz=0.0;
    
    FLOOP
    P[FIJK] = (p->wd + d->eta(i,j) - p->ZSN[FIJK])*p->W1*fabs(p->W22) + d->P[FIJK];
    
    // VRANS force on the static structures, consistent with the NHFLOW momentum equation of the Darcy velocity U:
    // F = sum dV ( -grad p - rho dU/dt ), written term-wise (D = A U + B |U| U, n = n_s of the static structure):
    //
    // B 265 0:  (1 + c) dU/dt = -grad p/rho - D,         c = gamma (1-n)/n^2
    //           the full pressure gradient acts on U, the reaction on the structure is
    //           F = sum dV rho ( D + c dU/dt ),  no (1-n) grad p term (it is contained in c dU/dt)
    //
    // B 265 1:  (1 + c_A)/n dU/dt = -grad p/rho - D,     c_A = gamma (1-n)/n   (Jensen et al. 2014)
    //           F = sum dV [ n rho D + rho c_A dU/dt - (1-n) grad p ]
    //
    // all terms vanish with (1 - n) outside the structure.
    // Cells inside direct forcing solids (DF < 0, e.g. a pile inside a porous cover) are skipped:
    // the pressure on their surface is integrated by the P 81 force, their volume would be counted twice.
    // grad p in sigma coordinates: dp/dx = dp/dx|_sig + sigx dp/dsig,  dp/dz = sigz dp/dsig
    if(p->B200>0)
    LOOP
	{
        if(p->DF[IJK]<0)
        continue;
        
        porous_coeff(p,d);
        
        if(H<1.0e-12)
        continue;
        
        dV = p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*p->WL[IJ];
        
        rho = d->RO[IJK];
        
        if(p->B265==1)
        {
        cmval = p->B260*(1.0-porval)/porval;
        nfac  = porval;
        pfac  = 1.0-porval;
        }
        
        else
        {
        cmval = p->B260*(1.0-porval)/(porval*porval);
        nfac  = 1.0;
        pfac  = 0.0;
        }
        
        Du = Aporval*d->U[IJK] + Bporval*d->U[IJK]*fabs(d->U[IJK]);
        Dv = Aporval*d->V[IJK] + Bporval*d->V[IJK]*fabs(d->V[IJK]);
        Dw = Aporval*d->W[IJK] + Bporval*d->W[IJK]*fabs(d->W[IJK]);
        
        dpdx = (0.5*(P[FIp1JKp1]+P[FIp1JK])-0.5*(P[FIm1JKp1]+P[FIm1JK]))/(p->DXP[IP]+p->DXP[IM1])
        
             + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*((P[FIJKp1]-P[FIJK])/p->DZN[KP]);
        
        dpdy = (0.5*(P[FIJp1Kp1]+P[FIJp1K])-0.5*(P[FIJm1Kp1]+P[FIJm1K]))/(p->DYP[JP]+p->DYP[JM1])
        
             + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*((P[FIJKp1]-P[FIJK])/p->DZN[KP]);
        
        dpdz = p->sigz[IJ]*(P[FIJKp1]-P[FIJK])/p->DZN[KP];
        
        
        Fx += dV*(nfac*rho*Du + rho*cmval*(d->U[IJK]-UN[IJK])/(alpha*p->dt)*cmfac - pfac*dpdx);
        
        Fy += dV*(nfac*rho*Dv + rho*cmval*(d->V[IJK]-VN[IJK])/(alpha*p->dt)*cmfac - pfac*dpdy);
        
        Fz += dV*(nfac*rho*Dw + rho*cmval*(d->W[IJK]-WN[IJK])/(alpha*p->dt)*cmfac - pfac*dpdz);
	}
    
    Fx = pgc->globalsum(Fx);
    Fy = pgc->globalsum(Fy);
    Fz = pgc->globalsum(Fz);
    
    LOOP
    {
    UN[IJK] = d->U[IJK];
    VN[IJK] = d->V[IJK];
    WN[IJK] = d->W[IJK];
    }
    
    // no time derivative before UN holds a velocity (first call)
    cmfac = 1.0;
    
    // print
    if(val==1)
    print_force(p,d,pgc);

}



