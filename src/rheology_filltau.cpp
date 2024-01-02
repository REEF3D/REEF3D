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

#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h" 


void rheology_f::filltau(lexer *p, fdm *a, ghostcell *pgc)
{
    LOOP
    {
        phival = a->phi(i,j,k);
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));    
        

        if(p->W111==1)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        
        if(p->W111==2)
        pval=0.5*((a->press(i,j,k)-p->pressgage)+(a->press(i+1,j,k)-p->pressgage));
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            pval=0.5*((a->press(i,j,k)-p->pressgage)+(a->press(i+1,j,k)-p->pressgage));
            
            if(phival>=p->W112*p->DXM)
            pval=phival*0.5*((a->press(i,j,k)-p->pressgage)+(a->press(i+1,j,k)-p->pressgage))*fabs(p->W22);
        }
        
        if(p->W111==4)
        pval=0.25*((a->press(i,j,k)-p->pressgage)+(a->press(i+1,j,k)-p->pressgage)) + 0.5*phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        
        
        if(p->W111==5)
        {
            if(p->count>=10)
            pval=0.5*((a->press(i,j,k)-p->pressgage)+(a->press(i+1,j,k)-p->pressgage));
            
            if(p->count<10)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        }
        
        // Yield Stress
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        tau0 = tanphi*pval + p->W102_c;
        
        if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
        tau0 = (tanphi*pval + p->W102_c)*(1.0-exp(-p->W103*gamma));
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // rho_water = 1000.0, new input?
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_p is new input W 104 
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_u also use new input W 104

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==7)
        tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);   
        
        tau_x(i,j,k) = tau0;
        tau_y(i,j,k) = tau0;
        tau_z(i,j,k) = tau0;
    }

    pgc->start4(p,tau_x,1);
    pgc->start4(p,tau_y,1);
    pgc->start4(p,tau_z,1);
}
