/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/
#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h" 
#include<algorithm>

void rheology_f::u_source(lexer *p, fdm *a)
{
    // Force base F = A*tau
    double tau;
 
    count=0;
    if(p->W110==2 || p->W110==3)
    ULOOP
    {
        pressurePhi(p,a,1,0,0);
        
        tau = yield_stress(p,a);

         
        if(p->W110==3)
        tau += ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->u(i,j,k))>1.0e-20?(a->u(i,j,k)/fabs(a->u(i,j,k))):0.0;

        H = heaviside(phival);
        
        a->rhsvec.V[count] -= H*f*(tau/(p->DZN[KP]*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))));

        ++count;
    }
    
    // Gradient based
    double dudx,dudy,dudz;
    double fxx,fxy,fxz;
    double dpdx,dpdy,dpdz;
  
    count=0;
    if(p->W110==4 || p->W110==5)
    ULOOP
    {
        pressurePhiGradient(p,a,1,0,0);

        // Yield Stress
        yieldStressGradient(p,a,1,0,0);
         
        if(p->W110==5)
        {
            tau01 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
            tau02 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i+1,j,k);
        }
        
        f = fabs(a->u(i,j,k))>1.0e-20?(a->u(i,j,k)/fabs(a->u(i,j,k))):0.0;

        H = heaviside(phival);
        
        a->rhsvec.V[count] += H*((tau02-tau01)/(p->DXN[IP]*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)))); // *f ?

        ++count;
	}
    
    
    count=0;
    if(p->W110==6 || p->W110==7)
    ULOOP
    {
        // PFo: Adding direction to the yield stress gradient below. Should be *f for each term, where f = (pu_idx_j)/fabs(pu_idx_j) 
        // f(pwdx)*tau_x-term, f(pwdy)*tau_y-term, f(pwdz)*tau_z-term - PFo not sure how to write this
        // Try with u_x first instead: f = (a->u(i,j,k)/fabs(a->u(i,j,k)))

        // Velocity gradients at location of u(i,j,k):
        dudx = (a->u(i+1,j,k) - a->u(i-1,j,k))/(p->DXN[IP]+p->DXN[IP1]);
        dudy = (a->u(i,j+1,k) - a->u(i,j-1,k))/(p->DYN[JP]+p->DYN[JP1]);
        dudz = (a->u(i,j,k+1) - a->u(i,j,k-1))/(p->DZN[KP]+p->DZN[KP1]);
        
        // Sign of velocity gradients, which determines the sign of the deviatoric stress gradient contributions to the source term in x-direction
        fxx = fabs(dudx)>1.0e-20?(dudx/fabs(dudx)):0.0; // dudx
        fxy = fabs(dudy)>1.0e-20?(dudy/fabs(dudy)):0.0; // dudy
        fxz = fabs(dudz)>1.0e-20?(dudz/fabs(dudz)):0.0; // dudz

        // Pressure gradients at location of u(i,j,k):
        dpdx = (a->press(i+1,j,k) - a->press(i,j,k))/(p->DXM);
        dpdy = (0.5*(a->press(i,j+1,k)+a->press(i+1,j+1,k)) - 0.5*(a->press(i,j-1,k)+a->press(i+1,j-1,k)))/(2.0*p->DXM);
        dpdz = (0.5*(a->press(i,j,k+1)+a->press(i+1,j,k+1)) - 0.5*(a->press(i,j,k-1)+a->press(i+1,j,k-1)))/(2.0*p->DXM);

        phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        
        H = heaviside(phival);
                 
        a->rhsvec.V[count] += H*tanphi*(fxx*dpdx + fxy*dpdy + fxz*dpdz)/(0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)));

        ++count;
    }
}

void rheology_f::v_source(lexer *p, fdm *a)
{
    double dvdx,dvdy,dvdz;
    double fyx,fyy,fyz;
    double dpdx,dpdy,dpdz;
    
    count=0;
    if(p->W110==2 || p->W110==3)
    VLOOP
    {
        pressurePhi(p,a,0,1,0);
        
        tau0 = yield_stress(p,a);
        
        if(p->W110==5)
            tau0 += ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->v(i,j,k))>1.0e-20?(a->v(i,j,k)/fabs(a->v(i,j,k))):0.0;
        
        H = heaviside(phival);
        
        a->rhsvec.V[count] -= H*f*(tau0/(p->DYN[JP]*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))));
        
        ++count;
    }

        
    count=0;
    if(p->W110==6 || p->W110==7)
    VLOOP
    {
        // Velocity gradients at location of v(i,j,k):
        dvdx = (a->v(i+1,j,k) - a->v(i-1,j,k))/(p->DXN[IP]+p->DXN[IP1]);
        dvdy = (a->v(i,j+1,k) - a->v(i,j-1,k))/(p->DYN[JP]+p->DYN[JP1]);
        dvdz = (a->v(i,j,k+1) - a->v(i,j,k-1))/(p->DZN[KP]+p->DZN[KP1]);
        
        // Sign of velocity gradients, which determines the sign of the deviatoric stress gradient contributions to the source term in y-direction
        fyx = fabs(dvdx)>1.0e-20?(dvdx/fabs(dvdx)):0.0; // dvdx
        fyy = fabs(dvdy)>1.0e-20?(dvdy/fabs(dvdy)):0.0; // dvdy
        fyz = fabs(dvdz)>1.0e-20?(dvdz/fabs(dvdz)):0.0; // dvdz

        // Pressure gradients at location of v(i,j,k):
        dpdx = (0.5*(a->press(i+1,j,k)+a->press(i+1,j+1,k)) - 0.5*(a->press(i-1,j,k)+a->press(i-1,j+1,k)))/(2.0*p->DXM);
        dpdy = (a->press(i,j+1,k) - a->press(i,j,k))/(p->DXM);
        dpdz = (0.5*(a->press(i,j,k+1)+a->press(i,j+1,k+1)) - 0.5*(a->press(i,j,k-1)+a->press(i,j+1,k-1)))/(2.0*p->DXM);

        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        H = heaviside(phival);
                 
        a->rhsvec.V[count] += H*tanphi*(fyx*dpdx + fyy*dpdy + fyz*dpdz)/(0.5*(a->ro(i,j,k)+a->ro(i,j+1,k)));

        ++count;
    }
}

void rheology_f::w_source(lexer *p, fdm *a)
{   
    count=0;
    if(p->W110==2 || p->W110==3)
    WLOOP
    {
        pressurePhi(p,a,0,0,1);
        
        tau0 = yield_stress(p,a);
            
        if(p->W110==3)
            tau0 += ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->w(i,j,k))>1.0e-20?(a->w(i,j,k)/fabs(a->w(i,j,k))):0.0;

        H = heaviside(phival);
        
        a->rhsvec.V[count] -= H*f*(tau0/(0.5*(p->DZN[KP]*a->ro(i,j,k)+a->ro(i,j,k+1))));

        ++count;
	}
    
    // Gradient Based
    double dwdx,dwdy,dwdz;
    double fzx,fzy,fzz;
    double dpdx,dpdy,dpdz;
    
    count=0;
    if(p->W110==4 || p->W110==5)
    WLOOP
    {
        pressurePhiGradient(p,a,0,0,1);

        // Yield Stress
        yieldStressGradient(p,a,0,0,1);
         
        if(p->W110==5)
        {
            tau01 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
            tau02 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k+1);
        }
        
        f = fabs(a->w(i,j,k))>1.0e-20?(a->w(i,j,k)/fabs(a->w(i,j,k))):0.0;

        H = heaviside(phival);
        
        a->rhsvec.V[count] += H*((tau02-tau01)/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1)))); // PFo: Missing "*f" ?

        ++count;
	}
    
    
    count=0;
    if(p->W110==6 || p->W110==7)
    WLOOP
    {
        // Velocity gradients at location of w(i,j,k):
        dwdx = (a->w(i+1,j,k) - a->w(i-1,j,k))/(p->DXN[IP]+p->DXN[IP1]);
        dwdy = (a->w(i,j+1,k) - a->w(i,j-1,k))/(p->DYN[JP]+p->DYN[JP1]);
        dwdz = (a->w(i,j,k+1) - a->w(i,j,k-1))/(p->DZN[KP]+p->DZN[KP1]);
    
        
        // Sign of velocity gradients, which determines the sign of the deviatoric stress gradient contributions to the source term in w-direction
        fzx = fabs(dwdx)>1.0e-20?(dwdx/fabs(dwdx)):0.0; // dwdx
        fzy = fabs(dwdy)>1.0e-20?(dwdy/fabs(dwdy)):0.0; // dwdy
        fzz = fabs(dwdz)>1.0e-20?(dwdz/fabs(dwdz)):0.0; // dwdz
       

        // Pressure gradients at location of w(i,j,k):
        dpdx = (0.5*(a->press(i+1,j,k)+a->press(i+1,j,k+1)) - 0.5*(a->press(i-1,j,k)+a->press(i-1,j,k+1)))/(2.0*p->DXM);
        dpdy = (0.5*(a->press(i,j+1,k)+a->press(i,j+1,k+1)) - 0.5*(a->press(i,j-1,k)+a->press(i,j-1,k+1)))/(2.0*p->DXM);
        dpdz = (a->press(i,j,k+1) - a->press(i,j,k))/(p->DXM);

        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        H = heaviside(phival);
                         
        a->rhsvec.V[count] += H*tanphi*(fzx*dpdx + fzy*dpdy + fzz*dpdz)/(0.5*(a->ro(i,j,k)+a->ro(i,j,k+1)));

        ++count;
    }
}

