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


void rheology_f::u_source(lexer *p, fdm *a)
{
    double pval;
    
    // Force base F = A*tau
 
    count=0;
    if(p->W110==2 || p->W110==3)
    ULOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));    
        

        if(p->W111==1)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        
        if(p->W111==2)
        pval=0.5*(a->press(i,j,k)+a->press(i+1,j,k));
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            pval=0.5*(a->press(i,j,k)+a->press(i+1,j,k));
            
            if(phival>=p->W112*p->DXM)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        }
        
        if(p->W111==4)
        pval=0.25*(a->press(i,j,k)+a->press(i+1,j,k)) + 0.5*phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        
        
        if(p->W111==5)
        {
            if(p->count>=10)
            pval=0.5*(a->press(i,j,k)+a->press(i+1,j,k));
            
            if(p->count<10)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        }
        
        // Yield Stress
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        tau0=tanphi*pval + p->W102_c;
        
        if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
        tau0 = (tanphi*pval + p->W102_c);
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // rho_water = 1000.0, new input?
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // m_p is new input W 104 
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c);    // m_u also use new input W 104

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==3)
        tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->u(i,j,k))>1.0e-20?(a->u(i,j,k)/fabs(a->u(i,j,k))):0.0;
        
        a->rhsvec.V[count] -= H*(tau0/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))))*f;

        ++count;
    
	}
    
    
    
    // Gradient based
    
    double pval1,pval2;
    double tau01,tau02;
    double dudx,dudy,dudz;
    double fxx,fxy,fxz;
    double dpdx,dpdy,dpdz;
  
    
    count=0;
    if(p->W110==4 || p->W110==5)
    ULOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 
       
        

        if(p->W111==1)
        {
        pval1 = phival*a->ro(i,j,k)*fabs(p->W22);
        pval2 = phival*a->ro(i+1,j,k)*fabs(p->W22);
        }
        
        if(p->W111==2)
        {
        pval1 = a->press(i,j,k);
        pval2 = a->press(i+1,j,k);
        }
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            {
            pval1 = a->press(i,j,k);
            pval2 = a->press(i+1,j,k);
            }
            
            if(phival>=p->W112*p->DXM)
            {
            pval1 = phival*a->ro(i,j,k)*fabs(p->W22);
            pval2 = phival*a->ro(i+1,j,k)*fabs(p->W22);
            }
        }


        
        // Yield Stress
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        { 
        tau01 = tanphi*pval1 + p->W102_c;
        tau02 = tanphi*pval2 + p->W102_c;
        }
        
        if(p->W101==2)  // HB-C dry sand, without MAX 
        {
        tau01 = (tanphi*pval1 + p->W102_c);
        tau02 = (tanphi*pval2 + p->W102_c);
        }
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        {
        tau01 = MAX(0.0,tanphi*pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    
        tau02 = MAX(0.0,tanphi*pval2*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k) + p->W102_c);   
        }
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        {
        tau01 = MAX(0.0,tanphi*pval1*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);  
        tau02 = MAX(0.0,tanphi*pval2*exp(-p->W104*gamma)*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k) + p->W102_c);   
        }
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        {
        tau01 = MAX(0.0,tanphi*MAX(0.0,pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c); 
        tau02 = MAX(0.0,tanphi*MAX(0.0,pval2*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k)-p->W104*gamma) + p->W102_c); 
        }   

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==5)
        {
        tau01 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        tau02 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i+1,j,k);
        }
        
        f = fabs(a->u(i,j,k))>1.0e-20?(a->u(i,j,k)/fabs(a->u(i,j,k))):0.0;
        
        a->rhsvec.V[count] += H*((tau02-tau01)/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)))); // *f ?

        ++count;
	}
    
    
    count=0;
    if(p->W110==6 || p->W110==7)
    ULOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 
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
                 
        a->rhsvec.V[count] += H*tanphi*(fxx*dpdx + fxy*dpdy + fxz*dpdz)/(0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)));

        ++count;
    }
}

void rheology_f::v_source(lexer *p, fdm *a)
{
    double pval;
    double dvdx,dvdy,dvdz;
    double fyx,fyy,fyz;
    double dpdx,dpdy,dpdz;
//    double pvalx1,pvalx2,pvaly1,pvaly2,pvalz1,pvalz2;
//    double pvaldxx1,pvaldxx2,pvaldxy1,pvaldxy2,pvaldxz1,pvaldxz2;
//    double pvaldyx1,pvaldyx2,pvaldyy1,pvaldyy2,pvaldyz1,pvaldyz2;    
    
    count=0;
    if(p->W110>=2)
    VLOOP
	{
    phival = 0.5*(a->phi(i,j,k)+a->phi(i,j+1,k));
    
    if(phival>epsi)
    H=1.0;

    if(phival<-epsi)
    H=0.0;

    if(fabs(phival)<=epsi)
    H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));  
    

    if(p->W111==1)
    pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22);
    
    if(p->W111==2)
    pval=0.5*(a->press(i,j,k)+a->press(i,j+1,k));
    
    if(p->W111==3)
    {
        if(phival<p->W112*p->DXM)
        pval=0.5*(a->press(i,j,k)+a->press(i,j+1,k));
        
        if(phival>=p->W112*p->DXM)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22);
    }
    
    if(p->W111==4)
    pval=0.25*(a->press(i,j,k)+a->press(i,j+1,k)) + 0.5*phival*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22);
    
    if(p->W111==5)
    {
        if(p->count>=10)
        pval=0.5*(a->press(i,j,k)+a->press(i,j+1,k));
        
        if(p->count<10)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22);
    }
    
    // Yield Stress
    if(p->W101==0)
    tau0=p->W96;
    
    if(p->W101==1)  // HB-C dry sand
    tau0=tanphi*pval + p->W102_c;
    
    if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
    tau0 = (tanphi*pval + p->W102_c);
        
    if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
    tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // rho_water = 1000.0, new input?
        
    if(p->W101==4)  // HB-C shear rate generated excess pore pressure
    tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // m_p is new input W 104 
        
    if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
    tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c);    // m_u also use new input W 104

        if(p->count==0)
        tau0=p->W96;
    
    if(p->W110==5)
    tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
    
	f = fabs(a->v(i,j,k))>1.0e-20?(a->v(i,j,k)/fabs(a->v(i,j,k))):0.0;
    
    a->rhsvec.V[count] -= H*(tau0/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))))*f;
    
	++count;
	}
    
    
    count=0;
    if(p->W110==6 || p->W110==7)
    VLOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 

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
                 
        a->rhsvec.V[count] += H*tanphi*(fyx*dpdx + fyy*dpdy + fyz*dpdz)/(0.5*(a->ro(i,j,k)+a->ro(i,j+1,k)));

        ++count;
    }
}

void rheology_f::w_source(lexer *p, fdm *a)
{
    double pval;
    
    count=0;
    if(p->W110==2 || p->W110==3)
    WLOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));  
        

        if(p->W111==1)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22);
        
        if(p->W111==2)
        pval=0.5*(a->press(i,j,k)+a->press(i,j,k+1));
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            pval=0.5*(a->press(i,j,k)+a->press(i,j,k+1));
            
            if(phival>=p->W112*p->DXM)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22);
        }
        
        if(p->W111==4)
        pval=0.5*0.5*(a->press(i,j,k)+a->press(i,j,k+1)) + 0.5*phival*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22);
        
        if(p->W111==5)
        {
            if(p->count>=10)
            pval=0.5*(a->press(i,j,k)+a->press(i,j,k+1));
            
            if(p->count<10)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22);
        }
        
        
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        tau0=tanphi*pval + p->W102_c;
        
        if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
        tau0 = (tanphi*pval + p->W102_c);
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // rho_water = 1000.0, new input?
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // m_p is new input W 104 
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c);    // m_u also use new input W 104
        
            if(p->count==0)
            tau0=p->W96;
            
        if(p->W110==3)
        tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->w(i,j,k))>1.0e-20?(a->w(i,j,k)/fabs(a->w(i,j,k))):0.0;
        
        a->rhsvec.V[count] -= H*(tau0/(0.5*(p->DXM*a->ro(i,j,k)+a->ro(i,j,k+1))))*f;

        ++count;
	}
    
    
    
    // Gradient Based
    
    double pval1,pval2;
    double tau01,tau02;
    double dwdx,dwdy,dwdz;
    double fzx,fzy,fzz;
    double dpdx,dpdy,dpdz;
//    double pvalx1,pvalx2,pvaly1,pvaly2,pvalz1,pvalz2;
//    double pvaldxx1,pvaldxx2,pvaldxy1,pvaldxy2,pvaldxz1,pvaldxz2;
//    double pvaldyx1,pvaldyx2,pvaldyy1,pvaldyy2,pvaldyz1,pvaldyz2;
    
    count=0;
    if(p->W110==4 || p->W110==5)
    WLOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 
       
        

        if(p->W111==1)
        {
        pval1 = phival*a->ro(i,j,k)*fabs(p->W22);
        pval2 = phival*a->ro(i,j,k+1)*fabs(p->W22);
        }
        
        if(p->W111==2)
        {
        pval1 = a->press(i,j,k);
        pval2 = a->press(i,j,k+1);
        }
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            {
            pval1 = a->press(i,j,k);
            pval2 = a->press(i,j,k+1);
            }
            
            if(phival>=p->W112*p->DXM)
            {
            pval1 = phival*a->ro(i,j,k)*fabs(p->W22);
            pval2 = phival*a->ro(i,j,k+1)*fabs(p->W22);
            }
        }


        
        // Yield Stress
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        { 
        tau01 = tanphi*pval1 + p->W102_c;
        tau02 = tanphi*pval2 + p->W102_c;
        }
        
        if(p->W101==2)  // HB-C dry sand, without MAX 
        {
        tau01 = (tanphi*pval1 + p->W102_c);
        tau02 = (tanphi*pval2 + p->W102_c);
        }
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        {
        tau01 = MAX(0.0,tanphi*pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    
        tau02 = MAX(0.0,tanphi*pval2*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1) + p->W102_c);   
        }
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        {
        tau01 = MAX(0.0,tanphi*pval1*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);  
        tau02 = MAX(0.0,tanphi*pval2*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1) + p->W102_c);   
        }
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        {
        tau01 = MAX(0.0,tanphi*MAX(0.0,pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c); 
        tau02 = MAX(0.0,tanphi*MAX(0.0,pval2*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1)-p->W104*gamma) + p->W102_c); 
        }   

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==5)
        {
        tau01 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        tau02 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k+1);
        }
        
        f = fabs(a->w(i,j,k))>1.0e-20?(a->w(i,j,k)/fabs(a->w(i,j,k))):0.0;
        
        a->rhsvec.V[count] += H*((tau02-tau01)/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1)))); // PFo: Missing "*f" ?

        ++count;
	}
    
    
    count=0;
    if(p->W110==6 || p->W110==7)
    WLOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 

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
                         
        a->rhsvec.V[count] += H*tanphi*(fzx*dpdx + fzy*dpdy + fzz*dpdz)/(0.5*(a->ro(i,j,k)+a->ro(i,j,k+1)));

        ++count;
    }
    
}


/*


void rheology_f::u_source(lexer *p, fdm *a)
{
    double pval;
    
    // Force base F = A*tau
 
    count=0;
    if(p->W110==2 || p->W110==3)
    ULOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));    
        

        if(p->W111==1)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        
        if(p->W111==2)
        pval=0.5*(a->press(i,j,k)+a->press(i+1,j,k));
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            pval=0.5*(a->press(i,j,k)+a->press(i+1,j,k));
            
            if(phival>=p->W112*p->DXM)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        }
        
        if(p->W111==4)
        pval=0.25*(a->press(i,j,k)+a->press(i+1,j,k)) + 0.5*phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        
        
        if(p->W111==5)
        {
            if(p->count>=10)
            pval=0.5*(a->press(i,j,k)+a->press(i+1,j,k));
            
            if(p->count<10)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22);
        }
        
        // Yield Stress
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        tau0=tanphi*pval + p->W102_c;
        
        if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
        tau0 = (tanphi*pval + p->W102_c);
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // rho_water = 1000.0, new input?
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // m_p is new input W 104 
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c);    // m_u also use new input W 104

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==3)
        tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->u(i,j,k))>1.0e-20?(a->u(i,j,k)/fabs(a->u(i,j,k))):0.0;
        
        a->rhsvec.V[count] -= H*(tau0/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))))*f;

        ++count;
    
	}
    
    
    
    // Gradient based
    
    double pval1,pval2;
    double tau01,tau02;
    double fxx,fxy,fxz;
    double pvalx1,pvalx2,pvaly1,pvaly2,pvalz1,pvalz2;
    double pvaldxx1,pvaldxx2,pvaldxy1,pvaldxy2,pvaldxz1,pvaldxz2;
    double pvaldyx1,pvaldyx2,pvaldyy1,pvaldyy2,pvaldyz1,pvaldyz2;
    
    
    count=0;
    if(p->W110==4 || p->W110==5)
    ULOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 
       
        

        if(p->W111==1)
        {
        pval1 = phival*a->ro(i,j,k)*fabs(p->W22);
        pval2 = phival*a->ro(i+1,j,k)*fabs(p->W22);
        }
        
        if(p->W111==2)
        {
        pval1 = a->press(i,j,k);
        pval2 = a->press(i+1,j,k);
        }
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            {
            pval1 = a->press(i,j,k);
            pval2 = a->press(i+1,j,k);
            }
            
            if(phival>=p->W112*p->DXM)
            {
            pval1 = phival*a->ro(i,j,k)*fabs(p->W22);
            pval2 = phival*a->ro(i+1,j,k)*fabs(p->W22);
            }
        }


        
        // Yield Stress
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        { 
        tau01 = tanphi*pval1 + p->W102_c;
        tau02 = tanphi*pval2 + p->W102_c;
        }
        
        if(p->W101==2)  // HB-C dry sand, without MAX 
        {
        tau01 = (tanphi*pval1 + p->W102_c);
        tau02 = (tanphi*pval2 + p->W102_c);
        }
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        {
        tau01 = MAX(0.0,tanphi*pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    
        tau02 = MAX(0.0,tanphi*pval2*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k) + p->W102_c);   
        }
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        {
        tau01 = MAX(0.0,tanphi*pval1*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);  
        tau02 = MAX(0.0,tanphi*pval2*exp(-p->W104*gamma)*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k) + p->W102_c);   
        }
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        {
        tau01 = MAX(0.0,tanphi*MAX(0.0,pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c); 
        tau02 = MAX(0.0,tanphi*MAX(0.0,pval2*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k)-p->W104*gamma) + p->W102_c); 
        }   

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==5)
        {
        tau01 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        tau02 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i+1,j,k);
        }
        
        f = fabs(a->u(i,j,k))>1.0e-20?(a->u(i,j,k)/fabs(a->u(i,j,k))):0.0;
        
        a->rhsvec.V[count] += H*((tau02-tau01)/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)))); // *f ?

        ++count;
	}
    
    
    count=0;
    if(p->W110==6 || p->W110==7)
    ULOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i+1,j,k));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 
        // PFo: Adding direction to the yield stress gradient below. Should be *f for each term, where f = (pu_idx_j)/fabs(pu_idx_j) 
        // f(pwdx)*tau_x-term, f(pwdy)*tau_y-term, f(pwdz)*tau_z-term - PFo not sure how to write this
        // Try with u_x first instead: f = (a->u(i,j,k)/fabs(a->u(i,j,k)))
        
        fxx = fabs(pudx(p,a))>1.0e-20?(pudx(p,a)/fabs(pudx(p,a))):0.0; // pudx
        fxy = fabs(pudy(p,a))>1.0e-20?(pudy(p,a)/sqrt(pudy(p,a)*pudy(p,a)+pwdy(p,a)*pwdy(p,a))):0.0; // pudy
        fxz = fabs(pudz(p,a))>1.0e-20?(pudz(p,a)/sqrt(pudz(p,a)*pudz(p,a)+pvdz(p,a)*pvdz(p,a))):0.0; // pudz

        // (Hydrostatic) pressure at cell surfaces, based on level set function:

        //pvalx1 = 0.5*(a->phi(i,j,k)+a->phi(i-1,j,k))*0.5*(a->ro(i,j,k)+a->ro(i-1,j,k))*fabs(p->W22);
        //pvalx2 = 0.5*(a->phi(i+1,j,k)+a->phi(i,j,k))*0.5*(a->ro(i+1,j,k)+a->ro(i,j,k))*fabs(p->W22);
//        pvaly1 = 0.5*(a->phi(i,j,k)+a->phi(i,j-1,k))*0.5*(a->ro(i,j,k)+a->ro(i,j-1,k))*fabs(p->W22);
//        pvaly2 = 0.5*(a->phi(i,j+1,k)+a->phi(i,j,k))*0.5*(a->ro(i,j+1,k)+a->ro(i,j,k))*fabs(p->W22);
//        pvalz1 = 0.5*(a->phi(i,j,k)+a->phi(i,j,k-1))*0.5*(a->ro(i,j,k)+a->ro(i,j,k-1))*fabs(p->W22);
//        pvalz2 = 0.5*(a->phi(i,j,k+1)+a->phi(i,j,k))*0.5*(a->ro(i,j,k+1)+a->ro(i,j,k))*fabs(p->W22);
        
        // (Hydrostatic) pressure X-gradient (dp/dx) at cell surfaces, based on level set function:

        //pvaldxx1 = fabs(p->W22)*(a->phi(i,j,k)*a->ro(i,j,k) - a->phi(i-1,j,k)*a->ro(i-1,j,k))/(p->DXM);
        //pvaldxx2 = fabs(p->W22)*(a->phi(i+1,j,k)*a->ro(i+1,j,k) - a->phi(i,j,k)*a->ro(i,j,k))/(p->DXM);
//        pvaldxy1 = fabs(p->W22)*(0.5*(a->phi(i+1,j-1,k)+a->phi(i+1,j,k))*0.5*(a->ro(i+1,j-1,k)+a->ro(i+1,j,k)) 
//        - 0.5*(a->phi(i-1,j-1,k)+a->phi(i-1,j,k))*0.5*(a->ro(i-1,j-1,k)+a->ro(i-1,j,k)))/(2*p->DXM);
//        pvaldxy2 = fabs(p->W22)*(0.5*(a->phi(i+1,j,k)+a->phi(i+1,j+1,k))*0.5*(a->ro(i+1,j,k)+a->ro(i+1,j+1,k)) 
//        - 0.5*(a->phi(i-1,j,k)+a->phi(i-1,j+1,k))*0.5*(a->ro(i-1,j,k)+a->ro(i-1,j+1,k)))/(2*p->DXM);
//        pvaldxz1 = fabs(p->W22)*(0.5*(a->phi(i+1,j,k-1)+a->phi(i+1,j,k))*0.5*(a->ro(i+1,j,k-1)+a->ro(i+1,j,k)) 
//        - 0.5*(a->phi(i-1,j,k-1)+a->phi(i-1,j,k))*0.5*(a->ro(i-1,j,k-1)+a->ro(i-1,j,k)))/(2*p->DXM);
//        pvaldxz2 = fabs(p->W22)*(0.5*(a->phi(i+1,j,k)+a->phi(i+1,j,k+1))*0.5*(a->ro(i+1,j,k)+a->ro(i+1,j,k+1)) 
//        - 0.5*(a->phi(i-1,j,k)+a->phi(i-1,j,k+1))*0.5*(a->ro(i-1,j,k)+a->ro(i-1,j,k+1)))/(2*p->DXM);
        
        // (Hydrostatic) pressure Y-gradient (dp/dy) at cell surfaces, based on level set function:

        //pvaldyx1 = fabs(p->W22)*(0.5*(a->phi(i-1,j+1,k)+a->phi(i,j+1,k))*0.5*(a->ro(i-1,j+1,k)+a->ro(i,j+1,k)) 
        //- 0.5*(a->phi(i-1,j-1,k)+a->phi(i,j-1,k))*0.5*(a->ro(i-1,j-1,k)+a->ro(i,j-1,k)))/(2*p->DXM);
        //pvaldyx2 = fabs(p->W22)*(0.5*(a->phi(i,j+1,k)+a->phi(i+1,j+1,k))*0.5*(a->ro(i,j+1,k)+a->ro(i+1,j+1,k)) 
        //- 0.5*(a->phi(i,j-1,k)+a->phi(i+1,j-1,k))*0.5*(a->ro(i,j-1,k)+a->ro(i+1,j-1,k)))/(2*p->DXM);
//        pvaldyy1 = fabs(p->W22)*(a->phi(i,j,k)*a->ro(i,j,k) - a->phi(i,j-1,k)*a->ro(i,j-1,k))/(p->DXM);
//        pvaldyy2 = fabs(p->W22)*(a->phi(i,j+1,k)*a->ro(i,j+1,k) - a->phi(i,j,k)*a->ro(i,j,k))/(p->DXM);
        //pvaldyz1 = fabs(p->W22)*(0.5*(a->phi(i,j+1,k-1)+a->phi(i,j+1,k))*0.5*(a->ro(i,j+1,k-1)+a->ro(i,j+1,k)) 
        //- 0.5*(a->phi(i,j-1,k-1)+a->phi(i,j-1,k))*0.5*(a->ro(i,j-1,k-1)+a->ro(i,j-1,k)))/(2*p->DXM);
        //pvaldyz2 = fabs(p->W22)*(0.5*(a->phi(i,j+1,k)+a->phi(i,j+1,k+1))*0.5*(a->ro(i,j+1,k)+a->ro(i,j+1,k+1)) 
        //- 0.5*(a->phi(i,j-1,k)+a->phi(i,j-1,k+1))*0.5*(a->ro(i,j-1,k)+a->ro(i,j-1,k+1)))/(2*p->DXM);

        // Pressure at cell surfaces, based on dynamic pressure:

        pvalx1 = 0.5*(a->press(i,j,k)+a->press(i-1,j,k));
        pvalx2 = 0.5*(a->press(i+1,j,k)+a->press(i,j,k));
        pvaly1 = 0.5*(a->press(i,j,k)+a->press(i,j-1,k));
        pvaly2 = 0.5*(a->press(i,j+1,k)+a->press(i,j,k));
        pvalz1 = 0.5*(a->press(i,j,k)+a->press(i,j,k-1));
        pvalz2 = 0.5*(a->press(i,j,k+1)+a->press(i,j,k));
        
        // Pressure X-gradient (dp/dx) at cell surfaces, based on dynamic pressure:

        //pvaldxx1 = ;
        //pvaldxx2 = ;
        //pvaldxy1 = (0.5*(a->press(i+1,j-1,k)+a->press(i+1,j,k))  
        //- 0.5*(a->press(i-1,j-1,k)+a->press(i-1,j,k)))/(2*p->DXM);
        //pvaldxy2 = (0.5*(a->press(i+1,j,k)+a->press(i+1,j+1,k)) 
        //- 0.5*(a->press(i-1,j,k)+a->press(i-1,j+1,k)))/(2*p->DXM);
        //pvaldxz1 = (0.5*(a->press(i+1,j,k-1)+a->press(i+1,j,k)) 
        //- 0.5*(a->press(i-1,j,k-1)+a->press(i-1,j,k)))/(2*p->DXM);
        //pvaldxz2 = (0.5*(a->press(i+1,j,k)+a->press(i+1,j,k+1)) 
        //- 0.5*(a->press(i-1,j,k)+a->press(i-1,j,k+1)))/(2*p->DXM);
        
        // Pressure Y-gradient (dp/dy) at cell surfaces, based on dynamic pressure:

        //pvaldyx1 = ;
        //pvaldyx2 = ;
        //pvaldyy1 = (a->press(i,j,k) - a->press(i,j-1,k))/(p->DXM);
        //pvaldyy2 = (a->press(i,j+1,k) - a->press(i,j,k))/(p->DXM);
        //pvaldyz1 = ;
        //pvaldyz2 = ;

        
        //x-direction: dtau_zx (k+1 - k) + dtau_yx (j+1 - j)
        //y-direction: dtau_zy (k+1 - k) + dtau_xy (i+1 - i)
        //z-direction: dtau_xz (i+1 - i) + dtau_yz (j+1 - j)
        
        // Gradient based, uten tau_fill, bruk p-values directly uten å gå via tau_fill. Comparable? correct??
        
//        a->rhsvec.V[count] += H*( 
//                                 (-(pvaldyy1+pvaldxy2)*pvaly2/(0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22)) + (pvaldyy1+pvaldxy1)*pvaly1/(0.5*(a->ro(i,j,k)+a->ro(i,j-1,k))*fabs(p->W22))) //dtau_yx (j+1 - j) - temp: ro(i,j,k) instead of average, sum dp/dx and dp/dy instead of better term
//                                 +(-(pvaldxz2*pvalz2/(0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22))) + (pvaldxz1*pvalz1/(0.5*(a->ro(i,j,k)+a->ro(i,j,k-1))*fabs(p->W22))))  //dtau_zx (k+1 - k) - temp: ro(i,j,k) instead of average
//                                )/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)));
       
//        a->rhsvec.V[count] += H* (MIN(fxx*...,-pvaldyy1...)+MIN(fxx*...,-pvaldyy2...)
//                                fxx*tau_x(i+1,j,k)      - fxx*tau_x(i,j,k)
//                                
//                                -(pvaldyy1+pvaldxy2)*pvaly2/(0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22))  + (pvaldyy1+pvaldxy1)*pvaly1/(0.5*(a->ro(i,j,k)+a->ro(i,j-1,k))*fabs(p->W22))
//                                +fxy*0.5*(tau_y(i,j+1,k)+tau_y(i+1,j+1,k)) - fxy*0.5*(tau_y(i,j-1,k)+tau_y(i+1,j-1,k))
//                                
//                                -(pvaldxz2*pvalz2/(0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22))) + (pvaldxz1*pvalz1/(0.5*(a->ro(i,j,k)+a->ro(i,j,k-1))*fabs(p->W22)))
//                                +fxz*0.5*(tau_z(i,j,k+1)+tau_z(i+1,j,k+1)) - fxz*0.5*(tau_z(i,j,k-1)+tau_z(i+1,j,k-1))
//                                +)
//                                /(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k)));
                         
        a->rhsvec.V[count] += H*tanphi*(fxx*(pvalx2 - pvalx1) 
                                + fxy*(pvaly2 - pvaly1)
                                + fxz*(pvalz2 - pvalz1)   )
                                /(p->DXM*(a->ro(i,j,k)));
        //a->rhsvec.V[count] += H*((tau_x(i+1,j,k)-tau_x(i,j,k) 
        //                        + 0.5*(tau_y(i,j+1,k)+tau_y(i+1,j+1,k)) - 0.5*(tau_y(i,j-1,k)+tau_y(i+1,j-1,k))
        //                        + 0.5*(tau_z(i,j,k+1)+tau_z(i+1,j,k+1)) - 0.5*(tau_z(i,j,k-1)+tau_z(i+1,j,k-1))   )
        //                        /(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))));

        ++count;
    }
}

void rheology_f::v_source(lexer *p, fdm *a)
{
    double pval;
    double fyx,fyy,fyz;
    double pvalx1,pvalx2,pvaly1,pvaly2,pvalz1,pvalz2;
    double pvaldxx1,pvaldxx2,pvaldxy1,pvaldxy2,pvaldxz1,pvaldxz2;
    double pvaldyx1,pvaldyx2,pvaldyy1,pvaldyy2,pvaldyz1,pvaldyz2;    
    
    count=0;
    if(p->W110>=2)
    VLOOP
	{
    phival = 0.5*(a->phi(i,j,k)+a->phi(i,j+1,k));
    
    if(phival>epsi)
    H=1.0;

    if(phival<-epsi)
    H=0.0;

    if(fabs(phival)<=epsi)
    H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));  
    

    if(p->W111==1)
    pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22);
    
    if(p->W111==2)
    pval=0.5*(a->press(i,j,k)+a->press(i,j+1,k));
    
    if(p->W111==3)
    {
        if(phival<p->W112*p->DXM)
        pval=0.5*(a->press(i,j,k)+a->press(i,j+1,k));
        
        if(phival>=p->W112*p->DXM)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22);
    }
    
    if(p->W111==4)
    pval=0.25*(a->press(i,j,k)+a->press(i,j+1,k)) + 0.5*phival*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22);
    
    if(p->W111==5)
    {
        if(p->count>=10)
        pval=0.5*(a->press(i,j,k)+a->press(i,j+1,k));
        
        if(p->count<10)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22);
    }
    
    // Yield Stress
    if(p->W101==0)
    tau0=p->W96;
    
    if(p->W101==1)  // HB-C dry sand
    tau0=tanphi*pval + p->W102_c;
    
    if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
    tau0 = (tanphi*pval + p->W102_c);
        
    if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
    tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // rho_water = 1000.0, new input?
        
    if(p->W101==4)  // HB-C shear rate generated excess pore pressure
    tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // m_p is new input W 104 
        
    if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
    tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c);    // m_u also use new input W 104

        if(p->count==0)
        tau0=p->W96;
    
    if(p->W110==5)
    tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
    
	f = fabs(a->v(i,j,k))>1.0e-20?(a->v(i,j,k)/fabs(a->v(i,j,k))):0.0;
    
    a->rhsvec.V[count] -= H*(tau0/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))))*f;
    
	++count;
	}
    
    
    count=0;
    if(p->W110==6 || p->W110==7)
    VLOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 
        
        // PFo: f = (a->v(i,j,k)/fabs(a->v(i,j,k)))        
        fyx = fabs(pvdx(p,a))>1.0e-20?(pvdx(p,a)/sqrt(pvdx(p,a)*pvdx(p,a)+pwdx(p,a)*pwdx(p,a))):0.0; // pvdx pudx(p,a)
        fyy = fabs(pvdy(p,a))>1.0e-20?(pvdy(p,a)/fabs(pvdy(p,a))):0.0; // pvdy
        fyz = fabs(pvdz(p,a))>1.0e-20?(pvdz(p,a)/sqrt(pvdz(p,a)*pvdz(p,a)+pudz(p,a)*pudz(p,a))):0.0; // pvdz

        // (Hydrostatic) pressure at cell surfaces, based on level set function:

        pvalx1 = 0.5*(a->press(i,j,k)+a->press(i-1,j,k));
        pvalx2 = 0.5*(a->press(i+1,j,k)+a->press(i,j,k));
        pvaly1 = 0.5*(a->press(i,j,k)+a->press(i,j-1,k));
        pvaly2 = 0.5*(a->press(i,j+1,k)+a->press(i,j,k));
        pvalz1 = 0.5*(a->press(i,j,k)+a->press(i,j,k-1));
        pvalz2 = 0.5*(a->press(i,j,k+1)+a->press(i,j,k));
        
        // (Hydrostatic) pressure X-gradient (dp/dx) at cell surfaces, based on level set function:

        //pvaldxx1 = (a->press(i,j,k) - a->press(i-1,j,k))/(p->DXM);
        //pvaldxx2 = (a->press(i+1,j,k) - a->press(i,j,k))/(p->DXM);
        //pvaldxy1 = ;
        //pvaldxy2 = ;
        //pvaldxz1 = ;
        //pvaldxz2 = ;
        
        // (Hydrostatic) pressure Y-gradient (dp/dy) at cell surfaces, based on level set function:

        //pvaldyx1 = (0.5*(a->press(i-1,j+1,k)+a->press(i,j+1,k)) 
        //- 0.5*(a->press(i-1,j-1,k)+a->press(i,j-1,k)))/(2*p->DXM);
        //pvaldyx2 = (0.5*(a->press(i,j+1,k)+a->press(i+1,j+1,k)) 
        //- 0.5*(a->press(i,j-1,k)+a->press(i+1,j-1,k)))/(2*p->DXM);
        //pvaldyy1 = ;
        //pvaldyy2 = ;
        //pvaldyz1 = (0.5*(a->press(i,j+1,k-1)+a->press(i,j+1,k)) 
        //- 0.5*(a->press(i,j-1,k-1)+a->press(i,j-1,k)))/(2*p->DXM);
        //pvaldyz2 = (0.5*(a->press(i,j+1,k)+a->press(i,j+1,k+1)) 
        //- 0.5*(a->press(i,j-1,k)+a->press(i,j-1,k+1)))/(2*p->DXM);
        
        //x-direction: dtau_zx (k+1 - k) + dtau_yx (j+1 - j)
        //y-direction: dtau_zy (k+1 - k) + dtau_xy (i+1 - i)
        //z-direction: dtau_xz (i+1 - i) + dtau_yz (j+1 - j)

//        a->rhsvec.V[count] += H*( 
//                                 (-(pvaldyx2+pvaldxx2)*pvalx2/(0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22)) + (pvaldyx1+pvaldxx1)*pvalx1/(0.5*(a->ro(i,j,k)+a->ro(i-1,j,k))*fabs(p->W22))) //dtau_xy (i+1 - i) - temp: ro(i,j,k) instead of average, sum dp/dx and dp/dy instead of better term
//                                 +(-(pvaldyz2*pvalz2/(0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22))) + (pvaldyz1*pvalz1/(0.5*(a->ro(i,j,k)+a->ro(i,j,k-1))*fabs(p->W22))))  //dtau_zy (k+1 - k) - temp: ro(i,j,k) instead of average
//                                )/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k)));

        a->rhsvec.V[count] += H*tanphi*(fyx*(pvalx2 - pvalx1) 
                                + fyy*(pvaly2 - pvaly1)
                                + fyz*(pvalz2 - pvalz1)   )
                                /(p->DXM*(a->ro(i,j,k)));
      
        //a->rhsvec.V[count] += H*((fyx*(0.5*(tau_x(i+1,j,k)+tau_x(i+1,j+1,k)) - 0.5*(tau_x(i-1,j,k)+tau_x(i-1,j+1,k)))
        //                        + fyy*(tau_y(i,j+1,k)-tau_y(i,j,k)) 
        //                        + fyz*(0.5*(tau_z(i,j,k+1)+tau_z(i,j+1,k+1)) - 0.5*(tau_z(i,j,k-1)+tau_z(i,j+1,k-1)))    )
        //                        /(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))));
        //a->rhsvec.V[count] += H*((0.5*(tau_x(i+1,j,k)+tau_x(i+1,j+1,k)) - 0.5*(tau_x(i-1,j,k)+tau_x(i-1,j+1,k))
        //                        + tau_y(i,j+1,k)-tau_y(i,j,k) 
        //                        + 0.5*(tau_z(i,j,k+1)+tau_z(i,j+1,k+1)) - 0.5*(tau_z(i,j,k-1)+tau_z(i,j+1,k-1))    )
        //                        /(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))));

        ++count;
    }
}

void rheology_f::w_source(lexer *p, fdm *a)
{
    double pval;
    
    count=0;
    if(p->W110==2 || p->W110==3)
    WLOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));  
        

        if(p->W111==1)
        pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22);
        
        if(p->W111==2)
        pval=0.5*(a->press(i,j,k)+a->press(i,j,k+1));
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            pval=0.5*(a->press(i,j,k)+a->press(i,j,k+1));
            
            if(phival>=p->W112*p->DXM)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22);
        }
        
        if(p->W111==4)
        pval=0.5*0.5*(a->press(i,j,k)+a->press(i,j,k+1)) + 0.5*phival*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22);
        
        if(p->W111==5)
        {
            if(p->count>=10)
            pval=0.5*(a->press(i,j,k)+a->press(i,j,k+1));
            
            if(p->count<10)
            pval=phival*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*fabs(p->W22);
        }
        
        
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        tau0=tanphi*pval + p->W102_c;
        
        if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
        tau0 = (tanphi*pval + p->W102_c);
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // rho_water = 1000.0, new input?
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // m_p is new input W 104 
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c);    // m_u also use new input W 104
        
            if(p->count==0)
            tau0=p->W96;
            
        if(p->W110==3)
        tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->w(i,j,k))>1.0e-20?(a->w(i,j,k)/fabs(a->w(i,j,k))):0.0;
        
        a->rhsvec.V[count] -= H*(tau0/(0.5*(p->DXM*a->ro(i,j,k)+a->ro(i,j,k+1))))*f;

        ++count;
	}
    
    
    
    // Gradient Based
    
    double pval1,pval2;
    double tau01,tau02;
    double fzx,fzy,fzz;
    double pvalx1,pvalx2,pvaly1,pvaly2,pvalz1,pvalz2;
    double pvaldxx1,pvaldxx2,pvaldxy1,pvaldxy2,pvaldxz1,pvaldxz2;
    double pvaldyx1,pvaldyx2,pvaldyy1,pvaldyy2,pvaldyz1,pvaldyz2;
    
    count=0;
    if(p->W110==4 || p->W110==5)
    WLOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 
       
        

        if(p->W111==1)
        {
        pval1 = phival*a->ro(i,j,k)*fabs(p->W22);
        pval2 = phival*a->ro(i,j,k+1)*fabs(p->W22);
        }
        
        if(p->W111==2)
        {
        pval1 = a->press(i,j,k);
        pval2 = a->press(i,j,k+1);
        }
        
        if(p->W111==3)
        {
            if(phival<p->W112*p->DXM)
            {
            pval1 = a->press(i,j,k);
            pval2 = a->press(i,j,k+1);
            }
            
            if(phival>=p->W112*p->DXM)
            {
            pval1 = phival*a->ro(i,j,k)*fabs(p->W22);
            pval2 = phival*a->ro(i,j,k+1)*fabs(p->W22);
            }
        }


        
        // Yield Stress
        if(p->W101==0)
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        { 
        tau01 = tanphi*pval1 + p->W102_c;
        tau02 = tanphi*pval2 + p->W102_c;
        }
        
        if(p->W101==2)  // HB-C dry sand, without MAX 
        {
        tau01 = (tanphi*pval1 + p->W102_c);
        tau02 = (tanphi*pval2 + p->W102_c);
        }
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        {
        tau01 = MAX(0.0,tanphi*pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    
        tau02 = MAX(0.0,tanphi*pval2*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1) + p->W102_c);   
        }
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        {
        tau01 = MAX(0.0,tanphi*pval1*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);  
        tau02 = MAX(0.0,tanphi*pval2*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1) + p->W102_c);   
        }
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        {
        tau01 = MAX(0.0,tanphi*MAX(0.0,pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c); 
        tau02 = MAX(0.0,tanphi*MAX(0.0,pval2*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1)-p->W104*gamma) + p->W102_c); 
        }   

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==5)
        {
        tau01 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        tau02 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k+1);
        }
        
        f = fabs(a->w(i,j,k))>1.0e-20?(a->w(i,j,k)/fabs(a->w(i,j,k))):0.0;
        
        a->rhsvec.V[count] += H*((tau02-tau01)/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1)))); // PFo: Missing "*f" ?

        ++count;
	}
    
    
    count=0;
    if(p->W110==6 || p->W110==7)
    WLOOP
	{
        phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k+1));
        
        if(phival>epsi)
        H=1.0;

        if(phival<-epsi)
        H=0.0;

        if(fabs(phival)<=epsi)
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi)); 
        
        // f = (a->w(i,j,k)/fabs(a->w(i,j,k)))
        fzx = fabs(pwdx(p,a))>1.0e-20?(pwdx(p,a)/sqrt(pwdx(p,a)*pwdx(p,a)+pvdx(p,a)*pvdx(p,a))):0.0; // pwdx
        fzy = fabs(pwdy(p,a))>1.0e-20?(pwdy(p,a)/sqrt(pwdy(p,a)*pwdy(p,a)+pudy(p,a)*pudy(p,a))):0.0; // pwdy
        fzz = fabs(pwdz(p,a))>1.0e-20?(pwdz(p,a)/fabs(pwdz(p,a))):0.0; // pwdz
        
        // (Hydrostatic) pressure at cell surfaces, based on level set function:

        pvalx1 = 0.5*(a->press(i,j,k)+a->press(i-1,j,k));
        pvalx2 = 0.5*(a->press(i+1,j,k)+a->press(i,j,k));
        pvaly1 = 0.5*(a->press(i,j,k)+a->press(i,j-1,k));
        pvaly2 = 0.5*(a->press(i,j+1,k)+a->press(i,j,k));
        pvalz1 = 0.5*(a->press(i,j,k)+a->press(i,j,k-1));
        pvalz2 = 0.5*(a->press(i,j,k+1)+a->press(i,j,k));
        
        // (Hydrostatic) pressure X-gradient (dp/dx) at cell surfaces, based on level set function:

        //pvaldxx1 = (a->press(i,j,k) - a->press(i-1,j,k))/(p->DXM);
        //pvaldxx2 = (a->press(i+1,j,k) - a->press(i,j,k))/(p->DXM);
        //pvaldxy1 = ;
        //pvaldxy2 = ;
        //pvaldxz1 = ;
        //pvaldxz2 = ;
        
        // (Hydrostatic) pressure Y-gradient (dp/dy) at cell surfaces, based on level set function:

        //pvaldyx1 = ;
        //pvaldyx2 = ;
        //pvaldyy1 = (a->press(i,j,k) - a->press(i,j-1,k))/(p->DXM);
        //pvaldyy2 = (a->press(i,j+1,k) - a->press(i,j,k))/(p->DXM);
        //pvaldyz1 = ;
        //pvaldyz2 = ;
        
        //x-direction: dtau_zx (k+1 - k) + dtau_yx (j+1 - j)
        //y-direction: dtau_zy (k+1 - k) + dtau_xy (i+1 - i)
        //z-direction: dtau_xz (i+1 - i) + dtau_yz (j+1 - j)

//        a->rhsvec.V[count] += H*( 
//                                 (-(pvaldxx2*pvalx2/(0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*fabs(p->W22))) + (pvaldxx1*pvalx1/(0.5*(a->ro(i,j,k)+a->ro(i-1,j,k))*fabs(p->W22)))) //dtau_xz (i+1 - i) - temp: ro(i,j,k) instead of average, sum dp/dx and dp/dy instead of better term
//                                 +(-(pvaldyy2*pvaly2/(0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*fabs(p->W22))) + (pvaldyy1*pvaly1/(0.5*(a->ro(i,j,k)+a->ro(i,j-1,k))*fabs(p->W22))))  //dtau_yz (j+1 - j) - temp: ro(i,j,k) instead of average
//                                )/(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1)));

        a->rhsvec.V[count] += H*tanphi*(fzx*(pvalx2 - pvalx1) 
                                + fzy*(pvaly2 - pvaly1)
                                + fzz*(pvalz2 - pvalz1)   )
                                /(p->DXM*(a->ro(i,j,k)));
        
        //a->rhsvec.V[count] += H*((fzx*(0.5*(tau_x(i+1,j,k)+tau_x(i+1,j,k+1)) - 0.5*(tau_x(i-1,j,k)+tau_x(i-1,j,k+1)))
        //                        + fzy*(0.5*(tau_y(i,j+1,k)+tau_y(i,j+1,k+1)) - 0.5*(tau_y(i,j-1,k)+tau_y(i,j-1,k+1)))
        //                        + fzz*(tau_z(i,j,k+1) - tau_z(i,j,k))    ) 
        //                        /(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))));
        //a->rhsvec.V[count] += H*((0.5*(tau_x(i+1,j,k)+tau_x(i+1,j,k+1)) - 0.5*(tau_x(i-1,j,k)+tau_x(i-1,j,k+1))
        //                        + 0.5*(tau_y(i,j+1,k)+tau_y(i,j+1,k+1)) - 0.5*(tau_y(i,j-1,k)+tau_y(i,j-1,k+1))
        //                        + tau_z(i+1,j,k) - tau_z(i,j,k)    )
        //                        /(p->DXM*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))));

        ++count;
    }
    
}*/
