/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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
            if(phival<p->W112*p->dx)
            pval=0.5*(a->press(i,j,k)+a->press(i+1,j,k));
            
            if(phival>=p->W112*p->dx)
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
        tau0 = (tanphi*pval + p->W102_c)*(1.0-exp(-p->W103*gamma));
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // rho_water = 1000.0, new input?
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_p is new input W 104 
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_u also use new input W 104

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==3)
        tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->u(i,j,k))>1.0e-20?(a->u(i,j,k)/fabs(a->u(i,j,k))):0.0;
        
        a->rhsvec.V[count] -= H*(tau0/(p->dx*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))))*f;

        ++count;
    
	}
    
    
    
    // Gradient based
    
    double pval1,pval2;
    double tau01,tau02;
    
    
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
            if(phival<p->W112*p->dx)
            {
            pval1 = a->press(i,j,k);
            pval2 = a->press(i+1,j,k);
            }
            
            if(phival>=p->W112*p->dx)
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
        tau01 = (tanphi*pval1 + p->W102_c)*(1.0-exp(-p->W103*gamma));
        tau02 = (tanphi*pval2 + p->W102_c)*(1.0-exp(-p->W103*gamma));
        }
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        {
        tau01 = MAX(0.0,tanphi*pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    
        tau02 = MAX(0.0,tanphi*pval2*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));   
        }
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        {
        tau01 = MAX(0.0,tanphi*pval1*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));  
        tau02 = MAX(0.0,tanphi*pval2*exp(-p->W104*gamma)*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));   
        }
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        {
        tau01 = MAX(0.0,tanphi*MAX(0.0,pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma)); 
        tau02 = MAX(0.0,tanphi*MAX(0.0,pval2*MAX(0.0,a->ro(i+1,j,k)-1000.0)/a->ro(i+1,j,k)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma)); 
        }   

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==5)
        {
        tau01 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        tau02 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i+1,j,k);
        }
        
        f = fabs(a->u(i,j,k))>1.0e-20?(a->u(i,j,k)/fabs(a->u(i,j,k))):0.0;
        
        a->rhsvec.V[count] += H*((tau02-tau01)/(p->dx*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))));

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
        
        a->rhsvec.V[count] += H*((tau_x(i+1,j,k)-tau_x(i,j,k) 
                                + 0.5*(tau_y(i,j+1,k)+tau_y(i+1,j+1,k)) - 0.5*(tau_y(i,j-1,k)+tau_y(i+1,j-1,k))
                                + 0.5*(tau_z(i,j,k+1)+tau_z(i+1,j,k+1)) - 0.5*(tau_z(i,j,k-1)+tau_z(i+1,j,k-1))   )
                                /(p->dx*0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))));

        ++count;
    }
}

void rheology_f::v_source(lexer *p, fdm *a)
{
    double pval;
    
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
        if(phival<p->W112*p->dx)
        pval=0.5*(a->press(i,j,k)+a->press(i,j+1,k));
        
        if(phival>=p->W112*p->dx)
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
    tau0 = (tanphi*pval + p->W102_c)*(1.0-exp(-p->W103*gamma));
        
    if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
    tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // rho_water = 1000.0, new input?
        
    if(p->W101==4)  // HB-C shear rate generated excess pore pressure
    tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_p is new input W 104 
        
    if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
    tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_u also use new input W 104

        if(p->count==0)
        tau0=p->W96;
    
    if(p->W110==5)
    tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
    
	f = fabs(a->v(i,j,k))>1.0e-20?(a->v(i,j,k)/fabs(a->v(i,j,k))):0.0;
    
    a->rhsvec.V[count] -= H*(tau0/(p->dx*0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))))*f;
    
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
        
        a->rhsvec.V[count] += H*((0.5*(tau_x(i+1,j,k)+tau_x(i+1,j+1,k)) - 0.5*(tau_x(i-1,j,k)+tau_x(i-1,j+1,k))
                                + tau_y(i,j+1,k)-tau_y(i,j,k) 
                                + 0.5*(tau_z(i,j,k+1)+tau_z(i,j+1,k+1)) - 0.5*(tau_z(i,j,k-1)+tau_z(i,j+1,k-1))    )
                                /(p->dx*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))));

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
            if(phival<p->W112*p->dx)
            pval=0.5*(a->press(i,j,k)+a->press(i,j,k+1));
            
            if(phival>=p->W112*p->dx)
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
        tau0 = (tanphi*pval + p->W102_c)*(1.0-exp(-p->W103*gamma));
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // rho_water = 1000.0, new input?
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pval*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_p is new input W 104 
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pval*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_u also use new input W 104
        
            if(p->count==0)
            tau0=p->W96;
            
        if(p->W110==3)
        tau0 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        
        f = fabs(a->w(i,j,k))>1.0e-20?(a->w(i,j,k)/fabs(a->w(i,j,k))):0.0;
        
        a->rhsvec.V[count] -= H*(tau0/(0.5*(p->dx*a->ro(i,j,k)+a->ro(i,j,k+1))))*f;

        ++count;
	}
    
    
    
    // Gradient Based
    
    double pval1,pval2;
    double tau01,tau02;
    
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
            if(phival<p->W112*p->dx)
            {
            pval1 = a->press(i,j,k);
            pval2 = a->press(i,j,k+1);
            }
            
            if(phival>=p->W112*p->dx)
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
        tau01 = (tanphi*pval1 + p->W102_c)*(1.0-exp(-p->W103*gamma));
        tau02 = (tanphi*pval2 + p->W102_c)*(1.0-exp(-p->W103*gamma));
        }
            
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        {
        tau01 = MAX(0.0,tanphi*pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));    
        tau02 = MAX(0.0,tanphi*pval2*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1) + p->W102_c)*(1.0-exp(-p->W103*gamma));   
        }
            
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        {
        tau01 = MAX(0.0,tanphi*pval1*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c)*(1.0-exp(-p->W103*gamma));  
        tau02 = MAX(0.0,tanphi*pval2*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1) + p->W102_c)*(1.0-exp(-p->W103*gamma));   
        }
            
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        {
        tau01 = MAX(0.0,tanphi*MAX(0.0,pval1*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma)); 
        tau02 = MAX(0.0,tanphi*MAX(0.0,pval2*MAX(0.0,a->ro(i,j,k+1)-1000.0)/a->ro(i,j,k+1)-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma)); 
        }   

            if(p->count==0)
            tau0=p->W96;
         
        if(p->W110==5)
        {
        tau01 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
        tau02 +=  ((p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k+1);
        }
        
        f = fabs(a->w(i,j,k))>1.0e-20?(a->w(i,j,k)/fabs(a->w(i,j,k))):0.0;
        
        a->rhsvec.V[count] += H*((tau02-tau01)/(p->dx*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))));

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
        
        a->rhsvec.V[count] += H*((0.5*(tau_x(i+1,j,k)+tau_x(i+1,j,k+1)) - 0.5*(tau_x(i-1,j,k)+tau_x(i-1,j,k+1))
                                + 0.5*(tau_y(i,j+1,k)+tau_y(i,j+1,k+1)) - 0.5*(tau_y(i,j-1,k)+tau_y(i,j-1,k+1))
                                + tau_z(i+1,j,k) - tau_z(i,j,k)    )
                                /(p->dx*0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))));

        ++count;
    }
    
}