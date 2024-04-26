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
Authors: Tobias Martin, Hans Bihs, Ahmet Soydan
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::solid_forcing(lexer *p, fdm *a, double alpha, field& uvel, field &vvel, field& wvel,
                              field& fx, field &fy, field &fz)
{
    // ghostcell update
    gcdf_update(p,a);
    gcb_velflagio(p,a);
    
     // Reset heaviside field
    ULOOP
    a->fbh1(i,j,k) = 0.0;

    VLOOP
    a->fbh2(i,j,k) = 0.0;
    
    WLOOP
    a->fbh3(i,j,k) = 0.0;

    LOOP
    a->fbh4(i,j,k) = 0.0;

    start1(p,a->fbh1,10);
    start2(p,a->fbh2,11);
    start3(p,a->fbh3,12);
    start4(p,a->fbh4,40);

// Calculate forcing fields
    double H,Ht, uf, vf, wf;
	double nx, ny, nz,norm ;
	double psi, phival_sf;
    double dirac;
    
    if(p->B20==2)
    {
        
    ULOOP
    {
        uf = 0.0;
        H = Hsolidface(p,a,1,0,0);
       
        fx(i,j,k) += H*(uf - uvel(i,j,k))/(alpha*p->dt);   
        a->fbh1(i,j,k) = min(a->fbh1(i,j,k) + H, 1.0); 
    }
    
    VLOOP
    {
        vf = 0.0;
        H = Hsolidface(p,a,0,1,0);
       
        fy(i,j,k) += H*(vf - vvel(i,j,k))/(alpha*p->dt);
        a->fbh2(i,j,k) = min(a->fbh2(i,j,k) + H, 1.0); 
    }
    
    WLOOP
    {
        wf = 0.0;
        H = Hsolidface(p,a,0,0,1);

        fz(i,j,k) += H*(wf - wvel(i,j,k))/(alpha*p->dt);
        a->fbh3(i,j,k) = min(a->fbh3(i,j,k) + H, 1.0); 
    }
    
    LOOP
    {
        H = Hsolidface(p,a,0,0,0);
        a->fbh4(i,j,k) = min(a->fbh4(i,j,k) + H, 1.0); 
        //a->test(i,j,k) = a->fbh4(i,j,k) ;
    }
    	
    psi = 1.1*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if (p->j_dir==0)
    psi = 1.1*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 

    LOOP
    {
        dirac = 0.0;
        if(fabs(MIN(a->solid(i,j,k),a->topo(i,j,k)))<psi)
        dirac = (0.5/psi)*(1.0 + cos((PI*(MIN(a->solid(i,j,k),a->topo(i,j,k))))/psi));
        
        a->fbh5(i,j,k) = 1.0-MIN(dirac,1.0);
    }
    
    }

// Construct solid heaviside function	
    if(p->B20==1)
    {
        
    ULOOP
    {
        uf = 0.0;
        
		// Normal vectors calculation 
		nx = -(a->topo(i+1,j,k) - a->topo(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->topo(i,j+1,k) - a->topo(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->topo(i,j,k+1) - a->topo(i,j,k-1))/(2.0*p->DZN[KP]);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;

		H = Hsolidface(p,a,1,0,0);
	    Ht = Hsolidface_t(p,a,1,0,0);
	
		// Level set function
		phival_sf = MIN(0.5*(a->solid(i,j,k) + a->solid(i+1,j,k)), 0.5*(a->topo(i,j,k) + a->topo(i+1,j,k))); 
        

		// Construct the field around the solid body to adjust the tangential velocity and calculate forcing
		if (phival_sf < 0.0)
		{
			fx(i,j,k) += H*(uf - uvel(i,j,k))/(alpha*p->dt); 
		}
        
		else if (phival_sf >0 && phival_sf<psi )
		{
			fx(i,j,k) +=   fabs(nx)*H*(uf - uvel(i,j,k))/(alpha*p->dt);
		}
        
		else
		{
			fx(i,j,k) += 0.0;
		}
	
        a->fbh1(i,j,k) = min(a->fbh1(i,j,k) + H, 1.0); 
    }
    
    VLOOP
    {
        vf = 0.0;
    
		// Normal vectors calculation 
		nx = -(a->topo(i+1,j,k) - a->topo(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->topo(i,j+1,k) - a->topo(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->topo(i,j,k+1) - a->topo(i,j,k-1))/(2.0*p->DZN[KP]);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;

        
         H = Hsolidface(p,a,0,1,0);
		Ht = Hsolidface_t(p,a,0,1,0);
		
      
		//Level set function
		phival_sf = MIN(0.5*(a->solid(i,j,k) + a->solid(i,j+1,k)), 0.5*(a->topo(i,j,k) + a->topo(i,j+1,k)));
	  
		//Construct the field around the solid body to adjust the tangential velocity and calculate forcing
	    if (phival_sf < 0.0)
		{
			fy(i,j,k) += H*(vf - vvel(i,j,k))/(alpha*p->dt); 
		}
		else if (phival_sf >0 && phival_sf<psi )
		{
            fy(i,j,k) +=   fabs(ny)*H*(vf - vvel(i,j,k))/(alpha*p->dt);
            
		}
		else
		{
			fy(i,j,k) += 0.0;
		}
	  
        a->fbh2(i,j,k) = min(a->fbh2(i,j,k) + H , 1.0); 
    }
	
    WLOOP
    {
        wf = 0.0;
        
		// Normal vectors calculation 
		nx = -(a->topo(i+1,j,k) - a->topo(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->topo(i,j+1,k) - a->topo(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->topo(i,j,k+1) - a->topo(i,j,k-1))/(2.0*p->DZN[KP]);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;

        
         H = Hsolidface(p,a,0,0,1);
		Ht = Hsolidface_t(p,a,0,0,1);


		// Level set function
		phival_sf = MIN(0.5*(a->solid(i,j,k) + a->solid(i,j,k+1)), 0.5*(a->topo(i,j,k) + a->topo(i,j,k+1)));
		
		// Construct the field around the solid body to adjust the tangential velocity and calculate forcing

		if (phival_sf < 0.0)
		{
			fz(i,j,k) += H*(wf - wvel(i,j,k))/(alpha*p->dt); 
		}
		else if (phival_sf >0 && phival_sf<psi )
		{
            fz(i,j,k) +=   fabs(nz)*H*(wf - wvel(i,j,k))/(alpha*p->dt);
		}
		else
		{
			fz(i,j,k) += 0.0;
		}
	
        a->fbh3(i,j,k) = min(a->fbh3(i,j,k) + H , 1.0); 
    }
    
    LOOP
    {
        H = Hsolidface(p,a,0,0,0);
		Ht = Hsolidface_t(p,a,0,0,0);
        a->fbh4(i,j,k) = min(a->fbh4(i,j,k) + H, 1.0); 
    }
    
    //double psi;
	
    psi = 1.1*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if (p->j_dir==0)
    psi = 1.1*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 

    
    LOOP
    {
        dirac = 0.0;
        if(fabs(MIN(a->solid(i,j,k),a->topo(i,j,k)))<psi)
        dirac = (0.5/psi)*(1.0 + cos((PI*(MIN(a->solid(i,j,k),a->topo(i,j,k))))/psi));
        
        a->fbh5(i,j,k) =  1.0-MIN(dirac,1.0);
    }
	
	}

    start1(p,a->fbh1,10);
    start2(p,a->fbh2,11);
    start3(p,a->fbh3,12);
    start4(p,a->fbh4,40);

    start1(p,fx,10);
    start2(p,fy,11);
    start3(p,fz,12);         
}

