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
Authors: Tobias Martin, Ahmet Soydan, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::updateForcing(lexer *p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, field &fx, field &fy, field &fz,int iter)
{
    // Determine floating body velocities
    Eigen::Matrix<double, 6, 1> u_fb;
    
        // U
        if(p->X11_u==0)
        u_fb(0) = 0.0;
        
        if(p->X11_u==1)
        u_fb(0) = p_(0)/Mass_fb;
        
        if(p->X11_u==2)
        u_fb(0) = Uext*ramp_vel(p);
        
        // V
        if(p->X11_v==0 || p->j_dir==0)
        u_fb(1) = 0.0;
        
        if(p->X11_u==1 && p->j_dir==1)
        u_fb(1) = p_(1)/Mass_fb;
        
        if(p->X11_u==2 && p->j_dir==1)
        u_fb(1) = Vext*ramp_vel(p);
        
        // W
        if(p->X11_w==0)
        u_fb(2) = 0.0;
        
        if(p->X11_w==1)
        u_fb(2) = p_(2)/Mass_fb;
        
        if(p->X11_w==2)
        u_fb(2) = Wext*ramp_vel(p);
        
        // rotation
        if(p->j_dir==0)
        {
        u_fb(3) = 0.0;
        u_fb(4) = omega_I(1);
        u_fb(5) = 0.0;
        }
        
        if(p->j_dir==1)
        {
        u_fb(3) = omega_I(0);
        u_fb(4) = omega_I(1);
        u_fb(5) = omega_I(2);
        }
        

// Calculate forcing fields
    double H,Ht, uf, vf, wf;
	double nx, ny, nz,norm ;
	double psi, phival_fb;
    double dirac;
    
    if(p->X14==1)
    {
        
    ULOOP
    {
        uf = u_fb(0) + u_fb(4)*(p->pos1_z() - c_(2)) - u_fb(5)*(p->pos1_y() - c_(1));
        H = Hsolidface(p,a,1,0,0);
       
        fx(i,j,k) += H*(uf - uvel(i,j,k))/(alpha[iter]*p->dt);   
        a->fbh1(i,j,k) = min(a->fbh1(i,j,k) + H, 1.0); 
    }
    VLOOP
    {
        vf = u_fb(1) + u_fb(5)*(p->pos2_x() - c_(0)) - u_fb(3)*(p->pos2_z() - c_(2));
        H = Hsolidface(p,a,0,1,0);
       
        fy(i,j,k) += H*(vf - vvel(i,j,k))/(alpha[iter]*p->dt);
        a->fbh2(i,j,k) = min(a->fbh2(i,j,k) + H, 1.0); 
    }
    WLOOP
    {
        wf = u_fb(2) + u_fb(3)*(p->pos3_y() - c_(1)) - u_fb(4)*(p->pos3_x() - c_(0));
        H = Hsolidface(p,a,0,0,1);

        fz(i,j,k) += H*(wf - wvel(i,j,k))/(alpha[iter]*p->dt);
        a->fbh3(i,j,k) = min(a->fbh3(i,j,k) + H, 1.0); 
    }
    LOOP
    {
        H = Hsolidface(p,a,0,0,0);
        a->fbh4(i,j,k) = min(a->fbh4(i,j,k) + H, 1.0); 
    }
    	
    psi = 1.1*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if (p->j_dir==0)
    psi = 1.1*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 

    LOOP
    {
        dirac = 0.0;
        if(fabs(a->fb(i,j,k))<psi)
        dirac = (0.5/psi)*(1.0 + cos((PI*(a->fb(i,j,k)))/psi));
        
        a->fbh5(i,j,k) = 1.0-MIN(dirac,1.0);
    }
    
    }

// Construct solid heaviside function	
    if(p->X14>=2)
    {
        
    ULOOP
    {
        uf = u_fb(0) + u_fb(4)*(p->pos1_z() - c_(2)) - u_fb(5)*(p->pos1_y() - c_(1));
        
		// Normal vectors calculation 
		nx = -(a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;

		H = Hsolidface(p,a,1,0,0);
	    Ht = Hsolidface_t(p,a,1,0,0);
	
	   //cout<<"Htx: "<<Ht<<endl;
		
		// Level set function
		phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i+1,j,k));	
		
		// Construct the field around the solid body to adjust the tangential velocity and calculate forcing
		if (phival_fb < 0)
		{
			fx(i,j,k) += H*(uf - uvel(i,j,k))/(alpha[iter]*p->dt); 
		}
		else if (phival_fb >0 && phival_fb<psi )
		{
            if(p->X14==2)
			fx(i,j,k) +=   fabs(nx)*H*(uf - uvel(i,j,k))/(alpha[iter]*p->dt);
            
            if(p->X14==3)
			fx(i,j,k) +=   (fabs(nx)*H + ((1.0-fabs(nx))*Ht))*(uf - uvel(i,j,k))/(alpha[iter]*p->dt);
		}
		else
		{
			fx(i,j,k) += 0;
		}
	
        a->fbh1(i,j,k) = min(a->fbh1(i,j,k) + H, 1.0); 
    }
    VLOOP
    {
        vf = u_fb(1) + u_fb(5)*(p->pos2_x() - c_(0)) - u_fb(3)*(p->pos2_z() - c_(2));
        
    
		// Normal vectors calculation 
		nx = -(a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;

        
        H = Hsolidface(p,a,0,1,0);
		Ht = Hsolidface_t(p,a,0,1,0);
		
      
		//Level set function
		phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i,j+1,k));
	  
		//Construct the field around the solid body to adjust the tangential velocity and calculate forcing
	    if (phival_fb < 0)
		{
			fy(i,j,k) += H*(vf - vvel(i,j,k))/(alpha[iter]*p->dt); 
		}
		else if (phival_fb >0 && phival_fb<psi )
		{
            if(p->X14==2)
            fy(i,j,k) +=   fabs(ny)*H*(vf - vvel(i,j,k))/(alpha[iter]*p->dt);
            
            if(p->X14==3)
            fy(i,j,k) +=   (fabs(ny)*H + ((1.0-fabs(ny))*Ht))*(vf - vvel(i,j,k))/(alpha[iter]*p->dt);
		}
		else
		{
			fy(i,j,k) += 0;
		}
	  
        a->fbh2(i,j,k) = min(a->fbh2(i,j,k) + H , 1.0); 
    }
	
    WLOOP
    {
        wf = u_fb(2) + u_fb(3)*(p->pos3_y() - c_(1)) - u_fb(4)*(p->pos3_x() - c_(0));
        

		// Normal vectors calculation 
		nx = -(a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;

        
         H = Hsolidface(p,a,0,0,1);
		Ht = Hsolidface_t(p,a,0,0,1);


		// Level set function
		phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i,j,k+1));
		
		// Construct the field around the solid body to adjust the tangential velocity and calculate forcing

		if (phival_fb < 0)
		{
			fz(i,j,k) += H*(wf - wvel(i,j,k))/(alpha[iter]*p->dt); 
		}
		else if (phival_fb >0 && phival_fb<psi )
		{
            if(p->X14==2)
            fz(i,j,k) +=   fabs(nz)*H*(wf - wvel(i,j,k))/(alpha[iter]*p->dt);
            
            if(p->X14==3)
            fz(i,j,k) +=   (fabs(nz)*H + ((1.0-fabs(nz))*Ht))*(wf - wvel(i,j,k))/(alpha[iter]*p->dt);
		}
		else
		{
			fz(i,j,k) += 0;
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
        if(fabs(a->fb(i,j,k))<psi)
        dirac = (0.5/psi)*(1.0 + cos((PI*(a->fb(i,j,k)))/psi));
        
        a->fbh5(i,j,k) =   1.0-MIN(dirac,1.0);
    }
	
	}

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);

    pgc->start1(p,fx,10);
    pgc->start2(p,fy,11);
    pgc->start3(p,fz,12);         
};

double sixdof_obj::Hsolidface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi, H, phival_fb,dirac;
	
    psi = p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if (p->knoy == 1)
    {
        psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
    }

    // Construct solid heaviside function
    phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i+aa,j+bb,k+cc));
	
    if (-phival_fb > psi)
    {
        H = 1.0;
    }
    else if (-phival_fb < -psi)
    {
        H = 0.0;
    }
    else
    {
        H = 0.5*(1.0 + -phival_fb/psi + (1.0/PI)*sin((PI*-phival_fb)/psi));
    }
	
    return H;
}

double sixdof_obj::Hsolidface_t(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi, H, phival_fb,dirac;
	

    if(p->j_dir==0)
    psi = 0.5*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
    
    if(p->j_dir==1)
    psi = 0.5*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);


    // Construct solid heaviside function
    phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i+aa,j+bb,k+cc));
	
    if(-phival_fb > psi)
    H = 1.0;

    else if(-phival_fb < -psi)
    H = 0.0;

    else
    H = 0.5*(1.0 + -phival_fb/psi + (1.0/PI)*sin((PI*-phival_fb)/psi));
	
    return H;
}



