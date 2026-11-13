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
Authors: Tobias Martin, Ahmet Soydan, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"heaviside.h"
#include"forcing_share.h"

// the weight H of a face of the body, at most what the solids (hs) and the bodies before it (hb)
// leave of 1: the forcing terms of a face add up to a weight of 1 at most (overlaps)
static inline double Hlimit(double H, double hb, double hs)
{
    const double hm = 1.0 - hb - hs;
    
    if(H>hm)
    H = (hm>0.0) ? hm : 0.0;
    
    return H;
}

// X 14 2: the normal of the hull at face d (0 u, 1 v, 2 w) of cell (i,j,k), the gradient of fb as the mean
// of the central differences of the two cells of the face (X 14 1 takes that of cell (i,j,k) alone, half a
// cell off the face: a body symmetric about a plane through faces of the grid is forced asymmetrically)
void sixdof_obj_cfd::face_normal(lexer *p, fdm *a, int d, double &nx, double &ny, double &nz)
{
    const int ii = i + (d==0), jj = j + (d==1), kk = k + (d==2);
    
    nx = -0.5*((a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]) + (a->fb(ii+1,jj,kk) - a->fb(ii-1,jj,kk))/(2.0*p->DXN[ii+marge]));
    ny = -0.5*((a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]) + (a->fb(ii,jj+1,kk) - a->fb(ii,jj-1,kk))/(2.0*p->DYN[jj+marge]));
    nz = -0.5*((a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]) + (a->fb(ii,jj,kk+1) - a->fb(ii,jj,kk-1))/(2.0*p->DZN[kk+marge]));
}

void sixdof_obj_cfd::update_forcing(lexer *p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, field &fx, field &fy, field &fz,int iter)
{
    // X 198: the rates of the discrete recoil correction of a deforming body at the new pose
    if(dc_on_ && amr_loads==nullptr)
    momentum_deform(p,a,pgc);

    
// Calculate forcing fields
    double H,Ht, uf, vf, wf;
	double nx, ny, nz,norm ;
	double psi, phival_fb;
    double dirac;
    
    H=Ht=0.0;
  
    if(p->X14==0)
    {
    ULOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;

        uf = u_fb(0) + u_fb(4)*(p->pos1_z() - c_(2)) - u_fb(5)*(p->pos1_y() - c_(1));
        
        if(uf!=uf)
        {
        cout<<"UF "<<uf<<" "<<u_fb(0)<<" "<<u_fb(4)<<" "<<(p->pos1_z() - c_(2))<<" "<<u_fb(5)<<" "<<(p->pos1_y() - c_(1))<<endl;
        }
        
        H = Hsolidface(p,a,1,0,0);
        
        // deforming body (X 195): velocity of the deformation (where the body forces the face)
        if(fish_!=nullptr && H>0.0)
        uf += deform_vel(p->pos1_x(),p->pos1_y(),p->pos1_z())(0);
        H = Hlimit(H,a->fbh1(i,j,k),sfh[0]!=nullptr ? (*sfh[0])(i,j,k) : 0.0);
        
        if(H!=H)
        cout<<"H "<<uf<<endl;
       
        fx(i,j,k) += H*(uf - uvel(i,j,k))/(alpha[iter]*p->dt);   
        a->fbh1(i,j,k) = min(a->fbh1(i,j,k) + H, 1.0); 
    }
    VLOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;

        vf = u_fb(1) + u_fb(5)*(p->pos2_x() - c_(0)) - u_fb(3)*(p->pos2_z() - c_(2));
        
        if(vf!=vf)
        cout<<"VF "<<vf<<endl;
        
        H = Hsolidface(p,a,0,1,0);
        
        // deforming body (X 195): velocity of the deformation (where the body forces the face)
        if(fish_!=nullptr && H>0.0)
        vf += deform_vel(p->pos2_x(),p->pos2_y(),p->pos2_z())(1);
        H = Hlimit(H,a->fbh2(i,j,k),sfh[1]!=nullptr ? (*sfh[1])(i,j,k) : 0.0);
        
        if(H!=H)
        cout<<"H "<<uf<<endl;
       
        fy(i,j,k) += H*(vf - vvel(i,j,k))/(alpha[iter]*p->dt);
        a->fbh2(i,j,k) = min(a->fbh2(i,j,k) + H, 1.0); 
    }
    WLOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;

        wf = u_fb(2) + u_fb(3)*(p->pos3_y() - c_(1)) - u_fb(4)*(p->pos3_x() - c_(0));
        
        if(wf!=wf)
        cout<<"WF "<<wf<<endl;
        
        H = Hsolidface(p,a,0,0,1);
        
        // deforming body (X 195): velocity of the deformation (where the body forces the face)
        if(fish_!=nullptr && H>0.0)
        wf += deform_vel(p->pos3_x(),p->pos3_y(),p->pos3_z())(2);
        H = Hlimit(H,a->fbh3(i,j,k),sfh[2]!=nullptr ? (*sfh[2])(i,j,k) : 0.0);
        
        if(H!=H)
        cout<<"H "<<uf<<endl;

        fz(i,j,k) += H*(wf - wvel(i,j,k))/(alpha[iter]*p->dt);
        a->fbh3(i,j,k) = min(a->fbh3(i,j,k) + H, 1.0); 
    }
    LOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;

        H = Hsolidface(p,a,0,0,0);
        a->fbh4(i,j,k) = min(a->fbh4(i,j,k) + H, 1.0); 
    }
    	
    psi = 1.1*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if (p->j_dir==0)
    psi = 1.1*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 

    LOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        {
            a->fbh5(i,j,k) = 1.0;
            continue;
        }

        dirac = 0.0;
        if(fabs(a->fb(i,j,k))<psi)
        dirac = (0.5/psi)*(1.0 + cos((PI*(a->fb(i,j,k)))/psi));
        
        a->fbh5(i,j,k) = 1.0-MIN(dirac,1.0);
    }
    
    }


// Construct solid heaviside function	
    if(p->X14>=1)
    {
        
    ULOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;

        uf = u_fb(0) + u_fb(4)*(p->pos1_z() - c_(2)) - u_fb(5)*(p->pos1_y() - c_(1));
        
        if(uf!=uf)
        cout<<"UF "<<uf<<endl;
        
		// Normal vectors calculation (X 14 2: at the face, the mean of the cells on both sides)
		nx = -(a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]);
        if(p->X14==2)
        face_normal(p,a,0,nx,ny,nz);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;


		H = Hsolidface(p,a,1,0,0);
        
        // deforming body (X 195): velocity of the deformation (where the body forces the face)
        if(fish_!=nullptr && H>0.0)
        uf += deform_vel(p->pos1_x(),p->pos1_y(),p->pos1_z())(0);
        H = Hlimit(H,a->fbh1(i,j,k),sfh[0]!=nullptr ? (*sfh[0])(i,j,k) : 0.0);
	    Ht = Hsolidface_t(p,a,1,0,0);
	
    
		// Level set function
		phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i+1,j,k));	
		
		// Construct the field around the solid body to adjust the tangential velocity and calculate forcing
		// X 46 1: the share continuous across the surface (forcing_share.h)
		if(p->X46>=1)
		fx(i,j,k) += forcing_weight_smooth(H,fabs(nx),phival_fb,forcing_psi(p,i,j,k))*(uf - uvel(i,j,k))/(alpha[iter]*p->dt);

		if(p->X46<1 && phival_fb<=0.0)
		fx(i,j,k) += H*(uf - uvel(i,j,k))/(alpha[iter]*p->dt); 

		if(p->X46<1 && phival_fb>0.0)
		fx(i,j,k) +=   fabs(nx)*H*(uf - uvel(i,j,k))/(alpha[iter]*p->dt);
            

        a->fbh1(i,j,k) = min(a->fbh1(i,j,k) + H, 1.0); 
    }
    
    if(p->j_dir==1)
    VLOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;

        vf = u_fb(1) + u_fb(5)*(p->pos2_x() - c_(0)) - u_fb(3)*(p->pos2_z() - c_(2));
        
        if(vf!=vf)
        cout<<"VF "<<vf<<endl;
    
		// Normal vectors calculation (X 14 2: at the face, the mean of the cells on both sides)
		nx = -(a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]);
        if(p->X14==2)
        face_normal(p,a,1,nx,ny,nz);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;

        
         H = Hsolidface(p,a,0,1,0);
        
        // deforming body (X 195): velocity of the deformation (where the body forces the face)
        if(fish_!=nullptr && H>0.0)
        vf += deform_vel(p->pos2_x(),p->pos2_y(),p->pos2_z())(1);
        H = Hlimit(H,a->fbh2(i,j,k),sfh[1]!=nullptr ? (*sfh[1])(i,j,k) : 0.0);
		Ht = Hsolidface_t(p,a,0,1,0);
		
      
		//Level set function
		phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i,j+1,k));
	  
		//Construct the field around the solid body to adjust the tangential velocity and calculate forcing
		// X 46 1: the share continuous across the surface (forcing_share.h)
		if(p->X46>=1)
		fy(i,j,k) += forcing_weight_smooth(H,fabs(ny),phival_fb,forcing_psi(p,i,j,k))*(vf - vvel(i,j,k))/(alpha[iter]*p->dt);

		if(p->X46<1 && phival_fb<=0.0)
		fy(i,j,k) += H*(vf - vvel(i,j,k))/(alpha[iter]*p->dt); 

		if(p->X46<1 && phival_fb>0.0)
		fy(i,j,k) +=   fabs(ny)*H*(vf - vvel(i,j,k))/(alpha[iter]*p->dt);


        a->fbh2(i,j,k) = min(a->fbh2(i,j,k) + H , 1.0); 
    }
	
    WLOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;

        wf = u_fb(2) + u_fb(3)*(p->pos3_y() - c_(1)) - u_fb(4)*(p->pos3_x() - c_(0));
        
        if(wf!=wf)
        cout<<"WF "<<wf<<endl;
        
		// Normal vectors calculation (X 14 2: at the face, the mean of the cells on both sides)
		nx = -(a->fb(i+1,j,k) - a->fb(i-1,j,k))/(2.0*p->DXN[IP]);
		ny = -(a->fb(i,j+1,k) - a->fb(i,j-1,k))/(2.0*p->DYN[JP]);
		nz = -(a->fb(i,j,k+1) - a->fb(i,j,k-1))/(2.0*p->DZN[KP]);
        if(p->X14==2)
        face_normal(p,a,2,nx,ny,nz);

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;

        
         H = Hsolidface(p,a,0,0,1);
        
        // deforming body (X 195): velocity of the deformation (where the body forces the face)
        if(fish_!=nullptr && H>0.0)
        wf += deform_vel(p->pos3_x(),p->pos3_y(),p->pos3_z())(2);
        H = Hlimit(H,a->fbh3(i,j,k),sfh[2]!=nullptr ? (*sfh[2])(i,j,k) : 0.0);
		Ht = Hsolidface_t(p,a,0,0,1);


		// Level set function
		phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i,j,k+1));
		
		// Construct the field around the solid body to adjust the tangential velocity and calculate forcing

		// X 46 1: the share continuous across the surface (forcing_share.h)
		if(p->X46>=1)
		fz(i,j,k) += forcing_weight_smooth(H,fabs(nz),phival_fb,forcing_psi(p,i,j,k))*(wf - wvel(i,j,k))/(alpha[iter]*p->dt);

		if(p->X46<1 && phival_fb<=0.0)
		fz(i,j,k) += H*(wf - wvel(i,j,k))/(alpha[iter]*p->dt); 

		if(p->X46<1 && phival_fb>0.0)
		fz(i,j,k) +=   fabs(nz)*H*(wf - wvel(i,j,k))/(alpha[iter]*p->dt);


        a->fbh3(i,j,k) = min(a->fbh3(i,j,k) + H , 1.0); 
    }
    
    LOOP
    {
        if(amr_cover!=nullptr && amr_cover[IJK])
        continue;

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
        if(amr_cover!=nullptr && amr_cover[IJK])
        {
            a->fbh5(i,j,k) = 1.0;
            continue;
        }

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

double sixdof_obj_cfd::Hsolidface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi, H, phival_fb;
	
    psi = p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if(p->j_dir==0)
    psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 


    // Construct solid heaviside function
    phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i+aa,j+bb,k+cc));
	
    H = heaviside(-phival_fb,psi);
    
    return H;
}

double sixdof_obj_cfd::Hsolidface_t(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double psi, H, phival_fb;
	

    if(p->j_dir==0)
    psi = 0.5*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
    
    if(p->j_dir==1)
    psi = 0.5*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);


    // Construct solid heaviside function
    phival_fb = 0.5*(a->fb(i,j,k) + a->fb(i+aa,j+bb,k+cc));
	
    H = heaviside(-phival_fb,psi);
	
    return H;
}



