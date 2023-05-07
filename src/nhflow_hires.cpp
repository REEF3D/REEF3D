/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"nhflow_hires.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"nhflow_flux_face_cds2.h"
#include"nhflow_vanleer.h"

nhflow_hires::nhflow_hires (lexer *p, int limiter) 
{
    pflux = new nhflow_flux_face_cds2(p);
    
	//if(limiter==10)
	//plim = new minmod(p);
	
	//if(limiter==11)
	plim = new nhflow_vanleer(p);
	/*
	if(limiter==12)
	plim = new umist(p);
	
	if(limiter==13)
	plim = new vanalbada(p);
	
	if(limiter==14)
	plim = new superbee(p);
	
	if(limiter==15)
	plim = new smart(p);
	
	if(limiter==16)
	plim = new limo3(p);
	
	if(limiter==42)
	plim = new tvdvof(p);*/
}

nhflow_hires::~nhflow_hires()
{

}

void nhflow_hires::start(lexer* p, fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W)
{ 	
        if(ipol==1)
        LOOP
        d->F[IJK]+=aij(p,d,F,1,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==2)
        LOOP
        d->G[IJK]+=aij(p,d,F,2,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==3)
        LOOP
        d->H[IJK]+=aij(p,d,F,3,U,V,W,p->DXN,p->DYN,p->DZN);

        if(ipol==4)
        LOOP
        d->L[IJK]+=aij(p,d,F,4,U,V,W,p->DXN,p->DYN,p->DZN);
}

double nhflow_hires::aij(lexer* p,fdm_nhf* d, double *F, int ipol, double *U, double *V, double *W, double *DX,double *DY, double *DZ)
{

    udir=vdir=wdir=0.0;
		
        pflux->u_flux(d,ipol,U,ivel1,ivel2);
        pflux->v_flux(d,ipol,V,jvel1,jvel2);
        pflux->w_flux(d,ipol,d->omega,kvel1,kvel2);
		
		// x-dir
        if(0.5*(ivel1+ivel2)>=0.0)
        udir=1.0;

		dx = udir*(ivel2*(F[IJK] + 0.5*plim->iphi(F,0,-1,1,0)*(F[Ip1JK]-F[IJK]))
        
                - ivel1*(F[Im1JK] + 0.5*plim->iphi(F,-1,-2,0,-1)*(F[IJK]-F[Im1JK])))/DX[IM1] 
            
            
             + (1.0-udir)*(ivel2*(F[Ip1JK] - 0.5*plim->iphi(F,1,0,2,1)*(F[Ip2JK]-F[Ip1JK]))
          
             -            ivel1*(F[IJK] - 0.5*plim->iphi(F,0,-1,1,0)*(F[Ip1JK]-F[IJK])))/DX[IP]; 
             

		// y-dir
        if(0.5*(jvel1+jvel2)>=0.0)
        vdir=1.0;

		dy = vdir*(jvel2*(F[IJK] + 0.5*plim->jphi(F,0,-1,1,0)*(F[IJp1K]-F[IJK]))
        
                - jvel1*(F[IJm1K] + 0.5*plim->jphi(F,-1,-2,0,-1)*(F[IJK]-F[IJm1K])))/DY[JM1] 
        
        
             + (1.0-vdir)*(jvel2*(F[IJp1K] - 0.5*plim->jphi(F,1,0,2,1)*(F[IJp2K]-F[IJp1K]))

             -             jvel1*(F[IJK] - 0.5*plim->jphi(F,0,-1,1,0)*(F[IJK]-F[Ip1JK])))/DY[JP]; 


		// z-dir
        if(0.5*(kvel1+kvel2)>=0.0)
        wdir=1.0;

		dz = wdir*(kvel2*(F[IJK] + 0.5*plim->kphi(F,0,-1,1,0)*(F[IJp1K]-F[IJK]))
        
                -  kvel1*(F[IJm1K] + 0.5*plim->kphi(F,-1,-2,0,-1)*(F[IJK]-F[IJm1K])))/DZ[KM1] 
        
        
            + (1.0-wdir)*(kvel2*(F[IJp1K] - 0.5*plim->kphi(F,1,0,2,1)*(F[IJKp2]-F[IJp1K]))
          
             -      kvel1*(F[IJK] - 0.5*plim->kphi(F,0,-1,1,0)*(F[IJp1K]-F[IJK])))/DZ[KP];
		

		L = -dx-dy-dz;

		return L;
}
