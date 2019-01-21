/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"fnpf_wenoflux.h"
#include"lexer.h"
#include"vec.h"
#include"fnpf_discrete_weights.h"

fnpf_wenoflux::fnpf_wenoflux(lexer* p) :  ddweno_f_nug(p)
{
    p->Darray(ckz,p->knoz+1+4*marge,5);
    
    fnpf_discrete_weights dw(p);

    dw.ck_weights(p, ckz, p->ZN, p->knoz+1, 1, 4, 6);
    
    
}

fnpf_wenoflux::~fnpf_wenoflux()
{
}

double fnpf_wenoflux::fx(lexer *p, field &f, double ivel1, double ivel2)
{
    grad=0.0;
    
    if(0.5*(ivel1+ivel2)>0.0)
    grad=ddwenox(f,1.0);
    
    if(0.5*(ivel1+ivel2)<0.0)
    grad=ddwenox(f,-1.0);
    
    return grad;
}

double fnpf_wenoflux::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    grad=0.0;
    
    if(0.5*(jvel1+jvel2)>0.0)
    grad=ddwenoy(f,1.0);
    
    if(0.5*(jvel1+jvel2)<0.0)
    grad=ddwenoy(f,-1.0);
    
    return grad;
}

double fnpf_wenoflux::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    grad=0.0;
    
    if(0.5*(kvel1+kvel2)>0.0)
    grad=ddwenoz(f,1.0);
    
    if(0.5*(kvel1+kvel2)<0.0)
    grad=ddwenoz(f,-1.0);
    
    return grad;
}

double fnpf_wenoflux::sx(lexer *p, slice &f, double ivel)
{
    grad=0.0;
        
        ivel1 = 0.0;
        ivel2 = 0.0;
		/*
		i-=1;
		fu1 = fx(p,a,b,uvel,ipol,ivel1);
		i+=1;
		
		fu2 = fx(p,a,b,uvel,ipol,ivel2);


		
		
		L =   - ((ivel2*fu2-ivel1*fu1)/DX[IP]) 
		      - ((jvel2*fv2-jvel1*fv1)/DY[JP]) 
			  - ((kvel2*fw2-kvel1*fw1)/DZ[KP]);
			  */
    
    return grad;
}

double fnpf_wenoflux::sy(lexer *p, slice &f, double jvel)
{
    grad=0.0;
    
    if(jvel>0.0)
    grad=dswenoy(f,1.0);
    
    if(jvel<0.0)
    grad=dswenoy(f,-1.0);
    
    return grad;   
}

double fnpf_wenoflux::sz(lexer *p, double *f)
{

    grad = (ckz[p->knoz+marge][4]*f[FIJK] + ckz[p->knoz+marge][3]*f[FIJKm1] + ckz[p->knoz+marge][2]*f[FIJKm2] 
          + ckz[p->knoz+marge][1]*f[FIJKm3] + ckz[p->knoz+marge][0]*f[FIJKm4]);
    
    return grad;   
}
