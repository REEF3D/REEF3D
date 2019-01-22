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

fnpf_wenoflux::fnpf_wenoflux(lexer* p) :  weno_nug_func(p)
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

    return grad;
}

double fnpf_wenoflux::fy(lexer *p, field &f, double jvel1, double jvel2)
{
    grad=0.0;

    return grad;
}

double fnpf_wenoflux::fz(lexer *p, field &f, double kvel1, double kvel2)
{
    grad=0.0;
    
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
    
    
    return grad;   
}

double fnpf_wenoflux::sz(lexer *p, double *f)
{

    grad = (ckz[p->knoz+marge][4]*f[FIJK] + ckz[p->knoz+marge][3]*f[FIJKm1] + ckz[p->knoz+marge][2]*f[FIJKm2] 
          + ckz[p->knoz+marge][1]*f[FIJKm3] + ckz[p->knoz+marge][0]*f[FIJKm4]);
    
    return grad;   
}



double fnpf_wenoflux::ffx(lexer *p, slice &f, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	iqmin(p,f);
	is_min_x();
	weight_min_x();

	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	iqmax(p,f);
	is_max_x();
	weight_max_x();
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
	return grad;
}


void fnpf_wenoflux::iqmin(lexer *p, slice &f)
{	
	q1 = f(i-2,j);
	q2 = f(i-1,j);
	q3 = f(i,j);
	q4 = f(i+1,j);
	q5 = f(i+2,j);
}

void fnpf_wenoflux::jqmin(lexer *p, slice &f)
{
	q1 = f(i,j-2);
	q2 = f(i,j-1);
	q3 = f(i,j);
	q4 = f(i,j+1);
	q5 = f(i,j+2);
}

void fnpf_wenoflux::iqmax(lexer *p, slice &f)
{
    q1 = f(i-1,j);
	q2 = f(i,j);
	q3 = f(i+1,j);
	q4 = f(i+2,j);
	q5 = f(i+3,j);
}

void fnpf_wenoflux::jqmax(lexer *p, slice &f)
{
	q1 = f(i,j-1);
	q2 = f(i,j,k);
	q3 = f(i,j+1);
	q4 = f(i,j+2);
	q5 = f(i,j+3);
}

