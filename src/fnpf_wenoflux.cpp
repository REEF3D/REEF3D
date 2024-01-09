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

#include"fnpf_wenoflux.h"
#include"lexer.h"
#include"slice.h"
#include"vec.h"
#include"fnpf_discrete_weights.h"

fnpf_wenoflux::fnpf_wenoflux(lexer* p) :  weno_nug_func(p)
{
    p->Darray(ckz,p->knoz+1+4*marge,5);
    
    fnpf_discrete_weights dw(p);

    dw.ck_weights(p, ckz, p->ZN, p->knoz+1, 1, 4, 6);
    
    uf=vf=wf=0;
}

fnpf_wenoflux::~fnpf_wenoflux()
{
}

double fnpf_wenoflux::sx(lexer *p, slice &f, slice &Fifsf)
{
    grad=0.0;
        
        ivel1 = (Fifsf(i,j) - Fifsf(i-1,j))/(p->DXP[IM1]);
        ivel2 = (Fifsf(i+1,j) - Fifsf(i,j))/(p->DXP[IP]);
        
        
        //ivel1 = (-Fifsf(i+1,j) + 27.0*Fifsf(i,j) - 27.0*Fifsf(i-1,j) + Fifsf(i-2,j))/(-p->XP[IP1] + 27.0*p->XP[IP] - 27.0*p->XP[IM1] + p->XP[IM2]);
        //ivel2 = (-Fifsf(i+2,j) + 27.0*Fifsf(i+1,j) - 27.0*Fifsf(i,j) + Fifsf(i-1,j))/(-p->XP[IP2] + 27.0*p->XP[IP1] - 27.0*p->XP[IP] + p->XP[IM1]);

		
		i-=1;
		fu1 = ffx(p,f,ivel1);
		i+=1;
		
		fu2 = ffx(p,f,ivel2);
		
		grad = ((ivel2*fu2-ivel1*fu1)/p->DXN[IP]);
        
    return grad;
}

double fnpf_wenoflux::sy(lexer *p, slice &f, slice &Fifsf)
{
    grad=0.0;
        
        jvel1 = (Fifsf(i,j) - Fifsf(i,j-1))/(p->DYP[JM1]);
        jvel2 = (Fifsf(i,j+1) - Fifsf(i,j))/(p->DYP[JP]);
		
		j-=1;
		fv1 = ffy(p,f,jvel1);
		j+=1;
		
		fv2 = ffy(p,f,jvel2);
		
		grad = ((jvel2*fv2-jvel1*fv1)/p->DYN[JP]);
			  
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

double fnpf_wenoflux::ffy(lexer *p, slice &f, double advec)
{
    grad = 0.0;

	if(advec>0.0)
	{
	jqmin(p,f);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(advec<0.0)
	{
	jqmax(p,f);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
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
	q2 = f(i,j);
	q3 = f(i,j+1);
	q4 = f(i,j+2);
	q5 = f(i,j+3);
}

