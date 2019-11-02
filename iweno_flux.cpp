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

#include"iweno_flux.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

iweno_flux::iweno_flux(lexer *p)
			:tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(1.0e-6),deltin (1.0/p->DXM)
{
    if(p->B269==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2(p);
        
        if(p->D11==3)
        pflux = new flux_face_QOU(p);
    }
    
    if(p->B269>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans(p);
        
        if(p->D11==3)
        pflux = new flux_face_FOU_vrans(p);
    }
}

iweno_flux::~iweno_flux()
{
}

void iweno_flux::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    wenoloop1(p,a,b,ipol,uvel,vvel,wvel);

    if(ipol==2)
    wenoloop2(p,a,b,ipol,uvel,vvel,wvel);

    if(ipol==3)
    wenoloop3(p,a,b,ipol,uvel,vvel,wvel);

    if(ipol==4)
    wenoloop4(p,a,b,ipol,uvel,vvel,wvel);

    if(ipol==5)
    wenoloop4(p,a,b,ipol,uvel,vvel,wvel);
}

void iweno_flux::wenoloop1(lexer *p, fdm *a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
	count=0;

	ULOOP
	{	

        pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

			
			i-=1;
			if(ivel1>=0.0)
			{
			a1=1.0;
			iqmin(a,b,uvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(ivel1<0.0)
			{
			a1=0.0;
			iqmax(a,b,uvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			i+=1;
			
			if(ivel2>=0.0)
			{
			a2=1.0;
			iqmin(a,b,uvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(ivel2<0.0)
			{
			a2=0.0;
			iqmax(a,b,uvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
			
			aij_x(p,a,b,a->F);


	
			j-=1;
			if(jvel1>=0.0)
			{
			b1=1.0;
			jqmin(a,b,vvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(jvel1<0.0)
			{
			b1=0.0;
			jqmax(a,b,vvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			j+=1;
			
			if(jvel2>=0.0)
			{
			b2=1.0;
			jqmin(a,b,vvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(jvel2<0.0)
			{
			b2=0.0;
			jqmax(a,b,vvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
			
			aij_y(p,a,b,a->F);


		
			k-=1;
			if(kvel1>=0.0)
			{
			c1=1.0;
			kqmin(a,b,wvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(kvel1<0.0)
			{
			c1=0.0;
			kqmax(a,b,wvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			k+=1;
			
			if(kvel2>=0.0)
			{
			c2=1.0;
			kqmin(a,b,wvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(kvel2<0.0)
			{
			c2=0.0;
			kqmax(a,b,wvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
		
			aij_z(p,a,b,a->F);
			
	++count;
	}
}

void iweno_flux::wenoloop2(lexer *p, fdm *a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
	count=0;

	VLOOP
	{

        pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
			
			i-=1;
			if(ivel1>=0.0)
			{
			a1=1.0;
			iqmin(a,b,uvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(ivel1<0.0)
			{
			a1=0.0;
			iqmax(a,b,uvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			i+=1;
			
			if(ivel2>=0.0)
			{
			a2=1.0;
			iqmin(a,b,uvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(ivel2<0.0)
			{
			a2=0.0;
			iqmax(a,b,uvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
			
			aij_x(p,a,b,a->G);


			
			j-=1;
			if(jvel1>=0.0)
			{
			b1=1.0;
			jqmin(a,b,vvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(jvel1<0.0)
			{
			b1=0.0;
			jqmax(a,b,vvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			j+=1;
			
			if(jvel2>=0.0)
			{
			b2=1.0;
			jqmin(a,b,vvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(jvel2<0.0)
			{
			b2=0.0;
			jqmax(a,b,vvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
			
			aij_y(p,a,b,a->G);


			
			k-=1;
			if(kvel1>=0.0)
			{
			c1=1.0;
			kqmin(a,b,wvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(kvel1<0.0)
			{
			c1=0.0;
			kqmax(a,b,wvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			k+=1;
			
			if(kvel2>=0.0)
			{
			c2=1.0;
			kqmin(a,b,wvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(kvel2<0.0)
			{
			c2=0.0;
			kqmax(a,b,wvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
		
			aij_z(p,a,b,a->G);
			
     ++count;
	}
}

void iweno_flux::wenoloop3(lexer *p, fdm *a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
	count=0;

	WLOOP
	{

		pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
			
			i-=1;
			if(ivel1>=0.0)
			{
			a1=1.0;
			iqmin(a,b,uvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(ivel1<0.0)
			{
			a1=0.0;
			iqmax(a,b,uvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			i+=1;
			
			if(ivel2>=0.0)
			{
			a2=1.0;
			iqmin(a,b,uvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(ivel2<0.0)
			{
			a2=0.0;
			iqmax(a,b,uvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
			
			aij_x(p,a,b,a->H);



			j-=1;
			if(jvel1>=0.0)
			{
			b1=1.0;
			jqmin(a,b,vvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(jvel1<0.0)
			{
			b1=0.0;
			jqmax(a,b,vvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			j+=1;
			
			if(jvel2>=0.0)
			{
			b2=1.0;
			jqmin(a,b,vvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(jvel2<0.0)
			{
			b2=0.0;
			jqmax(a,b,vvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
			
			aij_y(p,a,b,a->H);


		
			k-=1;
			if(kvel1>=0.0)
			{
			c1=1.0;
			kqmin(a,b,wvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(kvel1<0.0)
			{
			c1=0.0;
			kqmax(a,b,wvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			k+=1;
			
			if(kvel2>=0.0)
			{
			c2=1.0;
			kqmin(a,b,wvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(kvel2<0.0)
			{
			c2=0.0;
			kqmax(a,b,wvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
		
			aij_z(p,a,b,a->H);
     ++count;
	}
}

void iweno_flux::wenoloop4(lexer *p, fdm *a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
	count=0;

	LOOP
	{

		pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
        pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
        pflux->w_flux(a,ipol,wvel,kvel1,kvel2);
			
			i-=1;
			if(ivel1>=0.0)
			{
			a1=1.0;
			iqmin(a,b,uvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(ivel1<0.0)
			{
			a1=0.0;
			iqmax(a,b,uvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			i+=1;
			
			if(ivel2>=0.0)
			{
			a2=1.0;
			iqmin(a,b,uvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(ivel2<0.0)
			{
			a2=0.0;
			iqmax(a,b,uvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
			
			aij_x(p,a,b,a->L);


	
			j-=1;
			if(jvel1>=0.0)
			{
			b1=1.0;
			jqmin(a,b,vvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(jvel1<0.0)
			{
			b1=0.0;
			jqmax(a,b,vvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			j+=1;
			
			if(jvel2>=0.0)
			{
			b2=1.0;
			jqmin(a,b,vvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(jvel2<0.0)
			{
			b2=0.0;
			jqmax(a,b,vvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
			
			aij_y(p,a,b,a->L);


		
			k-=1;
			if(kvel1>=0.0)
			{
			c1=1.0;
			kqmin(a,b,wvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}

			if(kvel1<0.0)
			{
			c1=0.0;
			kqmax(a,b,wvel,ipol);
			is_1(b);
			alpha_calc_1();
			weight_calc_1();
			}
			k+=1;
			
			if(kvel2>=0.0)
			{
			c2=1.0;
			kqmin(a,b,wvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}

			if(kvel2<0.0)
			{
			c2=0.0;
			kqmax(a,b,wvel,ipol);
			is_2(b);
			alpha_calc_2();
			weight_calc_2();
			}
		
			aij_z(p,a,b,a->L);	

     ++count;
	}
}

void iweno_flux::aij_x(lexer* p, fdm* a, field &b, field &F)
{
	F(i,j,k) -= (-((1.0/3.0)*w1_1*a1*ivel1)*b(i-3,j,k) 
		
			+ ((1.0/3.0)*w1_2*a2*ivel2
			+ (7.0/6.0)*w1_1*a1*ivel1 + (1.0/6.0)*w2_1*a1*ivel1
			+ (1.0/6.0)*w3_1*(1.0-a1)*ivel1)*b(i-2,j,k) 
		
			+(-(1.0/6.0)*w3_2*a2*ivel2
			- (7.0/6.0)*w1_2*(1.0-a2)*ivel2 - (1.0/6.0)*w2_2*(1.0-a2)*ivel2
			- (1.0/3.0)*w1_1*(1.0-a1)*ivel1)*b(i+2,j,k)
		
			+ ((1.0/3.0)*w1_2*(1.0-a2)*ivel2)*b(i+3,j,k))/p->DXM;

	
	a->M.p[count] = ((11.0/6.0)*w1_2*a2*ivel2 + (5.0/6.0)*w2_2*a2*ivel2 + (1.0/3.0)*w3_2*a2*ivel2
					+ (1.0/3.0)*w2_2*(1.0-a2)*ivel2 + (5.0/6.0)*w3_2*(1.0-a2)*ivel2
							
					- (1.0/3.0)*w2_1*a1*ivel1 - (5.0/6.0)*w3_1*a1*ivel1
					- (11.0/6.0)*w1_1*(1.0-a1)*ivel1 - (5.0/6.0)*w2_1*(1.0-a1)*ivel1 - (1.0/3.0)*w3_1*(1.0-a1)*ivel1) / p->DXM;
					 
	a->M.s[count] = (-(7.0/6.0)*w1_2*a2*ivel2 - (1.0/6.0)*w2_2*a2*ivel2 
					- (1.0/6.0)*w3_2*(1.0-a2)*ivel2 
		 
					- (11.0/6.0)*w1_1*a1*ivel1 - (5.0/6.0)*w2_1*a1*ivel1 - (1.0/3.0)*w3_1*a1*ivel1
					- (1.0/3.0)*w2_1*(1.0-a1)*ivel1 - (5.0/6.0)*w3_1*(1.0-a1)*ivel1) / p->DXM;
							
	a->M.n[count] =  ((1.0/3.0)*w2_2*a2*ivel2 + (5.0/6.0)*w3_2*a2*ivel2
					+ (11.0/6.0)*w1_2*(1.0-a2)*ivel2 + (5.0/6.0)*w2_2*(1.0-a2)*ivel2 + (1.0/3.0)*w3_2*(1.0-a2)*ivel2
							
					+ (1.0/6.0)*w3_1*a1*ivel1 
					+ (7.0/6.0)*w1_1*(1.0-a1)*ivel1 + (1.0/6.0)*w2_1*(1.0-a1)*ivel1) / p->DXM;
}


void iweno_flux::aij_y(lexer* p, fdm* a, field &b, field &F)
{
	F(i,j,k) -= (-((1.0/3.0)*w1_1*b1*jvel1)*b(i,j-3,k)
	
				+((1.0/3.0)*w1_2*b2*jvel2
				+ (7.0/6.0)*w1_1*b1*jvel1 + (1.0/6.0)*w2_1*b1*jvel1
				+ (1.0/6.0)*w3_1*(1.0-b1)*jvel1)*b(i,j-2,k)
				
				+(-(1.0/6.0)*w3_2*b2*jvel2
				- (7.0/6.0)*w1_2*(1.0-b2)*jvel2 - (1.0/6.0)*w2_2*(1.0-b2)*jvel2
				- (1.0/3.0)*w1_1*(1.0-b1)*jvel1)*b(i,j+2,k)
				
				+((1.0/3.0)*w1_2*(1.0-b2)*jvel2)*b(i,j+3,k))/p->DXM;
				
	
	a->M.p[count]   += ((11.0/6.0)*w1_2*b2*jvel2 + (5.0/6.0)*w2_2*b2*jvel2 + (1.0/3.0)*w3_2*b2*jvel2
					+ (1.0/3.0)*w2_2*(1.0-b2)*jvel2 + (5.0/6.0)*w3_2*(1.0-b2)*jvel2
							
					- (1.0/3.0)*w2_1*b1*jvel1 - (5.0/6.0)*w3_1*b1*jvel1
					- (11.0/6.0)*w1_1*(1.0-b1)*jvel1 - (5.0/6.0)*w2_1*(1.0-b1)*jvel1 - (1.0/3.0)*w3_1*(1.0-b1)*jvel1) / p->DXM;
					 
	a->M.e[count]  = (-(7.0/6.0)*w1_2*b2*jvel2 - (1.0/6.0)*w2_2*b2*jvel2 
					- (1.0/6.0)*w3_2*(1.0-b2)*jvel2 
		 
					- (11.0/6.0)*w1_1*b1*jvel1 - (5.0/6.0)*w2_1*b1*jvel1 - (1.0/3.0)*w3_1*b1*jvel1 
					- (1.0/3.0)*w2_1*(1.0-b1)*jvel1 - (5.0/6.0)*w3_1*(1.0-b1)*jvel1) / p->DXM;
							
	a->M.w[count]  = ((1.0/3.0)*w2_2*b2*jvel2 + (5.0/6.0)*w3_2*b2*jvel2
					+ (11.0/6.0)*w1_2*(1.0-b2)*jvel2 + (5.0/6.0)*w2_2*(1.0-b2)*jvel2 + (1.0/3.0)*w3_2*(1.0-b2)*jvel2
							
					+ (1.0/6.0)*w3_1*b1*jvel1 
					+ (7.0/6.0)*w1_1*(1.0-b1)*jvel1 + (1.0/6.0)*w2_1*(1.0-b1)*jvel1) / p->DXM;
}


void iweno_flux::aij_z(lexer* p, fdm* a, field &b, field &F)
{
	F(i,j,k) -= (-((1.0/3.0)*w1_1*c1*kvel1)*b(i,j,k-3)
	
				+((1.0/3.0)*w1_2*c2*kvel2
				+ (7.0/6.0)*w1_1*c1*kvel1 + (1.0/6.0)*w2_1*c1*kvel1
				+ (1.0/6.0)*w3_1*(1.0-c1)*kvel1)*b(i,j,k-2)
				
				+(-(1.0/6.0)*w3_2*c2*kvel2
				- (7.0/6.0)*w1_2*(1.0-c2)*kvel2 - (1.0/6.0)*w2_2*(1.0-c2)*kvel2
				- (1.0/3.0)*w1_1*(1.0-c1)*kvel1)*b(i,j,k+2)
				
				+((1.0/3.0)*w1_2*(1.0-c2)*kvel2)*b(i,j,k+3))/p->DXM;
	
	a->M.p[count]  += ((11.0/6.0)*w1_2*c2*kvel2 + (5.0/6.0)*w2_2*c2*kvel2 + (1.0/3.0)*w3_2*c2*kvel2
					+ (1.0/3.0)*w2_2*(1.0-c2)*kvel2 + (5.0/6.0)*w3_2*(1.0-c2)*kvel2
							
					- (1.0/3.0)*w2_1*c1*kvel1 - (5.0/6.0)*w3_1*c1*kvel1
					- (11.0/6.0)*w1_1*(1.0-c1)*kvel1 - (5.0/6.0)*w2_1*(1.0-c1)*kvel1 - (1.0/3.0)*w3_1*(1.0-c1)*kvel1) / p->DXM;
					 
	a->M.b[count] =  (-(7.0/6.0)*w1_2*c2*kvel2 - (1.0/6.0)*w2_2*c2*kvel2 
					- (1.0/6.0)*w3_2*(1.0-c2)*kvel2 
		 
					- (11.0/6.0)*w1_1*c1*kvel1 - (5.0/6.0)*w2_1*c1*kvel1 - (1.0/3.0)*w3_1*c1*kvel1 
					- (1.0/3.0)*w2_1*(1.0-c1)*kvel1 - (5.0/6.0)*w3_1*(1.0-c1)*kvel1) / p->DXM;
							
	a->M.t[count] =  ((1.0/3.0)*w2_2*c2*kvel2 + (5.0/6.0)*w3_2*c2*kvel2
					+ (11.0/6.0)*w1_2*(1.0-c2)*kvel2 + (5.0/6.0)*w2_2*(1.0-c2)*kvel2 + (1.0/3.0)*w3_2*(1.0-c2)*kvel2
							
					+ (1.0/6.0)*w3_1*c1*kvel1 
					+ (7.0/6.0)*w1_1*(1.0-c1)*kvel1 + (1.0/6.0)*w2_1*(1.0-c1)*kvel1) / p->DXM;
}

void iweno_flux::iqmin(fdm* a, field& f, field& uvel, int ipol)
{	
	q1 = f(i-2,j,k);
	q2 = f(i-1,j,k);
	q3 = f(i,j,k);
	q4 = f(i+1,j,k);
	q5 = f(i+2,j,k);
}

void iweno_flux::jqmin(fdm* a, field& f, field& vvel, int ipol)
{
	q1 = f(i,j-2,k);
	q2 = f(i,j-1,k);
	q3 = f(i,j,k);
	q4 = f(i,j+1,k);
	q5 = f(i,j+2,k);
}

void iweno_flux::kqmin(fdm* a, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k-2);
	q2 = f(i,j,k-1);
	q3 = f(i,j,k);
	q4 = f(i,j,k+1);
	q5 = f(i,j,k+2);
}

void iweno_flux::iqmax(fdm* a, field& f, field& uvel, int ipol)
{
	q1 = f(i+3,j,k);
	q2 = f(i+2,j,k);
	q3 = f(i+1,j,k);
	q4 = f(i,j,k);
	q5 = f(i-1,j,k);
}

void iweno_flux::jqmax(fdm* a, field& f, field& vvel, int ipol)
{
	q1 = f(i,j+3,k);
	q2 = f(i,j+2,k);
	q3 = f(i,j+1,k);
	q4 = f(i,j,k);
	q5 = f(i,j-1,k);
}

void iweno_flux::kqmax(fdm* a, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k+3);
	q2 = f(i,j,k+2);
	q3 = f(i,j,k+1);
	q4 = f(i,j,k);
	q5 = f(i,j,k-1);
}

void iweno_flux::is_1(field& b)
{
	is1_1 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
	is2_1 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
	is3_1 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
}

void iweno_flux::is_2(field& b)
{
	is1_2 = tttw*pow(q1-2.0*q2+q3, 2.0) + fourth*pow(q1-4.0*q2+3.0*q3, 2.0);
	is2_2 = tttw*pow(q2-2.0*q3+q4, 2.0) + fourth*pow(q2-q4, 2.0);
	is3_2 = tttw*pow(q3-2.0*q4+q5, 2.0) + fourth*pow(3.0*q3-4.0*q4+q5, 2.0);
}

void iweno_flux::alpha_calc_1()
{
	alpha1_1=0.1/pow(epsilon+is1_1,2.0);
	alpha2_1=0.6/pow(epsilon+is2_1,2.0);
	alpha3_1=0.3/pow(epsilon+is3_1,2.0);
}

void iweno_flux::alpha_calc_2()
{
	alpha1_2=0.1/pow(epsilon+is1_2,2.0);
	alpha2_2=0.6/pow(epsilon+is2_2,2.0);
	alpha3_2=0.3/pow(epsilon+is3_2,2.0);
}

void iweno_flux::weight_calc_1()
{
	w1_1=alpha1_1/(alpha1_1+alpha2_1+alpha3_1);
	w2_1=alpha2_1/(alpha1_1+alpha2_1+alpha3_1);
	w3_1=alpha3_1/(alpha1_1+alpha2_1+alpha3_1);
}

void iweno_flux::weight_calc_2()
{
	w1_2=alpha1_2/(alpha1_2+alpha2_2+alpha3_2);
	w2_2=alpha2_2/(alpha1_2+alpha2_2+alpha3_2);
	w3_2=alpha3_2/(alpha1_2+alpha2_2+alpha3_2);
}
