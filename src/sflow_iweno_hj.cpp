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

#include"sflow_iweno_hj.h"
#include"lexer.h"
#include"fdm2D.h"
#include"sflow_flux_face_FOU.h"
#include"sflow_flux_face_CDS.h"
#include"sflow_flux_face_HJ.h"

sflow_iweno_hj::sflow_iweno_hj(lexer *p)
			:tttw(13.0/12.0),fourth(1.0/4.0),third(1.0/3.0),
			sevsix(7.0/6.0),elvsix(11.0/6.0),sixth(1.0/6.0),fivsix(5.0/6.0),tenth(1.0/10.0),
			sixten(6.0/10.0),treten(3.0/10.0),epsilon(1.0e-6),deltin (1.0/p->DXM)
{
    if(p->A216==1)
    pflux = new sflow_flux_face_FOU(p);
        
    if(p->A216==2)
    pflux = new sflow_flux_face_CDS(p);
    
    if(p->A216==4)
    pflux = new sflow_flux_face_HJ(p);
}

sflow_iweno_hj::~sflow_iweno_hj()
{
}

void sflow_iweno_hj::start(lexer* p, fdm2D* b, slice& f, int ipol, slice& uvel, slice& vvel)
{
    if(ipol==1)
    wenoloop1(p,b,f,ipol,uvel,vvel);

    if(ipol==2)
    wenoloop2(p,b,f,ipol,uvel,vvel);

    if(ipol==4)
    wenoloop4(p,b,f,ipol,uvel,vvel);
}

void sflow_iweno_hj::wenoloop1(lexer *p, fdm2D *b, slice &f, int ipol, slice &uvel, slice &vvel)
{
}

void sflow_iweno_hj::wenoloop2(lexer *p, fdm2D *b, slice &f, int ipol, slice &uvel, slice &vvel)
{
}

void sflow_iweno_hj::wenoloop4(lexer *p, fdm2D *b, slice &f, int ipol, slice &uvel, slice &vvel)
{
	count=0;

	SLICELOOP4
	{
		
        pflux->u_flux(ipol,uvel,iadvec,ivel2);
		pflux->v_flux(ipol,vvel,jadvec,jvel2);

			if(iadvec>=0.0)
			{
			is_south(f);
			alpha_calc();
			weight_calc();
			aij_south(p,b,f,b->L);
			}

			if(iadvec<0.0)
			{
			is_north(f);
			alpha_calc();
			weight_calc();
			aij_north(p,b,f,b->L);
			}


			if(jadvec>=0.0 && p->j_dir==1)
			{
			is_east(f);
			alpha_calc();
			weight_calc();
			aij_east(p,b,f,b->L);
			}

			if(jadvec<0.0 && p->j_dir==1)
			{
			is_west(f);
			alpha_calc();
			weight_calc();
			aij_west(p,b,f,b->L);
			}

     ++count;
	}
}

void sflow_iweno_hj::aij_south(lexer* p, fdm2D *b, slice &f, slice &F)
{
	F(i,j)    += (w1*third)*iadvec*deltin*f(i-3,j)
				 - (w2*sixth + w1*1.5)*iadvec*deltin*f(i-2,j)
				 + (w3*sixth)*iadvec*deltin*f(i+2,j);
				 
	b->M.p[count] = (-w3*0.5 + w2*0.5 + w1*elvsix)*iadvec*deltin;
					 
	b->M.s[count] = (-w3*third -w2 - w1*3.0)*iadvec*deltin;
	b->M.n[count] = (w3 + w2*third)*iadvec*deltin;
}

void sflow_iweno_hj::aij_north(lexer* p, fdm2D *b, slice &f, slice &F)
{
	F(i,j)   += -(w3*sixth)*iadvec*deltin*f(i-2,j)
				-	(-w1*1.5 - w2*sixth)*iadvec*deltin*f(i+2,j)
				-   (w1*third)*iadvec*deltin*f(i+3,j);
					
	b->M.p[count] = (-w1*elvsix - w2*0.5 + w3*0.5)*iadvec*deltin;
					 
	b->M.s[count] = (-w2*third - w3)*iadvec*deltin;
	b->M.n[count] = (w1*3.0 + w2 + w3*third)*iadvec*deltin;
}

void sflow_iweno_hj::aij_east(lexer* p, fdm2D *b, slice &f, slice &F)
{
	F(i,j) 	+= (w1*third)*jadvec*deltin*f(i,j-3)
				-  (w2*sixth + w1*1.5)*jadvec*deltin*f(i,j-2)
				-  (-w3*sixth)*jadvec*deltin*f(i,j+2);
					
	b->M.p[count] += (-w3*0.5 +w2*0.5 + w1*elvsix)*jadvec*deltin;
					 
	b->M.e[count] = (-w3*third -w2 - w1*3.0)*jadvec*deltin;
	b->M.w[count] = (w3 + w2*third)*jadvec*deltin;
}

void sflow_iweno_hj::aij_west(lexer* p, fdm2D *b, slice &f, slice &F)
{
	F(i,j) 	+= -(w3*sixth)*jadvec*deltin*f(i,j-2)
				-	(-w1*1.5 - w2*sixth)*jadvec*deltin*f(i,j+2)
				-	 (w1*third)*jadvec*deltin*f(i,j+3);
				
	b->M.p[count] += (-w1*elvsix - w2*0.5 + w3*0.5)*jadvec*deltin;
					 
	b->M.e[count] = (-w2*third - w3)*jadvec*deltin;
	b->M.w[count] = (w1*3.0 + w2 + w3*third)*jadvec*deltin;
}

void sflow_iweno_hj::is_south(slice &b)
{

	is1 = tttw*pow( ( -b(i-3,j) + 3.0*b(i-2,j) -3.0*b(i-1,j) + b(i,j)),2.0)
		+ fourth*pow( (-b(i-3,j) + 5.0*b(i-2,j) - 7.0*b(i-1,j) + 3.0*b(i,j)),2.0);

	is2 = tttw*pow( ( -b(i-2,j) + 3.0*b(i-1,j) -3.0*b(i,j) + b(i+1,j)),2.0)
		+ fourth*pow( (-b(i-2,j) + b(i-1,j) + b(i,j) - b(i+1,j)),2.0);

	is3 = tttw*pow( (      -b(i-1,j) + 3.0*b(i,j) - 3.0*b(i+1,j) + b(i+2,j)),2.0)
		+ fourth*pow( (-3.0*b(i-1,j) + 7.0*b(i,j) - 5.0*b(i+1,j) + b(i+2,j)),2.0);
}

void sflow_iweno_hj::is_north(slice &b)
{

	is1 = tttw*pow( (b(i+3,j) - 3.0*b(i+2,j) + 3.0*b(i+1,j) - b(i,j)),2.0)
		+ fourth*pow( (b(i+3,j) - 5.0*b(i+2,j) + 7.0*b(i+1,j) - 3.0*b(i,j)),2.0);

	is2 = tttw*pow( (b(i+2,j) - 3.0*b(i+1,j) + 3.0*b(i,j) - b(i-1,j)),2.0)
		+ fourth*pow( (b(i+2,j) - b(i+1,j) - b(i,j) + b(i-1,j)),2.0);

	is3 = tttw*pow( (b(i+1,j) - 3.0*b(i,j) + 3.0*b(i-1,j) - b(i-2,j)),2.0)
		+ fourth*pow( (3.0*b(i+1,j) - 7.0*b(i,j) + 5.0*b(i-1,j) - b(i-2,j)),2.0);
}

void sflow_iweno_hj::is_east(slice &b)
{
	is1 = tttw*pow( (-b(i,j-3) + 3.0*b(i,j-2) - 3.0*b(i,j-1) + b(i,j)),2.0)
		+ fourth*pow( (-b(i,j-3) + 5.0*b(i,j-2) - 7.0*b(i,j-1) + 3.0*b(i,j)),2.0);

	is2 = tttw*pow( (-b(i,j-2) + 3.0*b(i,j-1) - 3.0*b(i,j) + b(i,j+1)),2.0)
		+ fourth*pow( (-b(i,j-2) + b(i,j-1) +b(i,j) - b(i,j+1)),2.0);

	is3 = tttw*pow( (-b(i,j-1) + 3.0*b(i,j) - 3.0*b(i,j+1) + b(i,j+2)),2.0)
		+ fourth*pow( (-3.0*b(i,j-1) + 7.0*b(i,j) - 5.0*b(i,j+1) + b(i,j+2)),2.0);
}

void sflow_iweno_hj::is_west(slice &b)
{
	is1 = tttw*pow( (b(i,j+3) - 3.0*b(i,j+2) + 3.0*b(i,j+1) - b(i,j)),2.0)
		+ fourth*pow( (b(i,j+3) - 5.0*b(i,j+2) + 7.0*b(i,j+1) - 3.0*b(i,j)),2.0);

	is2 = tttw*pow( (b(i,j+2) - 3.0*b(i,j+1) + 3.0*b(i,j) - b(i,j-1)),2.0)
		+ fourth*pow( (b(i,j+2) - b(i,j+1) - b(i,j) + b(i,j-1)),2.0);

	is3 = tttw*pow( (b(i,j+1) - 3.0*b(i,j) + 3.0*b(i,j-1) - b(i,j-2)),2.0)
		+ fourth*pow( (3.0*b(i,j+1) - 7.0*b(i,j) + 5.0*b(i,j-1) - b(i,j-2)),2.0);
}

void sflow_iweno_hj::alpha_calc()
{
	alpha1=tenth/pow(epsilon+is1,2.0);
	alpha2=sixten/pow(epsilon+is2,2.0);
	alpha3=treten/pow(epsilon+is3,2.0);
}

void sflow_iweno_hj::weight_calc()
{
	w1=alpha1/(alpha1+alpha2+alpha3);
	w2=alpha2/(alpha1+alpha2+alpha3);
	w3=alpha3/(alpha1+alpha2+alpha3);

}
