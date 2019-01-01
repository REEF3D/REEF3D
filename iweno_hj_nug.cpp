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
/*
#include"iweno_hj_nug.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_HJ_CDS2.h"
#include"flux_HJ_CDS2_vrans.h"

iweno_hj_nug::iweno_hj_nug(lexer *p):weno_nug_func(p)
{
    if(p->B269==0)
    pflux = new flux_HJ_CDS2(p);
    
    if(p->B269>=1)
    pflux = new flux_HJ_CDS2_vrans(p);
}

iweno_hj_nug::~iweno_hj_nug()
{
}

void iweno_hj_nug::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
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

void iweno_hj_nug::wenoloop1(lexer *p, fdm *a, field& f, int ipol, field& uvel, field& vvel, field& wvel)
{
	count=0;
    
    

	ULOOP
	{
        pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
    

			if(iadvec>=0.0)
			{
			is_south(f);
			alpha_calc();
			weight_calc();
			aij_south(p,a,f,a->F);
			}

			if(iadvec<0.0)
			{
			is_north(f);
			alpha_calc();
			weight_calc();
			aij_north(p,a,f,a->F);
			}


			
			if(jadvec>=0.0)
			{
			is_east(f);
			alpha_calc();
			weight_calc();
			aij_east(p,a,f,a->F);
			}

			if(jadvec<0.0)
			{
			is_west(f);
			alpha_calc();
			weight_calc();
			aij_west(p,a,f,a->F);
			}



			if(kadvec>=0.0)
			{
			is_bottom(f);
			alpha_calc();
			weight_calc();
			aij_bottom(p,a,f,a->F);
			}

			if(kadvec<0.0)
			{
			is_top(f);
			alpha_calc();
			weight_calc();
			aij_top(p,a,f,a->F);
			}
		 ++count;
	}
}

void iweno_hj_nug::wenoloop2(lexer *p, fdm *a, field& f, int ipol, field& uvel, field& vvel, field& wvel)
{
	count=0;

	VLOOP
	{
        pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);


			if(iadvec>=0.0)
			{
			is_south(f);
			alpha_calc();
			weight_calc();
			aij_south(p,a,f,a->G);
			}

			if(iadvec<0.0)
			{
			is_north(f);
			alpha_calc();
			weight_calc();
			aij_north(p,a,f,a->G);
			}



			if(jadvec>=0.0)
			{
			is_east(f);
			alpha_calc();
			weight_calc();
			aij_east(p,a,f,a->G);
			}

			if(jadvec<0.0)
			{
			is_west(f);
			alpha_calc();
			weight_calc();
			aij_west(p,a,f,a->G);
			}



			if(kadvec>=0.0)
			{
			is_bottom(f);
			alpha_calc();
			weight_calc();
			aij_bottom(p,a,f,a->G);
			}

			if(kadvec<0.0)
			{
			is_top(f);
			alpha_calc();
			weight_calc();
			aij_top(p,a,f,a->G);
			}

		 ++count;
	}
}

void iweno_hj_nug::wenoloop3(lexer *p, fdm *a, field& f, int ipol, field& uvel, field& vvel, field& wvel)
{
	count=0;

	WLOOP
	{
	
		pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
                
			if(iadvec>=0.0)
			{
			iqmin(p,a,b,uvel,ipol);
            is_min_x();
            weight_min_x();
			aij_south(p,a,f,a->H);
			}

			if(iadvec<0.0)
			{
			iqmax(p,a,b,uvel,ipol);
            is_max_x();
            weight_max_x();
			aij_north(p,a,f,a->H);
			}


		
			if(jadvec>=0.0)
			{
			jqmin(p,a,b,vvel,ipol);
            is_min_y();
            weight_min_y();
			aij_east(p,a,f,a->H);
			}

			if(jadvec<0.0)
			{
			jqmax(p,a,b,vvel,ipol);
            is_max_y();
            weight_max_y();
			aij_west(p,a,f,a->H);
			}


			
			if(kadvec>=0.0)
			{
			kqmin(p,a,b,wvel,ipol);
            is_min_z();
            weight_min_z();
			aij_bottom(p,a,f,a->H);
			}

			if(kadvec<0.0)
			{
			kqmax(p,a,b,wvel,ipol);
            is_max_z();
            weight_max_z();
			aij_top(p,a,f,a->H);
			}
			
		 ++count;
	}
}

void iweno_hj_nug::wenoloop4(lexer *p, fdm *a, field& f, int ipol, field& uvel, field& vvel, field& wvel)
{
	count=0;

	LOOP
	{
		
        pflux->u_flux(a,ipol,uvel,iadvec,ivel2);
        pflux->v_flux(a,ipol,vvel,jadvec,jvel2);
        pflux->w_flux(a,ipol,wvel,kadvec,kvel2);
			
			if(iadvec>=0.0)
			{
			is_south(f);
			alpha_calc();
			weight_calc();
			aij_south(p,a,f,a->L);
			}

			if(iadvec<0.0)
			{
			is_north(f);
			alpha_calc();
			weight_calc();
			aij_north(p,a,f,a->L);
			}


			if(jadvec>=0.0)
			{
			is_east(f);
			alpha_calc();
			weight_calc();
			aij_east(p,a,f,a->L);
			}

			if(jadvec<0.0)
			{
			is_west(f);
			alpha_calc();
			weight_calc();
			aij_west(p,a,f,a->L);
			}

	
			if(kadvec>=0.0)
			{
			is_bottom(f);
			alpha_calc();
			weight_calc();
			aij_bottom(p,a,f,a->L);
			}

			if(kadvec<0.0)
			{
			is_top(f);
			alpha_calc();
			weight_calc();
			aij_top(p,a,f,a->L);
			}
     ++count;
	}
}

void iweno_hj_nug::aij_south(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k)    += (w3x*qfx[IP][uf][2][0])*iadvec*DX[IM3]*f(i-3,j,k)
    
                qfx[IP][uf][2][0]*(q1
				 - (w2*sixth + w1*1.5)*iadvec*deltin*f(i-2,j,k)
                 
				 + (w3*sixth)*iadvec*deltin*f(i+2,j,k);
				 
	a->M.p[count] = (-w3*0.5 + w2*0.5 + w1*elvsix)*iadvec*deltin;
					 
	a->M.s[count] = (-w3*third -w2 - w1*3.0)*iadvec*deltin;
	a->M.n[count] = (w3 + w2*third)*iadvec*deltin;
    

}

void iweno_hj_nug::aij_north(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k)   += -(w3*sixth)*iadvec*deltin*f(i-2,j,k)
				-	(-w1*1.5 - w2*sixth)*iadvec*deltin*f(i+2,j,k)
				-   (w1*third)*iadvec*deltin*f(i+3,j,k);
					
	a->M.p[count] = (-w1*elvsix - w2*0.5 + w3*0.5)*iadvec*deltin;
					 
	a->M.s[count] = (-w2*third - w3)*iadvec*deltin;
	a->M.n[count] = (w1*3.0 + w2 + w3*third)*iadvec*deltin;
}

void iweno_hj_nug::aij_east(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k) 	+= (w1*third)*jadvec*deltin*f(i,j-3,k)
				-  (w2*sixth + w1*1.5)*jadvec*deltin*f(i,j-2,k)
				-  (-w3*sixth)*jadvec*deltin*f(i,j+2,k);
					
	a->M.p[count] += (-w3*0.5 +w2*0.5 + w1*elvsix)*jadvec*deltin;
					 
	a->M.e[count] = (-w3*third -w2 - w1*3.0)*jadvec*deltin;
	a->M.w[count] = (w3 + w2*third)*jadvec*deltin;
}

void iweno_hj_nug::aij_west(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k) 	+= -(w3*sixth)*jadvec*deltin*f(i,j-2,k)
				-	(-w1*1.5 - w2*sixth)*jadvec*deltin*f(i,j+2,k)
				-	 (w1*third)*jadvec*deltin*f(i,j+3,k);
				
	a->M.p[count] += (-w1*elvsix - w2*0.5 + w3*0.5)*jadvec*deltin;
					 
	a->M.e[count] = (-w2*third - w3)*jadvec*deltin;
	a->M.w[count] = (w1*3.0 + w2 + w3*third)*jadvec*deltin;
}

void iweno_hj_nug::aij_bottom(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k) 	+= (w1*third)*kadvec*deltin*f(i,j,k-3)
				-  (w2*sixth + w1*1.5)*kadvec*deltin*f(i,j,k-2)
				-  (-w3*sixth)*kadvec*deltin*f(i,j,k+2);
	
	a->M.p[count] += (-w3*0.5 +w2*0.5 + w1*elvsix)*kadvec*deltin;
					 
	a->M.b[count] = (-w3*third -w2 - w1*3.0)*kadvec*deltin;
	a->M.t[count] = (w3 + w2*third)*kadvec*deltin;
}

void iweno_hj_nug::aij_top(lexer* p, fdm* a, field &f, field &F)
{
	F(i,j,k) 	+= -(w3*sixth)*kadvec*deltin*f(i,j,k-2)
				-   (-w1*1.5 - w2*sixth)*kadvec*deltin*f(i,j,k+2)
				-   (w1*third)*kadvec*deltin*f(i,j,k+3);
				
	a->M.p[count] += (-w1*elvsix - w2*0.5 + w3*0.5)*kadvec*deltin;
					 
	a->M.b[count] = (-w2*third - w3)*kadvec*deltin;
	a->M.t[count] = (w1*3.0 + w2 + w3*third)*kadvec*deltin;
}


void weno_flux_nug::iqmin(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{	
	q1 = f(i-2,j,k);
	q2 = f(i-1,j,k);
	q3 = f(i,j,k);
	q4 = f(i+1,j,k);
	q5 = f(i+2,j,k);
}

void weno_flux_nug::jqmin(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = f(i,j-2,k);
	q2 = f(i,j-1,k);
	q3 = f(i,j,k);
	q4 = f(i,j+1,k);
	q5 = f(i,j+2,k);
}

void weno_flux_nug::kqmin(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k-2);
	q2 = f(i,j,k-1);
	q3 = f(i,j,k);
	q4 = f(i,j,k+1);
	q5 = f(i,j,k+2);
}

void weno_flux_nug::iqmax(lexer *p,fdm *a, field& f, field& uvel, int ipol)
{
    q1 = f(i-1,j,k);
	q2 = f(i,j,k);
	q3 = f(i+1,j,k);
	q4 = f(i+2,j,k);
	q5 = f(i+3,j,k);
}

void weno_flux_nug::jqmax(lexer *p,fdm *a, field& f, field& vvel, int ipol)
{
	q1 = f(i,j-1,k);
	q2 = f(i,j,k);
	q3 = f(i,j+1,k);
	q4 = f(i,j+2,k);
	q5 = f(i,j+3,k);
}

void weno_flux_nug::kqmax(lexer *p,fdm *a, field& f, field& wvel, int ipol)
{
	q1 = f(i,j,k-1);
	q2 = f(i,j,k);
	q3 = f(i,j,k+1);
	q4 = f(i,j,k+2);
	q5 = f(i,j,k+3);
}*/