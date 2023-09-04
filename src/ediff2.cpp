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

#include"ediff2.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

ediff2::ediff2(lexer* p):gradient(p)
{

    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_scalar=80;
}

ediff2::~ediff2()
{
}

void ediff2::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &visc, field &eddyv, double sig, double alpha)
{	
    
    LOOP
	 a->L(i,j,k) += ((b(i+1,j,k)-b(i,j,k))*0.5*(visc(i+1,j,k)+eddyv(i+1,j,k)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DXP[IP])
		-(b(i,j,k)-b(i-1,j,k))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i-1,j,k)+eddyv(i-1,j,k)/sig)*(1.0/p->DXP[IM1]))/p->DXN[IP]

		+((b(i,j+1,k)-b(i,j,k))*0.5*(visc(i,j+1,k)+eddyv(i,j+1,k)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DYP[JP])
		-(b(i,j,k)-b(i,j-1,k))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i,j-1,k)+eddyv(i,j-1,k)/sig)*(1.0/p->DYP[JM1]))/p->DYN[JP]

		+((b(i,j,k+1)-b(i,j,k))*0.5*(visc(i,j,k+1)+eddyv(i,j,k+1)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DZP[KP])
		-(b(i,j,k)-b(i,j,k-1))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i,j,k-1)+eddyv(i,j,k-1)/sig)*(1.0/p->DZP[KM1]))/p->DZN[KP];
}

void ediff2::diff_scalar(lexer* p, fdm* a, ghostcell* pgc, solver* psolv, field &diff, field &b, field &visc, field &eddyv, double sig, double alpha)
{	
    LOOP
    {
        a->L(i,j,k) += ((b(i+1,j,k)-b(i,j,k))*0.5*(visc(i+1,j,k)+eddyv(i+1,j,k)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DXP[IP])
            -(b(i,j,k)-b(i-1,j,k))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i-1,j,k)+eddyv(i-1,j,k)/sig)*(1.0/p->DXP[IM1]))/p->DXN[IP]

            +((b(i,j+1,k)-b(i,j,k))*0.5*(visc(i,j+1,k)+eddyv(i,j+1,k)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DYP[JP])
            -(b(i,j,k)-b(i,j-1,k))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i,j-1,k)+eddyv(i,j-1,k)/sig)*(1.0/p->DYP[JM1]))/p->DYN[JP]

            +((b(i,j,k+1)-b(i,j,k))*0.5*(visc(i,j,k+1)+eddyv(i,j,k+1)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DZP[KP])
            -(b(i,j,k)-b(i,j,k-1))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i,j,k-1)+eddyv(i,j,k-1)/sig)*(1.0/p->DZP[KM1]))/p->DZN[KP];
        
        diff(i,j,k) = b(i,j,k);
    }

}

void ediff2::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	double visc_ddy_p,visc_ddy_m,visc_ddz_p,visc_ddz_m;

    ULOOP
	{
	u_ijk=u(i,j,k);
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=a->visc(i,j,k);

	visc_ddy_p = (visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k))*0.25;
	visc_ddy_m = (a->visc(i,j-1,k)+a->eddyv(i,j-1,k)  +a->visc(i+1,j-1,k)+a->eddyv(i+1,j-1,k) + visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k))*0.25;
	visc_ddz_p = (visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
	visc_ddz_m = (a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + a->visc(i+1,j,k-1)+a->eddyv(i+1,j,k-1) + visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k))*0.25;
	
    a->F(i,j,k) += 2.0*((u(i+1,j,k)-u_ijk)*((a->visc(i+1,j,k)+a->eddyv(i+1,j,k))/p->DXN[IP1])
                        -(u_ijk-u(i-1,j,k))*((visc_ijk+ev_ijk)/p->DXN[IP]))/p->DXP[IP]

        +   ((u(i,j+1,k)-u_ijk)*(visc_ddy_p/p->DYP[JP])

            -(u_ijk-u(i,j-1,k))*(visc_ddy_m/p->DYP[JM1]))/p->DYN[JP]

        +   ((u(i,j,k+1)-u_ijk)*(visc_ddz_p/p->DZP[KP])

            -(u_ijk-u(i,j,k-1))*(visc_ddz_m/p->DZP[KM1]))/p->DZN[KP]


        + ((v(i+1,j,k)-v(i,j,k))*(visc_ddy_p/p->DXP[IP]) - (v(i+1,j-1,k)-v(i,j-1,k))*(visc_ddy_m/p->DXP[IP]))/p->DYN[JP]

        + ((w(i+1,j,k)-w(i,j,k))*(visc_ddz_p/p->DXP[IP]) - (w(i+1,j,k-1)-w(i,j,k-1))*(visc_ddz_m/p->DXP[IP]))/p->DZN[KP];
		
	}
}

void ediff2::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	double visc_ddx_p,visc_ddx_m,visc_ddz_p,visc_ddz_m;
	
    VLOOP
	{
	v_ijk=v(i,j,k);
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=a->visc(i,j,k);
	visc_ddx_p = (visc_ijk+ev_ijk + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k))*0.25;
	visc_ddx_m = (a->visc(i-1,j,k)+a->eddyv(i-1,j,k) + a->visc(i-1,j+1,k)+a->eddyv(i-1,j+1,k) + visc_ijk+ev_ijk + a->visc(i,j+1,k)+a->eddyv(i,j+1,k))*0.25;
	visc_ddz_p = (visc_ijk+ev_ijk + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1))*0.25;
	visc_ddz_m = (a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + a->visc(i,j+1,k-1)+a->eddyv(i,j+1,k-1) + visc_ijk+ev_ijk + a->visc(i,j+1,k)+a->eddyv(i,j+1,k))*0.25;
	
	a->G(i,j,k) +=	   ((v(i+1,j,k)-v_ijk)*(visc_ddx_p/p->DXP[IP])

				- (v_ijk-v(i-1,j,k))*(visc_ddx_m/p->DXP[IM1]))/p->DXN[IP]

		+ 2.0*((v(i,j+1,k)-v_ijk)*((a->visc(i,j+1,k)+a->eddyv(i,j+1,k))/p->DYN[JP1])
		     -(v_ijk-v(i,j-1,k))*((visc_ijk+ev_ijk)/p->DYN[JP]))/p->DYP[JP]

		+ ((v(i,j,k+1)-v_ijk)*(visc_ddz_p/p->DZP[KP])

		  -(v_ijk-v(i,j,k-1))*(visc_ddz_m/p->DZP[KM1]))/p->DZN[KP]
          

		+ ((u(i,j+1,k)-u(i,j,k))*(visc_ddx_p/p->DYP[JP]) - (u(i-1,j+1,k)-u(i-1,j,k))*(visc_ddx_m/p->DYP[JP]))/p->DXN[IP]

		+ ((w(i,j+1,k)-w(i,j,k))*(visc_ddz_p/p->DYP[JP]) - (w(i,j+1,k-1)-w(i,j,k-1))*(visc_ddz_m/p->DYP[JP]))/p->DZN[KP];
	}
}

void ediff2::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &u, field &v, field &w, double alpha)
{
	double visc_ddx_p,visc_ddx_m,visc_ddy_p,visc_ddy_m;
	

    WLOOP
	{
	w_ijk=w(i,j,k);
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=a->visc(i,j,k);
	visc_ddx_p = (visc_ijk+ev_ijk + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
	visc_ddx_m = (a->visc(i-1,j,k)+a->eddyv(i-1,j,k) + a->visc(i-1,j,k+1)+a->eddyv(i-1,j,k+1) + visc_ijk+ev_ijk + a->visc(i,j,k+1)+a->eddyv(i,j,k+1))*0.25;
	visc_ddy_p = (visc_ijk+ev_ijk + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1))*0.25;
	visc_ddy_m = (a->visc(i,j-1,k)+a->eddyv(i,j-1,k) + a->visc(i,j-1,k+1)+a->eddyv(i,j-1,k+1) + visc_ijk+ev_ijk + a->visc(i,j,k+1)+a->eddyv(i,j,k+1))*0.25;
	
    a->H(i,j,k) +=  ((w(i+1,j,k)-w_ijk)*(visc_ddx_p/p->DXP[IP])
            
					-(w_ijk-w(i-1,j,k))*(visc_ddx_m/p->DXP[IM1]))/p->DXN[IP]

				+   ((w(i,j+1,k)-w_ijk)*(visc_ddy_p/p->DYP[JP])

					-(w_ijk-w(i,j-1,k))*(visc_ddy_m/p->DYP[JM1]))/p->DYN[JP]

        + 2.0*((w(i,j,k+1)-w_ijk)*((a->visc(i,j,k+1)+a->eddyv(i,j,k+1))/p->DZN[KP1])
        -(w_ijk-w(i,j,k-1))*((visc_ijk+ev_ijk)/p->DZN[KP]))/p->DZP[KP]
        

        + ((u(i,j,k+1)-u(i,j,k))*(visc_ddx_p/p->DZP[KP]) - (u(i-1,j,k+1)-u(i-1,j,k))*(visc_ddx_m/p->DZP[KP]))/p->DXN[IP]

     	+ ((v(i,j,k+1)-v(i,j,k))*(visc_ddy_p/p->DZP[KP]) - (v(i,j-1,k+1)-v(i,j-1,k))*(visc_ddy_m/p->DZP[KP]))/p->DYN[JP];
	}
}

void ediff2::idiff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &eddyv, double sig, double alpha)
{	
}

void ediff2::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &u_in, field &u, field &v, field &w, double alpha)
{
	double visc_ddy_p,visc_ddy_m,visc_ddz_p,visc_ddz_m;
    
    ULOOP
    diff(i,j,k) = u_in(i,j,k);
    
    ULOOP
	{
	u_ijk=u(i,j,k);
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=a->visc(i,j,k);

	visc_ddy_p = (visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k))*0.25;
	visc_ddy_m = (a->visc(i,j-1,k)+a->eddyv(i,j-1,k)  +a->visc(i+1,j-1,k)+a->eddyv(i+1,j-1,k) + visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k))*0.25;
	visc_ddz_p = (visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
	visc_ddz_m = (a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + a->visc(i+1,j,k-1)+a->eddyv(i+1,j,k-1) + visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k))*0.25;
	
    a->F(i,j,k) += 2.0*((u(i+1,j,k)-u_ijk)*((a->visc(i+1,j,k)+a->eddyv(i+1,j,k))/p->DXN[IP1])
                        -(u_ijk-u(i-1,j,k))*((visc_ijk+ev_ijk)/p->DXN[IP]))/p->DXP[IP]

        +   ((u(i,j+1,k)-u_ijk)*(visc_ddy_p/p->DYP[JP])

            -(u_ijk-u(i,j-1,k))*(visc_ddy_m/p->DYP[JM1]))/p->DYN[JP]

        +   ((u(i,j,k+1)-u_ijk)*(visc_ddz_p/p->DZP[KP])

            -(u_ijk-u(i,j,k-1))*(visc_ddz_m/p->DZP[KM1]))/p->DZN[KP]


        + ((v(i+1,j,k)-v(i,j,k))*(visc_ddy_p/p->DXP[IP]) - (v(i+1,j-1,k)-v(i,j-1,k))*(visc_ddy_m/p->DXP[IP]))/p->DYN[JP]

        + ((w(i+1,j,k)-w(i,j,k))*(visc_ddz_p/p->DXP[IP]) - (w(i+1,j,k-1)-w(i,j,k-1))*(visc_ddz_m/p->DXP[IP]))/p->DZN[KP];
	}
    
    pgc->start1(p,diff,gcval_u);
}

void ediff2::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &v_in, field &u, field &v, field &w, double alpha)
{
	double visc_ddx_p,visc_ddx_m,visc_ddz_p,visc_ddz_m;
	
    VLOOP
    diff(i,j,k) = v_in(i,j,k);
    
    VLOOP
	{
	v_ijk=v(i,j,k);
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=a->visc(i,j,k);
	visc_ddx_p = (visc_ijk+ev_ijk + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i+1,j+1,k)+a->eddyv(i+1,j+1,k))*0.25;
	visc_ddx_m = (a->visc(i-1,j,k)+a->eddyv(i-1,j,k) + a->visc(i-1,j+1,k)+a->eddyv(i-1,j+1,k) + visc_ijk+ev_ijk + a->visc(i,j+1,k)+a->eddyv(i,j+1,k))*0.25;
	visc_ddz_p = (visc_ijk+ev_ijk + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1))*0.25;
	visc_ddz_m = (a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + a->visc(i,j+1,k-1)+a->eddyv(i,j+1,k-1) + visc_ijk+ev_ijk + a->visc(i,j+1,k)+a->eddyv(i,j+1,k))*0.25;
	
	a->G(i,j,k) +=	   ((v(i+1,j,k)-v_ijk)*(visc_ddx_p/p->DXP[IP])

				- (v_ijk-v(i-1,j,k))*(visc_ddx_m/p->DXP[IM1]))/p->DXN[IP]

		+ 2.0*((v(i,j+1,k)-v_ijk)*((a->visc(i,j+1,k)+a->eddyv(i,j+1,k))/p->DYN[JP1])
		     -(v_ijk-v(i,j-1,k))*((visc_ijk+ev_ijk)/p->DYN[JP]))/p->DYP[JP]

		+ ((v(i,j,k+1)-v_ijk)*(visc_ddz_p/p->DZP[KP])

		  -(v_ijk-v(i,j,k-1))*(visc_ddz_m/p->DZP[KM1]))/p->DZN[KP]
          

		+ ((u(i,j+1,k)-u(i,j,k))*(visc_ddx_p/p->DYP[JP]) - (u(i-1,j+1,k)-u(i-1,j,k))*(visc_ddx_m/p->DYP[JP]))/p->DXN[IP]

		+ ((w(i,j+1,k)-w(i,j,k))*(visc_ddz_p/p->DYP[JP]) - (w(i,j+1,k-1)-w(i,j,k-1))*(visc_ddz_m/p->DYP[JP]))/p->DZN[KP];
	}
    
    pgc->start2(p,diff,gcval_v);
}

void ediff2::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &w_in, field &u, field &v, field &w, double alpha)
{
	double visc_ddx_p,visc_ddx_m,visc_ddy_p,visc_ddy_m;
	
    WLOOP
    diff(i,j,k) = w_in(i,j,k);
    
    WLOOP
	{
	w_ijk=w(i,j,k);
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=a->visc(i,j,k);
	visc_ddx_p = (visc_ijk+ev_ijk + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
	visc_ddx_m = (a->visc(i-1,j,k)+a->eddyv(i-1,j,k) + a->visc(i-1,j,k+1)+a->eddyv(i-1,j,k+1) + visc_ijk+ev_ijk + a->visc(i,j,k+1)+a->eddyv(i,j,k+1))*0.25;
	visc_ddy_p = (visc_ijk+ev_ijk + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i,j+1,k)+a->eddyv(i,j+1,k) + a->visc(i,j+1,k+1)+a->eddyv(i,j+1,k+1))*0.25;
	visc_ddy_m = (a->visc(i,j-1,k)+a->eddyv(i,j-1,k) + a->visc(i,j-1,k+1)+a->eddyv(i,j-1,k+1) + visc_ijk+ev_ijk + a->visc(i,j,k+1)+a->eddyv(i,j,k+1))*0.25;
	
    a->H(i,j,k) +=  ((w(i+1,j,k)-w_ijk)*(visc_ddx_p/p->DXP[IP])
            
					-(w_ijk-w(i-1,j,k))*(visc_ddx_m/p->DXP[IM1]))/p->DXN[IP]

				+   ((w(i,j+1,k)-w_ijk)*(visc_ddy_p/p->DYP[JP])

					-(w_ijk-w(i,j-1,k))*(visc_ddy_m/p->DYP[JM1]))/p->DYN[JP]

        + 2.0*((w(i,j,k+1)-w_ijk)*((a->visc(i,j,k+1)+a->eddyv(i,j,k+1))/p->DZN[KP1])
        -(w_ijk-w(i,j,k-1))*((visc_ijk+ev_ijk)/p->DZN[KP]))/p->DZP[KP]
        

        + ((u(i,j,k+1)-u(i,j,k))*(visc_ddx_p/p->DZP[KP]) - (u(i-1,j,k+1)-u(i-1,j,k))*(visc_ddx_m/p->DZP[KP]))/p->DXN[IP]

     	+ ((v(i,j,k+1)-v(i,j,k))*(visc_ddy_p/p->DZP[KP]) - (v(i,j-1,k+1)-v(i,j-1,k))*(visc_ddy_m/p->DZP[KP]))/p->DYN[JP];
	
	}
    
    pgc->start3(p,diff,gcval_w);
}

