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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ediff2_2D.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"

ediff2_2D::ediff2_2D(lexer* p):gradient(p)
{

    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_scalar=80;
}

ediff2_2D::~ediff2_2D()
{
}

void ediff2_2D::diff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &visc, field &eddyv, double sig, double alpha)
{	
    
    LOOP
	 a->L(i,j,k) += ((b(i+1,j,k)-b(i,j,k))*0.5*(visc(i+1,j,k)+eddyv(i+1,j,k)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DXP[IP])
		-(b(i,j,k)-b(i-1,j,k))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i-1,j,k)+eddyv(i-1,j,k)/sig)*(1.0/p->DXP[IM1]))/p->DXN[IP]

		+((b(i,j,k+1)-b(i,j,k))*0.5*(visc(i,j,k+1)+eddyv(i,j,k+1)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DZP[KP])
		-(b(i,j,k)-b(i,j,k-1))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i,j,k-1)+eddyv(i,j,k-1)/sig)*(1.0/p->DZP[KM1]))/p->DZN[KP];
}

void ediff2_2D::diff_scalar(lexer* p, fdm* a, ghostcell* pgc, solver* psolv, field &diff, field &b, field &visc, field &eddyv, double sig, double alpha)
{	
    LOOP
    {
        a->L(i,j,k) += ((b(i+1,j,k)-b(i,j,k))*0.5*(visc(i+1,j,k)+eddyv(i+1,j,k)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DXP[IP])
            -(b(i,j,k)-b(i-1,j,k))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i-1,j,k)+eddyv(i-1,j,k)/sig)*(1.0/p->DXP[IM1]))/p->DXN[IP]

            +((b(i,j,k+1)-b(i,j,k))*0.5*(visc(i,j,k+1)+eddyv(i,j,k+1)/sig+visc(i,j,k)+eddyv(i,j,k)/sig)*(1.0/p->DZP[KP])
            -(b(i,j,k)-b(i,j,k-1))*0.5*(visc(i,j,k)+eddyv(i,j,k)/sig+visc(i,j,k-1)+eddyv(i,j,k-1)/sig)*(1.0/p->DZP[KM1]))/p->DZN[KP];
        
        diff(i,j,k) = b(i,j,k);
    }

}

void ediff2_2D::idiff_scalar(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &b, field &eddyv, double sig, double alpha)
{	
}

void ediff2_2D::diff_u(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &u_in, field &u, field &v, field &w, double alpha)
{
	double visc_ddy_p,visc_ddy_m,visc_ddz_p,visc_ddz_m;
    
    ULOOP
    diff(i,j,k) = u_in(i,j,k);
    
    ULOOP
	{
	u_ijk=u(i,j,k);
	ev_ijk=a->eddyv(i,j,k);
	visc_ijk=a->visc(i,j,k);

	visc_ddz_p = (visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k) + a->visc(i,j,k+1)+a->eddyv(i,j,k+1) + a->visc(i+1,j,k+1)+a->eddyv(i+1,j,k+1))*0.25;
	visc_ddz_m = (a->visc(i,j,k-1)+a->eddyv(i,j,k-1) + a->visc(i+1,j,k-1)+a->eddyv(i+1,j,k-1) + visc_ijk+ev_ijk + a->visc(i+1,j,k)+a->eddyv(i+1,j,k))*0.25;
	
    a->F(i,j,k) += 2.0*((u(i+1,j,k)-u_ijk)*((a->visc(i+1,j,k)+a->eddyv(i+1,j,k))/p->DXN[IP1])
                        -(u_ijk-u(i-1,j,k))*((visc_ijk+ev_ijk)/p->DXN[IP]))/p->DXP[IP]


        +   ((u(i,j,k+1)-u_ijk)*(visc_ddz_p/p->DZP[KP])

            -(u_ijk-u(i,j,k-1))*(visc_ddz_m/p->DZP[KM1]))/p->DZN[KP]


        + ((w(i+1,j,k)-w(i,j,k))*(visc_ddz_p/p->DXP[IP]) - (w(i+1,j,k-1)-w(i,j,k-1))*(visc_ddz_m/p->DXP[IP]))/p->DZN[KP];
	}
    
    pgc->start1(p,diff,gcval_u);
}

void ediff2_2D::diff_v(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &v_in, field &u, field &v, field &w, double alpha)
{
}

void ediff2_2D::diff_w(lexer* p, fdm* a, ghostcell *pgc, solver *psolv, field &diff, field &w_in, field &u, field &v, field &w, double alpha)
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
	
    a->H(i,j,k) +=  ((w(i+1,j,k)-w_ijk)*(visc_ddx_p/p->DXP[IP])
            
					-(w_ijk-w(i-1,j,k))*(visc_ddx_m/p->DXP[IM1]))/p->DXN[IP]


        + 2.0*((w(i,j,k+1)-w_ijk)*((a->visc(i,j,k+1)+a->eddyv(i,j,k+1))/p->DZN[KP1])
        -(w_ijk-w(i,j,k-1))*((visc_ijk+ev_ijk)/p->DZN[KP]))/p->DZP[KP]
        

        + ((u(i,j,k+1)-u(i,j,k))*(visc_ddx_p/p->DZP[KP]) - (u(i-1,j,k+1)-u(i-1,j,k))*(visc_ddx_m/p->DZP[KP]))/p->DXN[IP];
	
	}
    
    pgc->start3(p,diff,gcval_w);
}

