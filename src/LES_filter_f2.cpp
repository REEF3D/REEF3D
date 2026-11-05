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

#include"LES_filter_f2.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"

// T21 2: second-order high-pass filter, u' = (I-G)^2 u (Stolz, Schlatter & Kleiser 2005), with the same
// discrete filter G as T21 1 (LES_filter_f1: (1,2,1)/4 in x, y and z, u' = (I-G) u). The Smagorinsky or
// WALE eddy viscosity is evaluated from the strain of u', so it acts only on the smallest resolved scales.
LES_filter_f2::LES_filter_f2(lexer* p, fdm* a) : strain(p), u1(p), ut1(p), ut2(p), ut3(p),
                                                 v1(p), vt1(p), vt2(p), vt3(p),
                                                 w1(p), wt1(p), wt2(p), wt3(p)
{
}

LES_filter_f2::~LES_filter_f2()
{
}

void LES_filter_f2::start(lexer *p, fdm *a, ghostcell *pgc, field &uprime, field &vprime, field &wprime, int gcval)
{
	if(gcval==10)
	{
        // u1 = (I-G) u
        filter_u(p,pgc,a->u,ut3,gcval);
        
        ULOOP
        u1(i,j,k) = a->u(i,j,k) - ut3(i,j,k);
        
        pgc->start1(p,u1,gcval);
        
        // u' = (I-G) u1
        filter_u(p,pgc,u1,ut3,gcval);
        
        ULOOP
        uprime(i,j,k) = u1(i,j,k) - ut3(i,j,k);
        
        pgc->start1(p,uprime,gcval);
	}

	if(gcval==11)
	{
        filter_v(p,pgc,a->v,vt3,gcval);
        
        VLOOP
        v1(i,j,k) = a->v(i,j,k) - vt3(i,j,k);
        
        pgc->start2(p,v1,gcval);
        
        filter_v(p,pgc,v1,vt3,gcval);
        
        VLOOP
        vprime(i,j,k) = v1(i,j,k) - vt3(i,j,k);
        
        pgc->start2(p,vprime,gcval);
	}

	if(gcval==12)
	{
        filter_w(p,pgc,a->w,wt3,gcval);
        
        WLOOP
        w1(i,j,k) = a->w(i,j,k) - wt3(i,j,k);
        
        pgc->start3(p,w1,gcval);
        
        filter_w(p,pgc,w1,wt3,gcval);
        
        WLOOP
        wprime(i,j,k) = w1(i,j,k) - wt3(i,j,k);
        
        pgc->start3(p,wprime,gcval);
	}
}

// the stencils are those of LES_filter_f1 (x, then y, then z), applied to the field f
void LES_filter_f2::filter_u(lexer *p, ghostcell *pgc, field &f, field &out, int gcval)
{
    ULOOP
    ut1(i,j,k) = f(i,j,k) + 0.5*(p->DXN[IP]*(f(i+1,j,k) - f(i,j,k)) - p->DXN[IP1]*(f(i,j,k) - f(i-1,j,k)))/(p->DXN[IP]+p->DXN[IP1]);
    
    pgc->start1(p,ut1,gcval);
    
    ULOOP
    ut2(i,j,k) = ut1(i,j,k) + 0.5*(p->DYP[JM1]*(ut1(i,j+1,k) - ut1(i,j,k)) - p->DYP[JP]*(ut1(i,j,k) - ut1(i,j-1,k)))/(p->DYP[JP]+p->DYP[JM1]);
    
    pgc->start1(p,ut2,gcval);
    
    ULOOP
    out(i,j,k) = ut2(i,j,k) + 0.5*(p->DZP[KM1]*(ut2(i,j,k+1) - ut2(i,j,k)) - p->DZP[KP]*(ut2(i,j,k) - ut2(i,j,k-1)))/(p->DZP[KP]+p->DZP[KM1]);
    
    pgc->start1(p,out,gcval);
}

void LES_filter_f2::filter_v(lexer *p, ghostcell *pgc, field &f, field &out, int gcval)
{
    VLOOP
    vt1(i,j,k) = f(i,j,k) + 0.5*(p->DXP[IM1]*(f(i+1,j,k) - f(i,j,k)) - p->DXP[IP]*(f(i,j,k) - f(i-1,j,k)))/(p->DXP[IP]+p->DXP[IM1]);
    
    pgc->start2(p,vt1,gcval);
    
    VLOOP
    vt2(i,j,k) = vt1(i,j,k) + 0.5*(p->DYN[JP]*(vt1(i,j+1,k) - vt1(i,j,k)) - p->DYN[JP1]*(vt1(i,j,k) - vt1(i,j-1,k)))/(p->DYN[JP]+p->DYN[JP1]);
    
    pgc->start2(p,vt2,gcval);
    
    VLOOP
    out(i,j,k) = vt2(i,j,k) + 0.5*(p->DZP[KM1]*(vt2(i,j,k+1) - vt2(i,j,k)) - p->DZP[KP]*(vt2(i,j,k) - vt2(i,j,k-1)))/(p->DZP[KP]+p->DZP[KM1]);
    
    pgc->start2(p,out,gcval);
}

void LES_filter_f2::filter_w(lexer *p, ghostcell *pgc, field &f, field &out, int gcval)
{
    WLOOP
    wt1(i,j,k) = f(i,j,k) + 0.5*(p->DXP[IM1]*(f(i+1,j,k) - f(i,j,k)) - p->DXP[IP]*(f(i,j,k) - f(i-1,j,k)))/(p->DXP[IP]+p->DXP[IM1]);
    
    pgc->start3(p,wt1,gcval);
    
    WLOOP
    wt2(i,j,k) = wt1(i,j,k) + 0.5*(p->DYP[JM1]*(wt1(i,j+1,k) - wt1(i,j,k)) - p->DYP[JP]*(wt1(i,j,k) - wt1(i,j-1,k)))/(p->DYP[JP]+p->DYP[JM1]);
    
    pgc->start3(p,wt2,gcval);
    
    WLOOP
    out(i,j,k) = wt2(i,j,k) + 0.5*(p->DZN[KP]*(wt2(i,j,k+1) - wt2(i,j,k)) - p->DZN[KP1]*(wt2(i,j,k) - wt2(i,j,k-1)))/(p->DZN[KP]+p->DZN[KP1]);
    
    pgc->start3(p,out,gcval);
}
