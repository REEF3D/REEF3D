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

#include"vrans_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field.h"

// Darcy-Forchheimer resistance on the Darcy velocity u (van Gent 1995, Jensen et al. 2014):
//   f_i = (A + B |U|) u_i,  A = porA nu,  B = porB  (cell values, see set_cell)
// A and B are averaged from the cells to the face (not alpha, beta, d50 and n separately, which
// underestimates the resistance at the faces next to the open flow), |U| is the magnitude of
// the velocity vector at the face.
// B 268 0: explicit source with u(n)
// B 268 1: point-implicit in the stages of momentum_rk (implicit_drag), no explicit source
double vrans_f::drag_rate(lexer *p, fdm *a, int c, field &U, field &V, field &W)
{
    double A,B,uu,vv,ww;
    
    if(c==0)
    {
    A  = 0.5*(a->porA(i,j,k)*a->visc(i,j,k) + a->porA(i+1,j,k)*a->visc(i+1,j,k));
    B  = 0.5*(a->porB(i,j,k) + a->porB(i+1,j,k));
    uu = U(i,j,k);
    vv = p->j_dir==1 ? 0.25*(V(i,j,k) + V(i+1,j,k) + V(i,j-1,k) + V(i+1,j-1,k)) : 0.0;
    ww = 0.25*(W(i,j,k) + W(i+1,j,k) + W(i,j,k-1) + W(i+1,j,k-1));
    }
    
    else if(c==1)
    {
    A  = 0.5*(a->porA(i,j,k)*a->visc(i,j,k) + a->porA(i,j+1,k)*a->visc(i,j+1,k));
    B  = 0.5*(a->porB(i,j,k) + a->porB(i,j+1,k));
    uu = 0.25*(U(i,j,k) + U(i,j+1,k) + U(i-1,j,k) + U(i-1,j+1,k));
    vv = V(i,j,k);
    ww = 0.25*(W(i,j,k) + W(i,j+1,k) + W(i,j,k-1) + W(i,j+1,k-1));
    }
    
    else
    {
    A  = 0.5*(a->porA(i,j,k)*a->visc(i,j,k) + a->porA(i,j,k+1)*a->visc(i,j,k+1));
    B  = 0.5*(a->porB(i,j,k) + a->porB(i,j,k+1));
    uu = 0.25*(U(i,j,k) + U(i,j,k+1) + U(i-1,j,k) + U(i-1,j,k+1));
    vv = p->j_dir==1 ? 0.25*(V(i,j,k) + V(i,j,k+1) + V(i,j-1,k) + V(i,j-1,k+1)) : 0.0;
    ww = W(i,j,k);
    }
    
    return A + B*sqrt(uu*uu + vv*vv + ww*ww);
}

void vrans_f::u_source(lexer *p, fdm *a)
{
    if(p->B268==1)
    return;
    
    count=0;
    ULOOP
	{
        if(PORVAL1<1.0)
        a->rhsvec.V[count] -= drag_rate(p,a,0,a->u,a->v,a->w)*a->u(i,j,k);
        
	++count;
	}
}

void vrans_f::v_source(lexer *p, fdm *a)
{
    if(p->B268==1)
    return;
    
    count=0;
    VLOOP
	{
        if(PORVAL2<1.0)
        a->rhsvec.V[count] -= drag_rate(p,a,1,a->u,a->v,a->w)*a->v(i,j,k);
        
	++count;
	}
}

void vrans_f::w_source(lexer *p, fdm *a)
{
    if(p->B268==1)
    return;
    
    count=0;
    WLOOP
	{
        if(PORVAL3<1.0)
        a->rhsvec.V[count] -= drag_rate(p,a,2,a->u,a->v,a->w)*a->w(i,j,k);
        
	++count;
	}
}

// point-implicit resistance in a Runge-Kutta stage of momentum_rk: with the explicit terms of the
// stage combined in f (weight w on dt F), the resistance -n (A + B|U|) u(s+1) of the stage gives
//   f <- f / (1 + w dt CPOR n (A + B|U|)),
// K = A + B|U| from the stage velocity (U,V,W). Unconditionally stable, and the steady state is the
// Darcy-Forchheimer balance independent of dt (explicit: dt < 2/(CPOR n (A + 2B|U|))).
void vrans_f::implicit_drag(lexer *p, fdm *a, int c, double w, field &f, field &U, field &V, field &W)
{
    if(p->B268!=1)
    return;
    
    if(c==0)
    ULOOP
    if(PORVAL1<1.0)
    f(i,j,k) /= 1.0 + w*p->dt*CPOR1*PORVAL1*drag_rate(p,a,0,U,V,W);
    
    if(c==1)
    VLOOP
    if(PORVAL2<1.0)
    f(i,j,k) /= 1.0 + w*p->dt*CPOR2*PORVAL2*drag_rate(p,a,1,U,V,W);
    
    if(c==2)
    WLOOP
    if(PORVAL3<1.0)
    f(i,j,k) /= 1.0 + w*p->dt*CPOR3*PORVAL3*drag_rate(p,a,2,U,V,W);
}
