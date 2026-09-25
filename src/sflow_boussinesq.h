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

#ifndef SFLOW_BOUSSINESQ_H_
#define SFLOW_BOUSSINESQ_H_

#include"increment.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;
class slice;
class solver2D;

using namespace std;

// Fully nonlinear Boussinesq equations of FUNWAVE-TVD (A 220 4)
// Shi, Kirby, Harris, Geiman & Grilli (2012), Ocean Modelling 43-44, 36-51;
// Chen (2006) with the conservative form of Shi et al.
//
//   eta_t + div(M) = 0,                 M = H (u_a + u2)
//   V_t + div(M M/H) + grad(g(eta^2 + 2 h eta)/2) = g eta grad(h) + psi + H R
//   V = H (u_a + V1')
//
//   u_a at the reference level z_a = zeta h + (1+zeta) eta, zeta = -0.53
//   A = div(h u_a), B = div(u_a)
//   u2  = (z_a^2/2 - (h^2 - h eta + eta^2)/6) grad(B) + (z_a + (h - eta)/2) grad(A)
//   V1' = z_a^2/2 grad(B) + z_a grad(A) - grad(eta^2/2 B + eta A)
//   V1''= grad(eta_t (A + eta B))
//   V2  = grad((z_a - eta)(u_a.grad)A + (z_a^2 - eta^2)/2 (u_a.grad)B + (A + eta B)^2/2)
//   V3  = omega_0 i_z x u2 + omega_2 i_z x u_a
//   psi = eta_t (V1' - u2) + H (u_a.grad u2 + u2.grad u_a - V1'' - V2 - V3)
//
// In SFLOW: UH,VH hold V (prognostic, RK combination), b->UA,VA hold u_a,
// b->MX,MY hold M (flux variable for HLL), b->U,V = M/H (depth-averaged velocity).
// u_a is recovered from V with a line-implicit solve per component (x-x and y-y
// derivatives implicit, cross derivatives lagged). Dispersion is switched off
// (nonlinear shallow water equations) in breaking cells, for eta/h >= 0.8, in
// strong drawdown (eta/h <= -0.5), next to the shoreline and next to in- and
// outflow boundaries. The switch m is widened
// by one cell, smoothed over two cells and enters as A -> m A, B -> m B with compact
// face-based gradients, so that u2 and V1' share one symmetric operator (stable
// for any spatial variation of m). Where m changes, V is kept and u_a follows.

class sflow_boussinesq : public increment
{
public:
	sflow_boussinesq(lexer*, fdm2D*, ghostcell*);
	virtual ~sflow_boussinesq();
    
    void mask_update(lexer*, fdm2D*, ghostcell*, slice&);
    void source(lexer*, fdm2D*, ghostcell*, slice&);
    void invert(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, slice&, slice&);
    void flux(lexer*, fdm2D*, ghostcell*, slice&);
    void forward(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, int);
    void save(lexer*, fdm2D*);
    
private:
    void cellterms(lexer*, fdm2D*, ghostcell*);
    void operators(lexer*, fdm2D*, ghostcell*);
    void explicitterms(lexer*, fdm2D*, ghostcell*);
    
    double za(lexer*, fdm2D*, int, int);
    double ddx(lexer*, slice&);
    double ddy(lexer*, slice&);
    double lx(lexer*, fdm2D*, slice&);
    double ly(lexer*, fdm2D*, slice&);
    
    void coef_x(lexer*, fdm2D*, double&, double&, double&);
    void coef_y(lexer*, fdm2D*, double&, double&, double&);
    
    const double zeta, breakratio;
    int q;
    
    slice4 A,B,Ax,Ay,Bx,By;
    slice4 U4,V4,U1p,V1p,Cx,Cy;
    slice4 T1,T2,etat;
    slice4 vy,hvy,ux,hux;
    slice4 mask,io;
    slice4 ua_n,va_n,f;
};

#endif
