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

#ifndef SFLOW_PJM_QUAD_H_
#define SFLOW_PJM_QUAD_H_

#include"sflow_pressure.h"
#include"increment.h"
#include"slice4.h"

using namespace std;

// depth-integrated non-hydrostatic pressure, quadratic vertical profile
// (Jeschke et al. 2017), collocated for the HLL scheme.
//   press = depth-averaged non-hydrostatic pressure q
//   bottom value q_b = 1.5 q + rho*h*phi/4,  phi = Dw_b/Dt, w_b = -u.grad(d)
//     phi = -(Du/Dt).grad(d) - (u^2 d_xx + 2uv d_xy + v^2 d_yy)
//   UH -= alpha*dt/rho * ( d(h q)/dx - q_b d(d)/dx )
//   WH += alpha*dt/rho * q_b
//   constraint: h div(u) + 2(w + u.grad(d)) = 0
// linear dispersion: omega^2 = g h k^2 / (1 + k^2 h^2/3)
//
// A 220 3: improved dispersion (Chazel, Lannes & Marche 2011, J. Sci. Comput. 48)
//   (I + a T)(du/dt + g grad eta) = T(g grad eta),  T: Green-Naghdi dispersive operator,
//   a = A 224 (default 1.159, a=1 recovers A 220 2)
//   -> projected pressure q = a*q_GN, bottom pressure q_b = 1.5 q/a + rho h phi/4 in the
//      w-equation and a*q_b in the momentum, plus the explicit source
//      d(UH)/dt += B g grad(h^3 div(grad(eta))),  B = (a-1)/3
//      discretised with central differences only (consistent with the collocated
//      projection), explicit: dt <= 1.8 dx^2/sqrt(g B h^3), see sflow_etimestep
//   linear dispersion: omega^2 = g h k^2 (1 + B k^2 h^2) / (1 + a k^2 h^2/3)

class sflow_pjm_quad final : public sflow_pressure, public increment
{
public:
    sflow_pjm_quad(lexer*, fdm2D*, ghostcell*, patchBC_interface*);
	virtual ~sflow_pjm_quad();

	void start(lexer*, fdm2D*, ghostcell*, solver2D*, ioflow*, slice&, slice&, slice&, slice&, slice&, slice&, double) override final;
	void upgrad(lexer*, fdm2D*, slice&) override final;
	void vpgrad(lexer*, fdm2D*, slice&) override final;
    void wpgrad(lexer*, fdm2D*, slice&) override final;

private:
    void quad_calc(lexer*, fdm2D*, ghostcell*, slice&, slice&, slice&, slice&, slice&, double);
    void rhs(lexer*, fdm2D*, slice&, double);
    void poisson(lexer*, fdm2D*, slice&, double);
    void ucorr(lexer*, fdm2D*, slice&, slice&, double);
	void vcorr(lexer*, fdm2D*, slice&, slice&, double);
	void wcorr(lexer*, fdm2D*, slice&, slice&, double);
    int active(lexer*, fdm2D*);
    
	double starttime;
    int gcval_press,q;
	double solvtime,ptime;
    const double cb;
    double adisp,cw,Bdisp;
    
    slice4 phi,Uest,Vest,Ld,Gx,Gy;
    
    patchBC_interface *pBC;
    ghostcell *pgc;
};

#endif
