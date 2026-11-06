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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef GODUNOV_H_
#define GODUNOV_H_

#include"increment.h"
#include"convection.h"
#include<vector>

class lexer;
class fdm;
class field;

using namespace std;

// The convection of the one-step scheme (N 40 5, P9 of the CFD AMR performance plan): an unsplit
// second-order Godunov method (corner transport upwind, the predictor of Bell, Colella and Glaz).
// The value on every face of a control volume is extrapolated from the upwind side to the half time
// step: f + (dx/2 - dt/2 u) f_x - dt/2 (v f_y + w f_z), with the slopes f_x of the monotonised
// central limiter and the transverse derivatives upwind. The convection term is then
// -u (f_e - f_w)/dx - v (f_n - f_s)/dy - w (f_t - f_b)/dz with the velocities of the control volume.
// One evaluation per time step; stable up to a Courant number of 1 in two dimensions.
// Without the forces of the momentum in the face states (the half step dt/2 F of Bell, Colella and
// Glaz): in the air next to the interface the pressure gradient divided by the air density gives
// accelerations of hundreds of m/s^2, which made spurious air jets in the dam break (validation 73);
// on the waves the gauges are the same with and without.
// N 56 (default 2): the momentum first-order upwind in the air further than N56 cells from the free
// surface. The air away from the interface is driven by it and its convection need not be second
// order; the dissipation keeps the spurious air velocities down (validation 74: the dam break with
// refinement after the first impact).
class godunov final : public convection, public increment
{
public:
	godunov(lexer*);
	virtual ~godunov();

	void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

private:
    // the velocities on the two faces of the control volume (i,j,k) of the variable ipol in x, y, z
    void vx(int,int,int,int,field&,double&,double&);
    void vy(int,int,int,int,field&,double&,double&);
    void vz(int,int,int,int,field&,double&,double&);

    // position of the centre of the control volume and its width in direction d
    double pos(lexer*,int,int,int) const;
    double width(lexer*,int,int,int) const;

    static double mc(double,double);

    // the stored values of the control volumes of the extended range
    std::vector<double> sf, sx, sy, sz, uc, vc, wc, tx, ty, tz;
    int i0, j0, k0, ni, nj, nk;
    inline int id(int ii, int jj, int kk) const {return ((ii-i0)*nj + (jj-j0))*nk + (kk-k0);}
};

#endif
