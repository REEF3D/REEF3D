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
// N 59 1 (with the level set, the momentum of N 40 5, 6): the momentum per volume rho u carried with
// the mass flux. The density of the control volume of the velocity is that of the projection (the
// face density of the level set, density_f), its face states come from the same Godunov
// reconstruction as those of the velocity, and on every face of the control volume
//   m = U rho_face (the mass flux),   the momentum flux m u_face,
//   rho* = rho - dt [div(m) - rho div(U)],   (rho u)* = rho u - dt [div(m u) - rho u div(U)],
//   u* = (rho u)*/rho*  (div(U): the divergence of the advecting velocities of the control volume,
//   taken out of both, so that a constant density stays constant),
// so that mass and momentum move with the same fluxes (Rudman 1998, Raessi and Pitsch 2012, Nangia
// et al. 2019). rho* is only used for u*: the next step starts again from the density of the level
// set. The second-order update is taken where rho* keeps at least half the density of the control
// volume and u* stays within the range of the velocity of the control volume and its face states.
// Otherwise (a light control volume next to a heavy one, whose reconstructed outflow density can
// exceed its own; a large relative loss of mass with outflow states other than u, which would
// extrapolate u* = (rho u - m u_out)/(rho - m)) the same update with first-order upwind states,
// u* = u + dt sum m_in (u_in - u)/rho*, a convex combination; where that fails too, the advective
// form.
// N 59 2: the face density of the mass flux from the level set at the half step of its own transport,
// rho((phi_start + phi)/2) on the faces (validation 84), where the level set has already been moved
// with the advecting velocity of the update: N 40 5 with the level set kept before its last
// transport (phi_start, momentum_rk), the full step of N 40 6 (phi(n), phi(n+1)). The mass flux then
// moves the water as the level set (with N 63 its WENO stages) did, time-centred and without a
// reconstruction of the density; in the normal direction the faces of the control volume are the
// centres of the cells, in the others the mean of the four cells of the edge. Elsewhere (the first
// step, a new patch, the half step of N 40 6) the Godunov states of N 59 1.
class godunov final : public convection, public increment
{
public:
	godunov(lexer*);
	virtual ~godunov();

	void start(lexer*,fdm*,field&,int,field&,field&,field&) override final;

    // N 59: the level set of the density at the start of the update (nullptr: a->phi) and the time
    // of the update (<= 0: dt); set by momentum_rk for the steps of N 40 6
    void consistent_set(field *phi, double dtu) {cphi=phi; cdt=dtu;}

private:
    // the velocities on the two faces of the control volume (i,j,k) of the variable ipol in x, y, z
    void vx(int,int,int,int,field&,double&,double&);
    void vy(int,int,int,int,field&,double&,double&);
    void vz(int,int,int,int,field&,double&,double&);

    // position of the centre of the control volume and its width in direction d
    double pos(lexer*,int,int,int) const;
    double width(lexer*,int,int,int) const;

    static double mc(double,double);

    // N 59: the density of the control volume (ii,jj,kk) of the variable ipol from the level set
    double rho_cv(lexer*,field&,int,int,int,int) const;

    // the stored values of the control volumes of the extended range
    std::vector<double> sf, sx, sy, sz, uc, vc, wc, tx, ty, tz;

    // N 59: the density of the control volumes, its slopes and transverse terms
    std::vector<double> rf, rx, ry, rz, rtx, rty, rtz;
    field *cphi = nullptr;
    double cdt = 0.0;
    int i0, j0, k0, ni, nj, nk;
    inline int id(int ii, int jj, int kk) const {return ((ii-i0)*nj + (jj-j0))*nk + (kk-k0);}
};

#endif
