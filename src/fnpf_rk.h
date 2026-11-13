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

#ifndef FNPF_RK_H_
#define FNPF_RK_H_

#include"fnpf.h"
#include"fnpf_ini.h"
#include"fnpf_sigma.h"
#include"slice4.h"

class fnpf_laplace;
class fnpf_laplace_cds2;
class fnpf_fsf;
class fnpf_body;
class fnpf_amr;

using namespace std;

// The building blocks of an explicit Runge-Kutta step of FNPF, shared by fnpf_RK3 and fnpf_RK4.
// A stage consists of
//   stage_tendency  the kinematic and dynamic free-surface conditions of the current state: ek, fk
//   (the caller combines the stage values eta, Fifsf from ek, fk)
//   stage_surface   damping, relaxation zones, wetting-drying and coastline, breaking, the free-surface
//                   discretisation, the sigma grid and the boundary conditions of Fi for the stage values
//   stage_solve     the Laplace solve for Fi and the free-surface velocities
// and step_velocity forms the 3D velocities at the end of the step.
class fnpf_rk : public fnpf_ini, public fnpf_sigma
{
public:
	fnpf_rk(lexer*, fdm_fnpf*, ghostcell*);
	virtual ~fnpf_rk();

    void inidisc_step1(lexer*, fdm_fnpf*, ghostcell*, ioflow*, solver*) override;
    void inidisc_step2(lexer*, fdm_fnpf*, ghostcell*, ioflow*, solver*) override;
    fnpf_body* body() override final {return pbody;}

protected:
    // kinematic FSBC -> ek, dynamic FSBC with eta -> fk; then the patches' tendencies and the body stage s
    void stage_tendency(lexer*, fdm_fnpf*, ghostcell*, solver*, slice &eta, slice &ek, slice &fk, int s);

    // the stage values e, f: damping with weight alpha, relaxation zones (targets at time t, weight w) and coastline,
    // breaking (weight balpha; zones: once per step, in the first stage), fsfdisc, sigma grid, boundary conditions
    void stage_surface(lexer*, fdm_fnpf*, ghostcell*, ioflow*, slice &e, slice &f, double alpha, double t, double w,
                       double balpha, int s, bool zones);

    // the Laplace solve of the stage; svq > 0: start value plap0->start_value(svq,svc) (A 321), keep: keep the
    // coarse levels of REEFMG of the last solve
    void stage_solve(lexer*, fdm_fnpf*, ghostcell*, solver*, slice &e, slice &f, int svq, double svc, bool keep);

    // the 3D velocities at the end of the step
    void step_velocity(lexer*, fdm_fnpf*, ghostcell*);

    int gcval;
    int gcval_eta,gcval_fifsf;

    slice4 en;      // eta at t^n, for the breaking onset criterion

    fnpf_laplace *plap;
    fnpf_laplace_cds2 *plap0;
    fnpf_amr *pamr = nullptr;
    fnpf_fsf *pf;
    fnpf_body *pbody;
};

#endif
