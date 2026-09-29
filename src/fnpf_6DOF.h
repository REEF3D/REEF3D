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

#ifndef FNPF_6DOF_H_
#define FNPF_6DOF_H_

#include"increment.h"
#include"fnpf_body.h"
#include"slice4.h"
#include<vector>

class lexer;
class fdm_fnpf;
class ghostcell;
class solver;
class fnpf_laplace;
class fnpf_fsf;
class fnpf_bed_update;
class fnpf_fsf_update;
class sixdof_obj;
class slice;

using namespace std;

// Resolved rigid bodies in REEF3D::FNPF (X 10 1), driven by the shared 6DOF kernel.
//
// Geometry: vertical ray cast of the trimesh per sigma column -> body nodes (FBF),
// staircase Neumann faces in the Laplace matrix, body footprint where the
// free-surface node lies inside the body. In the footprint eta and Fifsf are a
// smooth (harmonic) extension of the surrounding free surface, so the sigma grid
// stays regular; they carry no physics there.
//
// Loads: psi = phi_t is solved with the same operator as phi,
//   psi_0: Dirichlet phi_t at the free surface (from the RK tendencies), Neumann
//          w x (w x r).n on the hull (m-terms neglected)
//   psi_j: Dirichlet 0, Neumann N_j (unit modes) -> added mass A_ij
// and the body is advanced with (M + A) a = F_0 + F_ext, stage-synchronous with
// fnpf_RK3 or fnpf_RK4. A is refreshed once per time step (first stage).
//
// Coupling to the time stepping only through the fnpf_body hooks: stage(), surface(),
// and a decorated Laplace solver (geometry before, body-band extrapolation after the
// phi solve). The psi solves use the undecorated solver.

class fnpf_6DOF : public fnpf_body, public increment
{
public:
    fnpf_6DOF(lexer*, fdm_fnpf*, ghostcell*);
    virtual ~fnpf_6DOF();
    
    // fnpf_body hooks
    void stage(lexer*, fdm_fnpf*, ghostcell*, solver*, fnpf_fsf*, slice&, slice&, int) override;
    void surface(lexer*, fdm_fnpf*, ghostcell*, slice&, slice&, int, int) override;
    fnpf_laplace* laplace(fnpf_laplace*) override;
    
    // used by the decorated Laplace solver around the phi solve:
    // body geometry on the new sigma grid, then extrapolation of the body band
    void pre_solve(lexer*, fdm_fnpf*, ghostcell*);
    void post_solve(lexer*, fdm_fnpf*, ghostcell*, double*);
    
    bool initialized;
    
private:
    void geometry(lexer*, fdm_fnpf*, ghostcell*);
    void extrapolate(lexer*, fdm_fnpf*, ghostcell*, double*);
    void ini(lexer*, fdm_fnpf*, ghostcell*);
    void forces(lexer*, fdm_fnpf*, ghostcell*, solver*, fnpf_fsf*, slice&, slice&, int);
    void motion(lexer*, fdm_fnpf*, ghostcell*, int);
    void footprint(lexer*, fdm_fnpf*, ghostcell*, slice&, slice&, int, int);
    
    void solve_psi(lexer*, fdm_fnpf*, ghostcell*, solver*, fnpf_laplace*, fnpf_fsf*, double*, slice&);
    void zero_face(lexer*, fdm_fnpf*);
    void exchange_face(lexer*, fdm_fnpf*, ghostcell*);
    
    vector<sixdof_obj*> fb_obj;
    int nbody;
    int gcval;
    int footcount;
    
    double *psi0;
    double *psi[6];
    double *zero;
    int *mark;
    
    slice4 foot,psiD,zeroslice;
    slice4 eta_ext,fi_ext;
    bool ext_ini;
    
    fnpf_bed_update *pbed;
    fnpf_fsf_update *pvel;
    fnpf_laplace *plap;
};

#endif
