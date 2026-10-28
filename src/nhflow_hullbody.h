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

#ifndef NHFLOW_HULLBODY_H_
#define NHFLOW_HULLBODY_H_

#include"increment.h"
#include<vector>
#include<fstream>
#include<functional>

class lexer;
class fdm_nhf;
class ghostcell;
class slice;
class sixdof_obj_nhflow;
class nhflow_thinbody;

using namespace std;

// The 6DOF hull as a sharp thin body (X 15 8, NHFLOW X 10 1/2): a client of nhflow_thinbody.
// The sides of the cell centres and nodes come from the body level set d->FB, the links between
// cells on opposite sides are blocked (wall fluxes, wall velocity in the projection, link mobility),
// the wall velocity is the normal part (u_b . n) n of the rigid-body velocity at the face (n the hull
// normal from the level set), the loads are the outer pressure on the blocked links, each face weighted
// with (n_face . n) n so that the loads and the wall fluxes do the same work. No direct forcing. See
// nhflow_hullbody.cpp.
//
//   stage()   every RK stage after the body advance and ray cast (nhflow_forcing::forcing)
//   forces()  from sixdof_obj_nhflow::hydrodynamic_forces_nhflow (after the projection)

class nhflow_hullbody : public increment
{
public:
    nhflow_hullbody(lexer*, fdm_nhf*, ghostcell*, sixdof_obj_nhflow*);
    virtual ~nhflow_hullbody();

    void stage(lexer*, fdm_nhf*, ghostcell*, double, double*, double*, double*, slice&);
    void forces(lexer*, fdm_nhf*, ghostcell*, slice&, bool);
    void print(lexer*, fdm_nhf*, ghostcell*, const double* = nullptr);
    
    // mesh refinement (nhflow_amr): every grid holds its own hull (the patches: amr_patch, no start message); the
    // loads are summed over the leaf part of the grids (loads, keep: the finest grid at the face centre) and then
    // given to the body by the hull of level 0 (finish_loads)
    void loads(lexer*, fdm_nhf*, slice&, const std::function<bool(double,double)>*, double*, double*);
    void finish_loads(lexer*, fdm_nhf*, ghostcell*, const double*, bool, const double*);
    nhflow_thinbody *thinbody() const { return tb; }
    bool amr_patch = false;
    // regrid: the state of column (ii,jj) from an old patch to a new one (the sides of the last stage with their
    // hysteresis, the state of the thin body, the face corrections of the Rhie-Chow continuity flux)
    int amr_nstate(lexer*) const;
    void amr_get(lexer*, fdm_nhf*, int, int, double*);
    void amr_set(lexer*, fdm_nhf*, int, int, const double*);

private:
    void body_velocity(double, double, double, double&, double&, double&) const;
    bool hull_normal(lexer*, fdm_nhf*, int, int, double&, double&, double&) const;

    sixdof_obj_nhflow *obj;
    nhflow_thinbody *tb;
    double *side0C;             // cell sides of the last stage (hysteresis)
    double *side0N;             // node sides of the last stage
    double hyst;                // hysteresis distance
    double Rn;                  // hydraulic resistance of the blocked links [m/s]
    int nbx,nby,nbz;            // blocked links on this rank (last stage)
    int nlink;                  // global
    bool told,marked;
    double Fp[6];               // last loads: Fx, Fy, Fz, Mx, My, Mz
    double tstage,tsum;

    ofstream out;
};

#endif
