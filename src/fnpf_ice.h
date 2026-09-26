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

#ifndef FNPF_ICE_H_
#define FNPF_ICE_H_

#include"increment.h"
#include"slice4.h"
#include"fnpf_ice_floe.h"
#include<vector>
#include<array>
#include<fstream>
#include<mpi.h>

class lexer;
class fdm_fnpf;
class ghostcell;
class slice;
class ice_contact;

using namespace std;

// REEF3D::FNPF rigid ice floes, compliant-lid coupling (A 380).
//
// Each floe acts on the free surface through the pressure
//
//   p = phi * max(0, rho_i*g*h + k*(z_s - z_b) + c*(deta/dt - dz_b/dt)),   z_s = wd + eta,
//
// added to the dynamic FSBC as dFi/dt = ... - p/rho_w (same slot as the wind pressure).
// phi is the footprint coverage of the cell, z_b(x,y) the height of the floe bottom plane.
// The coverage is tapered smoothly over +-A389 cells across the floe edge (area preserving): a sharp
// pressure edge makes the free surface a step of the floe draft, which the single-valued sigma-grid
// surface does not tolerate once the floe moves.
// The weight term makes the static draft exact, the unilateral spring-dashpot k = A381*rho_w*g
// enforces rigid-body kinematics up to O(1/A381) and lets floes lift off.
// The same pressure, integrated over the tilted bottom, plus ice-water skin drag gives the
// hydrodynamic force and moment on the floe: heave/roll/pitch from the lid, surge/sway from the
// bottom slope and the drag, yaw from the drag. Radiation damping and added mass come from FNPF.
//
// Step:  prestep   footprint coverage phi from the floe state at t^n (frozen over the step)
//        RK stages lid_forcing: z_b, dz_b/dt from the stage floe state, lid pressure with the stage eta
//                  and stage deta/dt (store_etat), stage loads on the floes, and the floe 6DOF state
//                  advanced with the same RK scheme as the free surface (A 310: SSP-RK3 or RK4).
//                  Integrating floes and surface as one explicit system avoids the added-mass type
//                  instability of a staggered coupling with a stiff lid and light floes.
//        poststep  planar non-smooth contact (ice_contact) as a velocity/position correction, output
//
// Floes are global: every rank sums its own cells, one MPI_Allreduce per stage, every rank advances
// the floes, rank 0 solves the contact and broadcasts the state.

class fnpf_ice : public increment
{
public:
    fnpf_ice(lexer*, fdm_fnpf*, ghostcell*);
    virtual ~fnpf_ice();

    void ini(lexer*, fdm_fnpf*, ghostcell*);
    void prestep(lexer*, fdm_fnpf*, ghostcell*);
    void store_etat(lexer*, fdm_fnpf*, slice&);
    void lid_forcing(lexer*, fdm_fnpf*, slice&, slice&);
    void poststep(lexer*, fdm_fnpf*, ghostcell*);
    void timestep(lexer*, fdm_fnpf*, ghostcell*);

private:
    // input and geometry
    void read(lexer*, ghostcell*);
    void read_file(lexer*);
    void add_polygon(lexer*, vector<double>&, vector<double>&, double, int);
    void floe_field(lexer*, double, double, double, double, double, double, double, double, int, int, int);
    void geometry(fnpf_ice_floe&);
    void place(lexer*, fnpf_ice_floe&, double, double, double, double, double, double);
    void serialize_geometry(vector<double>&);
    void deserialize_geometry(lexer*, vector<double>&);

    // footprint and lid
    void kinematics(fnpf_ice_floe&);
    void footprint(lexer*, fdm_fnpf*);
    double zbottom(const fnpf_ice_floe&, double, double) const;
    double zbottom_t(const fnpf_ice_floe&, double, double, double) const;
    double coverage(lexer*, const fnpf_ice_floe&, int, int) const;
    double edge_distance(const fnpf_ice_floe&, double, double) const;
    double ramp(double, double) const;
    void stage_forces(lexer*, fdm_fnpf*, slice&, slice&);

    // dynamics
    void derivative(const fnpf_ice_floe&, double*) const;
    void rk_stage(lexer*);
    void contact(lexer*);
    void broadcast_state(lexer*, ghostcell*);
    static void get_state(const fnpf_ice_floe&, double*);
    static void set_state(fnpf_ice_floe&, const double*);

    // breaking (A 390)
    struct split
    {
        int f;              // floe index
        int mech;           // 1 flexural, 2 contact splitting
        double nbx,nby,s;   // cut line in the body frame: nb.x = s
        double val;         // stress [Pa] or contact force [N]
        double lim;         // strength or splitting load
    };
    void breaking_moments(lexer*);
    void breaking_decide(lexer*);
    void breaking_apply(lexer*, ghostcell*);
    int split_floe(lexer*, size_t, double, double, double, int);
    static void clip_halfplane(const vector<double>&, const vector<double>&, double, double, double, double, vector<double>&, vector<double>&);
    static double chord(const vector<double>&, const vector<double>&, double, double, double);
    int breakflag,ndir,noff;
    double sigf,KIC,Csplit,Dmin;
    vector<double> Mcut;
    vector<int> cmap;       // contact body index -> floe index, last contact solve
    vector<split> splits;
    ofstream breakout;
    int nbreak;

    // output
    void print_ini(lexer*);
    void print(lexer*, fdm_fnpf*, ghostcell*);
    void print_vtp(lexer*, int);
    void print_log(lexer*);

    // helpers
    static void quat_to_matrix(const double*, double (*)[3]);
    static void quat_to_euler(const double*, double&, double&, double&);
    static void quat_normalize(double*);
    static void convex_hull(vector<double>&, vector<double>&);

    vector<fnpf_ice_floe> floe;
    int nfloe,nobst;

    // footprint cells, rebuilt every step
    struct lidcell
    {
        int i,j,f;
        double phi,area,xc,yc;
        double off;     // bottom offset in the taper and in shared cells: draft_f - sum_g phi_g*draft_g
        double pl;      // lid pressure of the last RK stage (for the bending moments)
    };
    vector<lidcell> cell;

    slice4 etat,dtot;
    ice_contact *pcontact;
    MPI_Comm comm;
    ghostcell *pgc_;

    // RK state of the floes: x[3], q[4], v[3], w[3]
    static constexpr int NY = 13;
    vector<array<double,NY>> Yn,D1,D2,D3;
    int stage,nstage,stagewarn;

    double rhow,g,wd;
    double alpha,zeta,Cd;
    int nsub;
    double taper,dxmean;
    int is2D;
    double dtfac;

    // output
    double printtime;
    int printcount;
    ofstream logout,obstout;
};

#endif
