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

/*--------------------------------------------------------------------
dem_f: REEF3D adapter of the DEM core.

Coupling modes (ctrl E 11):
 0 dry        no fluid forces
 1 unresolved particles smaller than the grid: drag (Haider-Levenspiel with
              Di Felice voidage correction, implicit in the DEM step), buoyancy from
              volume quadrature, added mass and fluid acceleration; the reaction is
              spread to the fluid momentum with a compact kernel
 2 resolved   particles larger than the grid: direct forcing of the fluid velocity
              inside the particle level set, hydrodynamic force from the volume
              integral of the forcing (Uhlmann 2005)
 3 hybrid     per particle, resolved if d_eq/dx >= E 12, else unresolved

Domain decomposition (dem_f_mpi.cpp): small particles are owned by the rank
holding their centroid and appear as ghosts on neighbouring ranks; large
particles (E 24) are replicated. Fluid data and forces are evaluated on the rank
that owns the evaluation point and returned to the particle's owner.
--------------------------------------------------------------------*/

#ifndef DEM_F_H_
#define DEM_F_H_

#include"dem.h"
#include"dem_core.h"
#include"increment.h"
#include<vector>
#include<string>
#include<unordered_map>

class lexer;
class fdm;
class fdm_nhf;
class ghostcell;
class field;
class field1;
class field2;
class field3;
class field4a;
class slice;

using namespace std;

// serialisation helpers for the particle exchange
struct dem_writer
{
    vector<double> &b;
    explicit dem_writer(vector<double> &buf) : b(buf) {}
    void d(double v) {b.push_back(v);}
    void v(const dem_vec &x) {b.push_back(x(0)); b.push_back(x(1)); b.push_back(x(2));}
    void q(const dem_quat &x) {b.push_back(x.w()); b.push_back(x.x()); b.push_back(x.y()); b.push_back(x.z());}
    void m(const dem_mat &x) {for(int r=0; r<3; ++r) for(int c=0; c<3; ++c) b.push_back(x(r,c));}
};

struct dem_reader
{
    const vector<double> &b;
    size_t pos;
    dem_reader(const vector<double> &buf, size_t start) : b(buf), pos(start) {}
    double d() {return b[pos++];}
    int i() {return int(llround(b[pos++]));}
    dem_vec v() {dem_vec x(b[pos],b[pos+1],b[pos+2]); pos+=3; return x;}
    dem_quat q() {dem_quat x(b[pos],b[pos+1],b[pos+2],b[pos+3]); pos+=4; return x;}
    dem_mat m() {dem_mat x; for(int r=0; r<3; ++r) for(int c=0; c<3; ++c) x(r,c)=b[pos++]; return x;}
};

class dem_f final : public dem, public increment
{
public:
    dem_f(lexer*, ghostcell*);
    virtual ~dem_f();

    void start_cfd(lexer*, fdm*, ghostcell*) override final;
    void start_nhflow(lexer*, fdm_nhf*, ghostcell*) override final;

    void forcing_cfd(lexer*, fdm*, ghostcell*, int, double, field&, field&, field&, bool) override final;
    void forcing_nhflow(lexer*, fdm_nhf*, ghostcell*, int, double, double*, double*, double*, slice&, bool, bool) override final;

private:
    // setup
    void read_input(lexer*, ghostcell*);
    void ini(lexer*, ghostcell*);
    void ini_cfd(lexer*, ghostcell*);
    void ini_nhflow(lexer*, ghostcell*);

    // time step
    void step(lexer*, ghostcell*);
    void set_loads(lexer*, ghostcell*);
    void sync(lexer*, ghostcell*);
    void deactivate(lexer*);
    int substeps(lexer*, ghostcell*);

    // point ownership in the domain decomposition
    bool owns(lexer*, double, double, double);

    // CFD coupling
    void fluid_cfd(lexer*, fdm*, ghostcell*);
    void feedback_cfd(lexer*, fdm*, ghostcell*);
    void walls_cfd(lexer*, fdm*, ghostcell*, double, vector<dem_contact>&);
    double wallphi_cfd(lexer*, fdm*, double, double, double);

    // NHFLOW coupling
    void fluid_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void feedback_nhflow(lexer*, fdm_nhf*, ghostcell*);
    void walls_nhflow(lexer*, fdm_nhf*, ghostcell*, double, vector<dem_contact>&);
    double wallphi_nhflow(lexer*, fdm_nhf*, double, double, double);

    // domain decomposition
    void decomp_ini(lexer*, ghostcell*);
    double halo(const dem_body&);
    bool in_box(int, const dem_vec&);
    int find_owner(const dem_vec&);
    double boxdist(int, const dem_vec&);
    void p2p(ghostcell*, vector<vector<double>>&, vector<vector<double>>&);
    void rebuild_maps();
    void erase_ghosts();
    void ghost_exchange(lexer*, ghostcell*);
    void pack_body(vector<double>&, const dem_body&);
    void unpack_body(dem_reader&, dem_body&);
    void migrate(lexer*, ghostcell*);
    void reduce_owner(ghostcell*, vector<double>&, int, bool, bool);
    void owner_to_ghosts(ghostcell*, vector<double>&, int);
    double sync_solver(ghostcell*, int, double, bool);
    void mass_split(ghostcell*);
    void route_walls(lexer*, ghostcell*, const vector<double>&, vector<dem_contact>&);
    void global_stats(ghostcell*, double&, double&);
    void gather_output(ghostcell*, vector<dem_body>&);
    void run_step(lexer*, ghostcell*, const dem_core::wallfunc&);

    int myrank = 0, nproc = 1;
    vector<double> boxes;
    double p_gmax[3] = {1.0e30,1.0e30,1.0e30};
    double rmax0 = 0.0;
    vector<int> partners, partner_index;
    vector<vector<int>> ghostdest;
    unordered_map<int,int> gid2loc;
    vector<int> tier1idx;
    int nglobal = 0;
    int ncolors = 1, mycolor = 0;

    // common
    void drag(lexer*, int, double, double, double);
    double heaviside(double, double);

    // grid point relative to a particle centre; in 2D the y offset is ignored (x-z slab)
    dem_vec relpos(lexer*, double, double, double, const dem_vec&);
    double kernel(double, double);
    double kradius(int);        // momentum kernel radius of a particle
    double vradius(int);        // solid fraction kernel radius
    void cellrange(lexer*, int, double, int&, int&, int&, int&, int&, int&);

    // output
    void print(lexer*, ghostcell*);
    void print_vtp(lexer*, const vector<dem_body>&);
    void print_state(lexer*, const vector<dem_body>&);
    void print_log(lexer*, ghostcell*);

    dem_core core;

    int solver;                 // 6 CFD, 5 NHFLOW
    int coupling;
    bool initialized;
    int nb;

    string infile;
    double node_spacing_factor;
    int sdf_res, quad_res;
    int gridwalls;

    // resolved: hydrodynamic force from the forcing, fluid mass inside the particle
    void combine_stages(lexer*, ghostcell*, int, double);
    bool rkwarn = false;
    double dxs = 1.0;           // horizontal grid length scale (NHFLOW: DXM mixes in the sigma spacing)

    // resolved: fluid momentum and angular momentum inside the particles (Kempe & Froehlich 2012)
    void internal_cfd(lexer*, fdm*, ghostcell*);
    void internal_nhflow(lexer*, fdm_nhf*, ghostcell*);
    double dt_old;


    // CFD source fields and solid fraction
    field1 *Sx;
    field2 *Sy;
    field3 *Sz;
    field4a *ALPHA;

    // NHFLOW source arrays and solid fraction
    double *SX,*SY,*SZ,*ALPHAV;
    // NHFLOW: point-implicit drag coefficient and the velocities it was evaluated with (E 27),
    // solid fraction for the porosity (E 28)
    double *SD,*UFB,*VFB,*WFB,*ALPHAP;
    int multipoint, fimplicit, pormode, fixeddrag, porinit;

    // NHFLOW: 4-point Peskin kernel spreading of a point load onto the sigma grid (E 26)
    double peskin(double) const;
    template<class F> void peskin_spread(lexer*, fdm_nhf*, const dem_vec&, F&&);
    bool fluidcoupled(const dem_body&) const;

    // parameters
    double hybrid_ratio, Ca, hs_factor, travel, kernel_cells, void_factor;
    int fluidacc;
    int minsub, nsub;

    // output
    double printtime;
    int printcount;
    double steptime;
    int maxiter_used;
};

#endif
