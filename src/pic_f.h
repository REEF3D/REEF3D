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

#ifndef PIC_F_H_
#define PIC_F_H_

#include"pic.h"
#include"increment.h"
#include<cstdint>
#include<fstream>
#include<vector>

// one particle; plain doubles, sent between the ranks as they are
struct pic_data
{
    double x,y,z;   // position
    double u,v,w;   // velocity (phase 0: the CFD velocity at the particle)
    double id;      // global id, the same for any decomposition: seeding cell * F 83^dim + sub-position
    double b[9] = {0,0,0,0,0,0,0,0,0};   // APIC (F 77 1, 2): the gradients of u, v, w at the particle
};

class field;
class field1;
class field2;
class field3;
class reini;
class ioflow;
class heat;
class concentration;
class fluid_update;

class pic_f : public pic, public increment
{
public:
    pic_f(lexer*, fdm*, ghostcell*, reini*, ioflow*, heat*, concentration*);
    virtual ~pic_f();

    void ini(lexer*, fdm*, ghostcell*) override;
    void start(lexer*, fdm*, ghostcell*) override;
    void step_begin(lexer*, fdm*, ghostcell*) override;
    bool flip() const override { return flipmode; }
    bool singlephase() const override { return single; }
    void air_velocity(lexer*, fdm*, ghostcell*, field&, field&, field&) override;
    void mask_convection(lexer*, fdm*, int, const double*, field&) override;

private:
    // seeding (pic_f_seed.cpp)
    void seed(lexer*, fdm*, ghostcell*);
    static double hash_uniform(uint64_t, uint64_t);

    // motion (pic_f_move.cpp)
    void move(lexer*, fdm*, ghostcell*);
    void velocity(lexer*, fdm*, double, double, double, double&, double&, double&);
    bool blocked(lexer*, fdm*, double, double, double);
    bool domain(lexer*, double&, double&, double&);
    void slide(lexer*, fdm*, pic_data&, double, double, double);
    void push_out(lexer*, fdm*, pic_data&);

    // interpolation, also in the ghost cells (pic_f_ipol.cpp)
    double ipol(lexer*, field&, int, double, double, double) const;
    int cell(const double*, int, double) const;
    int stencil(lexer*, int, double, double, double, int*, int*, int*, double*) const;
    int stencil_grad(lexer*, int, double, double, double, int*, int*, int*, double*, double*) const;
    void node(lexer*, int, int, int, int, double*) const;
    void ipol_grad(lexer*, field&, int, double, double, double, double&, double*) const;
    void affine(lexer*, pic_data&);   // APIC: the gradients of the water velocity at the particle

    // parallel (pic_f_xchange.cpp)
    void xchange(lexer*, ghostcell*);
    void halo(lexer*, ghostcell*, std::vector<double>&, bool vel=false);
    void sendrecv(ghostcell*, int, int, std::vector<double>&, std::vector<double>&, std::vector<double>&, int);

    // the free surface from the particles (pic_f_phi.cpp; F 81 1, 2)
    void phi_ini(lexer*, fdm*, ghostcell*);
    void phi_build(lexer*, fdm*, ghostcell*);
    void redistance(lexer*, fdm*, const std::vector<char>&);
    void support(lexer*, fdm*);
    void fraction_table(lexer*);
    double fraction_distance(double) const;
    int lid(lexer*, int, int, int) const;
    static double kernel(double d2, double R2) { const double q = 1.0 - d2/R2; return q>0.0 ? q*q*q : 0.0; }

    // the velocity of the water, extended into the air (pic_f_extend.cpp; F 81 1, 2)
    void water_velocity(lexer*, fdm*, ghostcell*);
    void extend(lexer*, fdm*, ghostcell*, field&, field&, int, int, bool);
    bool ballistic(lexer*, fdm*, const pic_data&) const;
    bool solid_distance(lexer*, fdm*, double, double, double, double&, double*);

    // the FLIP / PIC update of the particle velocities (pic_f_flip.cpp; F 81 1, 2)
    void g2p(lexer*, fdm*);

    // the redistribution of the particles (pic_f_redistribute.cpp; F 78, F 79)
    void redistribute(lexer*, fdm*, ghostcell*);

    // particles kept in the water (pic_f_reseed.cpp; F 86), and in the relaxation zones of the
    // wave generation and the numerical beach (B 98 2, B 99 1, 2) the water of the relaxed phi
    void reseed(lexer*, fdm*, ghostcell*, bool interior);
    bool inzone(lexer*, double, double);
    void zone_surface(lexer*, fdm*, ghostcell*, std::vector<double>&, const std::vector<double>&);
    void zone_stretch(lexer*, fdm*);
    std::vector<double> etap, etar, zbot;   // per column of the domain: surface of the particles, relaxed, bottom

    // the particles of the inflow in the ghost cells (pic_f_inflow.cpp; B 60)
    void inflow_seed(lexer*, fdm*, ghostcell*);

    // output (pic_f_print.cpp)
    void print(lexer*, ghostcell*);
    bool print_due(lexer*) const;
    void gather(lexer*, ghostcell*, std::vector<double>&);
    bool print_lagoon(lexer*, const std::vector<double>&);
    void print_vtp(lexer*, const std::vector<double>&);
    double print_interval(lexer*) const;

    std::vector<pic_data> P;

    double xloc_s,xloc_e,yloc_s,yloc_e,zloc_s,zloc_e;   // the subdomain
    double xglo_s,xglo_e,yglo_s,yglo_e,zglo_s,zglo_e;   // the domain
    int fmarg;                // the ghost cells of the fields (lexer::margin; increment::marge offsets the coordinates)
    int npd;                  // particles per cell and direction (F 83)
    int npc;                  // particles per cell
    bool surface;             // F 81 1, 2: the free surface (phi) from the particles; F 81 10: passive particles
    bool flipmode;            // F 81 1, 2: the momentum of the water from the particles
    bool single;              // F 81 1: single phase (PICIN), p = 0 at the surface, no air in the momentum
    bool apic;                // F 81 1, 2 with F 77 1, 2: affine velocities (Jiang et al. 2015)
    bool relaxzones;          // the relaxation zones of the waves (B 90 with B 98 2, B 99 1, 2)
    std::vector<char> zcol;   // per cell column of the subdomain: in a relaxation zone
    bool inflowsrc;           // B 60: particles in the ghost cells of the inflow (x-min)

    // the free surface from the particles
    reini *preini;
    ioflow *pflow;
    fluid_update *pupdate;
    int gcval_phi;
    field1 *uw, *uo, *km1;
    field2 *vw, *vo, *km2;
    field3 *ww, *wo, *km3;
    std::vector<double> ctab_d,ctab_c;   // water fraction across a flat surface (height above it in cells)
    double dfull;             // distance (in cells) of the edge of the band of the surface
    std::vector<double> cnorm,crad;   // per cell: kernel sum over the fluid, kernel radius
    long long added, removed, added_total, removed_total;
    double phitime;

    double printtime;
    int printcount;
    long long seeded, lost, lost_total, sent, sent_total;
    std::ofstream logout;

    // values per particle in the gathered output: x y z u v w id rank
    static const int NG = 8;
    static const int NF = sizeof(pic_data)/sizeof(double);
};

#endif
