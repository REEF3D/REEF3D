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
};

class field;

class pic_f : public pic, public increment
{
public:
    pic_f(lexer*, ghostcell*);
    virtual ~pic_f();

    void ini(lexer*, fdm*, ghostcell*) override;
    void start(lexer*, fdm*, ghostcell*) override;

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

    // interpolation, also in the ghost cells (pic_f_ipol.cpp)
    double ipol(lexer*, field&, int, double, double, double) const;
    int cell(const double*, int, double) const;

    // parallel (pic_f_xchange.cpp)
    void xchange(lexer*, ghostcell*);

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
    int npd;                  // particles per cell and direction (F 83)
    int npc;                  // particles per cell

    double printtime;
    int printcount;
    long long seeded, lost, lost_total, sent, sent_total;
    std::ofstream logout;

    // values per particle in the gathered output: x y z u v w id rank
    static const int NG = 8;
    static const int NF = sizeof(pic_data)/sizeof(double);
};

#endif
