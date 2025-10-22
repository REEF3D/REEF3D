/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Authora: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#ifndef PARTRES_H_
#define PARTRES_H_

#include"increment.h"
#include"part.h"
#include"slice4.h"
#include"field4a.h"
#include"boundarycheck.h"

class lexer;
class fdm;
class ghostcell;
class sediment_fdm;
class turbulence;
class part;
class vrans;

using namespace std;

class partres : public increment
{
public:
    partres(lexer*, ghostcell*);
    virtual ~partres() = default;

    void move_RK2(lexer*, fdm*, ghostcell*, sediment_fdm*, turbulence*);

    void update(lexer*, fdm*, ghostcell*, sediment_fdm*, field&, field&);

    void timestep(lexer*, ghostcell*);

    void seed_topo(lexer*, fdm*, ghostcell*, sediment_fdm*);

    void print_particles(lexer*,sediment_fdm*);
private:
    void advec_plain(lexer*, fdm*, part&, sediment_fdm*, turbulence*,
                        double*, double*, double*, double*, double*, double*,
                        double&, double&, double&, double);
    void advec_mppic(lexer*, fdm*, part&, sediment_fdm*, turbulence*,
                        double*, double*, double*, double*, double*, double*,
                        double&, double&, double&, double);

    // drag
    double drag_model(lexer *, double, double, double, double);

    void count_particles(lexer*, fdm*, ghostcell*, sediment_fdm*);

    void stress_tensor(lexer*, ghostcell*, sediment_fdm*);
    void stress_gradient(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void pressure_gradient(lexer*, fdm*, ghostcell*, sediment_fdm*);
    void cellSum_update(lexer*, ghostcell*, sediment_fdm*, int);
    void cellSum_full_update(lexer*, ghostcell*, sediment_fdm*, int);

    void bedchange(lexer*, fdm*, ghostcell*, sediment_fdm*, int);
    void bedchange_update(lexer*, ghostcell*, int);

    void boundcheck(lexer*, int);

    part P;

    slice4 bedch;

    field4a Tau,Ts;
    field4a cellSum;

    // relax
    void relax_ini(lexer*);
    void relax(lexer*, ghostcell*, sediment_fdm*);
    double rf(lexer*, double, double);
    double r1(lexer*, double, double);
    double distcalc(lexer*, double , double, double , double, double);

    void print_vtp(lexer*,sediment_fdm*);
    void pvtp(lexer*,int);

    boundarycheck boundaries;

    int printcount;
    double printtime;

    field4a dPx,dPy,dPz;
    field4a dTx,dTy,dTz;

    double *tan_betaQ73,*betaQ73,*dist_Q73;

    bool timestep_ini = true;
};

#endif
