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

#ifndef HYDRO_SFLOW_H_
#define HYDRO_SFLOW_H_

#include"increment.h"
#include"slice4.h"
#include"lexer.h"
#include"hydro_raster.h"
#include<fstream>
#include<functional>
#include<string>
#include<vector>

class fdm2D;
class ghostcell;
class slice;

/*--------------------------------------------------------------------
REEF3D::HYDRO, milestone M1: rain-on-grid, point sources, open-water
evaporation, infiltration and a Manning raster for SFLOW (A 10 2), level 0
of Cartesian and curvilinear (G 50) grids. Not with mesh refinement or
multi-block grids yet (stopped with a message).

Sources, in the continuity of every stage (explicit: sflow_eta::update;
semi-implicit: the right-hand side of the eta system and the depth update
of sflow_momentum_SI), on wet AND dry cells, evaluated once per step at
the time of the start of the step:
  R 10 1    rain R 11 [mm/h], uniform and constant
  R 10 2    rain time series rain.dat: lines t [s] rate [mm/h]
  R 10 3    rain rasters: rain_series.dat, lines t [s] file.asc, every
            file an ESRI ASCII grid of the rate [mm/h] in world coordinates
            (hydro_raster.h); a cell takes the value of the raster cell its
            centre lies in, 0 outside the raster and at NODATA
  R 12      in time: 0 a value holds until the next time (rain totals of
            an interval, default), 1 linear; before the first time the
            first value, after the last time the last value
  R 40 ID x y Q   point source of Q [m3/s] in the cell that holds (x,y)
            (world coordinates as the gauges P 50); repeatable
  R 41 ID   point source ID takes Q(t) from hydrograph_src_<ID>.dat
            (lines t [s] Q [m3/s], linear, clamped); repeatable
  Q > 0 enters the cell without momentum. Q < 0 is an abstraction (an
  intake), taken with the losses below.

Losses, after the continuity of the step (operator splitting), limited
to the water above the wet/dry depth A 244 (the last A 244 of a cell is
never taken, so the losses never dry a cell below the film that wetdry
keeps; with losses A 244 acts as a depression storage):
  R 20      open-water evaporation [mm/day] where there is water above A 244
  R 30 1    infiltration at the constant rate R 31 [mm/h]
  R 30 2    Green-Ampt, R 32 Ks [mm/h] psi [m] dtheta [-]: the potential
            infiltration of the step from the cumulative infiltration F of
            the cell, F1 - F0 - P ln((F1 + P)/(F0 + P)) = Ks dt, P = psi
            dtheta (exact for a ponded cell), the cell takes what it has
  Shares: when the water is short, evaporation, infiltration and
  abstraction are reduced in proportion.

Wetting (R 10 or R 40): a dry cell gets wet by its own depth (above
A 244 + 1e-6 m) as well, not only from a higher wet neighbour (A 243 2),
and an isolated wet cell deeper than 2 A 244 is not dried (a puddle that
rain or a source fills); see sflow_eta::wetdry. Use A 243 2: with A 243 1
the rain on a dry cell below 1e-6 m per step is reset with the film.

Roughness R 50 1: Manning n from manning.asc (ESRI ASCII, world
coordinates), ks = (20 n)^6 so that n = ks^(1/6)/20 of sflow_rough_manning
and the implicit friction of sflow_momentum_SI is that n; cells outside
the raster keep B 50. Needs A 218 1 for the friction.

Volume budget (R 90 interval [s], <= 0 every step), rank 0:
REEF3D_SFLOW_Log/REEF3D-SFLOW-Hydro-Budget.dat, cumulative volumes [m3]
since the start: rain, sources (Q > 0), abstraction, evaporation,
infiltration, the storage sum(WL A) of the fluid cells and its change, the
residual dS - (rain + sources - abstraction - evaporation - infiltration),
which is the net inflow over the open boundaries plus the mass error of
the scheme (on a closed domain: the mass error), and the time integrals
of the inflow and outflow discharge of ioflow (Qi, Qo).

Checkpoints (P 48, I 45): the cumulative infiltration and the budget.
--------------------------------------------------------------------*/

class hydro_sflow : public increment
{
public:
    hydro_sflow(lexer*, fdm2D*, ghostcell*);
    virtual ~hydro_sflow() = default;

    static bool active(lexer *p) {return p->A10==2 && (p->R10>0 || p->R20>0.0 || p->R30>0 || p->R40>0 || p->R50>0);}
    static bool wetting(lexer *p) {return p->A10==2 && (p->R10>0 || p->R40>0);}

    // inputs, the roughness raster, the cells of the point sources, the storage at the start
    void ini(lexer*, fdm2D*, ghostcell*);

    // source rate S [m/s] of the step at p->simtime (the start of the step)
    void prepare(lexer*, fdm2D*, ghostcell*);

    // K += S, all fluid cells (continuity of a stage)
    void mass_source(lexer*, fdm2D*, slice &K);

    // after the continuity of the step (same dt): losses and the volumes of the step
    void sinks(lexer*, fdm2D*, ghostcell*);

    // after the step (p->simtime at its end): the budget output (the storage of sinks)
    void budget(lexer*, fdm2D*, ghostcell*);

    // checkpoints
    void ckpt_fields(const std::function<void(const char*, char, double*)>&, const std::function<void(const char*, char, int*)>&);
    void ckpt_state(std::vector<double>&) const;
    bool ckpt_restore(const std::vector<double>&);

    slice4 S;       // source rate of the step [m/s]: rain + point sources (Q > 0)
    slice4 R;       // rain rate of the step [m/s]
    slice4 F;       // cumulative infiltration [m]

    // budget, cumulative volumes [m3] since the start (the same on all ranks)
    double Vrain=0.0, Vsrc=0.0, Vabs=0.0, Vevap=0.0, Vinf=0.0, Vqi=0.0, Vqo=0.0, Vstore0=0.0, Vstore=0.0;

private:
    struct source
    {
        int ID = 0;
        double x = 0.0, y = 0.0, Q = 0.0;
        bool own = false;
        int ci = 0, cj = 0;
        double area = 0.0;
        std::vector<double> ht, hq;     // hydrograph_src_<ID>.dat
        double q = 0.0;                 // discharge of the step
    };

    void stop(lexer*, ghostcell*, const std::string&) const;
    double xc(lexer*, int, int) const;
    double yc(lexer*, int, int) const;
    double rain_series(lexer*, double t) const;
    void rain_rasters(lexer*, ghostcell*, double t);
    void load_record(lexer*, ghostcell*, int n, slice4 &dst);
    double storage(lexer*, fdm2D*, ghostcell*);
    static double ga_step(double F0, double P, double Ksdt);
    static bool read_series(const std::string &file, std::vector<double> &t, std::vector<double> &v);
    static double ipol(const std::vector<double> &t, const std::vector<double> &v, double tt, int mode);

    std::vector<source> src;

    std::vector<double> rt, rv;             // R 10 2
    std::vector<double> st;                 // R 10 3: times and files
    std::vector<std::string> sf;
    slice4 Ra, Rb;                          // R 10 3: records k, k+1 at the cells
    int rk = -2;                            // record in Ra (-2: none)

    double qrain = 0.0, qsrc = 0.0;         // global rates of the step [m3/s]
    double tnext = 0.0;
    std::ofstream bout;
};

#endif
