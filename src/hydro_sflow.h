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
evaporation, infiltration, roughness and soil rasters for SFLOW (A 10 2),
on Cartesian and curvilinear (G 50) grids, with mesh refinement (G 1-G 40)
and on multi-block grids.

One object per grid (fdm2D::hydro): level 0 holds the inputs of all grids
(the rain series and rasters, the point sources, the rasters of the
roughness and the soil), the budget and its output; the patches of the
mesh refinement and the blocks 1.. of a multi-block grid have their own
object for the fields of their cells (sflow_amr_hydro.cpp), which reads
the inputs of level 0.

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
  R 33 1    the infiltration parameters per cell from rasters (cells
            outside or at NODATA keep R 31 / R 32): R 30 1 inf_rate.asc
            [mm/h]; R 30 2 ga_ks.asc [mm/h], ga_psi.asc [m], ga_dtheta.asc [-]
            (each file optional)
  R 34      sealed fraction s of every cell [0..1] (default 0), R 35 1:
            from sealed.asc; the infiltration of a cell is (1 - s) times
            that of its soil (F: the cumulative infiltration of the soil)
  Shares: when the water is short, evaporation, infiltration and
  abstraction are reduced in proportion.

Wetting (R 10 or R 40): a dry cell gets wet by its own depth (above
A 244 + 1e-6 m) as well, not only from a higher wet neighbour (A 243 2),
and an isolated wet cell is not dried; see sflow_eta::wetdry. Use
A 243 2: with A 243 1 the rain on a dry cell below 1e-6 m per step is
reset with the film.

Roughness R 50 1: Manning n from manning.asc (ESRI ASCII, world
coordinates), ks = (20 n)^6 so that n = ks^(1/6)/20 of sflow_rough_manning
and the implicit friction of sflow_momentum_SI is that n; cells outside
the raster keep B 50. Needs A 218 1 for the friction.

Mesh refinement and multi-block grids: every grid applies the sources and
losses on its cells (on the cells a finer patch covers they are replaced
by the restriction); the budget counts the leaf cells (sflow_amr::rsi_leaf).
A fresh patch takes the cumulative infiltration of the old patch where it
had one, else of its parent.

Volume budget (R 90 interval [s], <= 0 every step), rank 0:
REEF3D_SFLOW_Log/REEF3D-SFLOW-Hydro-Budget.dat, cumulative volumes [m3]
since the start: rain, sources (Q > 0), abstraction, evaporation,
infiltration, the storage sum(WL A) of the fluid leaf cells and its change,
the residual dS - (rain + sources - abstraction - evaporation -
infiltration) = net inflow over the boundaries + mass error, the inflow
and outflow over the open boundaries (semi-implicit: the discharges of the
eta system, exact; explicit: the time integrals of ioflow Qi, Qo) and the
error = residual - (inflow - outflow) (semi-implicit: the mass error of
the scheme).

Coupled with the catchment model (R 100 2, hydro_couple.h): the source rates
C of the catchment (runoff of the cells under the grid, discharge of the
flow paths that enter it) are added to S; they count as sources in the
budget (Vcpl: their part). Rain R 10 on top of it would count twice.

Output (with any HYDRO key): REEF3D_SFLOW_VTP_FSF rain [mm/h], rain_depth
[mm] (cumulative), infiltration [mm] (cumulative).

Checkpoints (P 48, I 45): the cumulative rain and infiltration and the
budget.
--------------------------------------------------------------------*/

class hydro_sflow : public increment
{
public:
    // level 0 (lv0 nullptr) or another grid with the inputs of level 0
    hydro_sflow(lexer*, fdm2D*, ghostcell*, const hydro_sflow *lv0=nullptr);
    virtual ~hydro_sflow() = default;

    static bool active(lexer *p) {return p->A10==2 && (p->R10>0 || p->R20>0.0 || p->R30>0 || p->R40>0 || p->R50>0 || p->R100==2);}
    static bool wetting(lexer *p) {return p->A10==2 && (p->R10>0 || p->R40>0 || p->R100==2);}

    // level 0: the inputs (collective), its grid, the storage at the start, the budget file
    void ini(lexer*, fdm2D*, ghostcell*);

    // every grid: the rasters at its cells (roughness, soil) and the cells of the point sources (no communication)
    void grid_ini(lexer*, fdm2D*);

    // level 0 without patches or blocks: the step (prepare, sinks after the momentum step, budget after the step)
    void prepare(lexer*, fdm2D*, ghostcell*);
    void sinks(lexer*, fdm2D*, ghostcell*);
    void budget(lexer*, fdm2D*, ghostcell*);

    // the parts of the step for several grids (sflow_amr_hydro.cpp)
    void inputs_step(lexer*, ghostcell*);                  // level 0: rain records and discharges at p->simtime
    void grid_prepare(lexer*);                             // every grid: S and R of the step
    double grid_rain(lexer*, const std::function<bool(int,int)>&);       // sum R A of the cells
    void grid_sinks(lexer*, fdm2D*, double dt, const std::function<bool(int,int)>&, double *v);     // v += evap, inf, abs
    double grid_storage(lexer*, fdm2D*, const std::function<bool(int,int)>&);
    void step_volumes(lexer*, ghostcell*, double ve, double vi, double va, double store);    // level 0: the volumes of the step (qrain of step_rain or prepare)
    void boundary(double qin, double qout) {bq_in = qin; bq_out = qout; bq_set = true;}

    // K += S, all fluid cells (continuity of a stage)
    void mass_source(lexer*, fdm2D*, slice &K);

    // output (sflow_vtp_fsf)
    void offset_vtp(lexer*, int*, int&);
    void name_vtp(std::ostream&, int*, int&);
    static void name_pvtp(std::ostream&);
    void print_vtp(lexer*, std::ostream&);

    // checkpoints
    void ckpt_fields(const std::function<void(const char*, char, double*)>&, const std::function<void(const char*, char, int*)>&);
    void ckpt_state(std::vector<double>&) const;
    bool ckpt_restore(const std::vector<double>&);

    slice4 S;       // source rate of the step [m/s]: rain + point sources (Q > 0)
    slice4 R;       // rain rate of the step [m/s]
    slice4 F;       // cumulative infiltration into the soil (the pervious part of the cell) [m]
    slice4 Rc;      // cumulative rain [m]
    slice4 fc;      // R 30 1: infiltration rate [m/s]; R 30 2: Ks [m/s]
    slice4 Pc;      // R 30 2: psi dtheta [m]
    slice4 sc;      // sealed fraction
    slice4 L;       // the losses of the step [m] (the semi-implicit water column: sflow_momentum_SI::water_column)
    slice4 C;       // R 100 2: source rate from the catchment model [m/s] (hydro_couple), with cpl
    slice4 Bq;      // R 100 2: outflow over the open boundary faces of a cell since the last exchange [m3] (semi-implicit)
    slice4 Fp;      // R 100 2: (1 - sealed) F at the last exchange [m]
    bool cpl = false;
    double qcpl = 0.0, Vcpl = 0.0;  // R 100 2: its global rate [m3/s] (part of the sources) and cumulative volume [m3]

    // budget, cumulative volumes [m3] since the start (the same on all ranks)
    double Vrain=0.0, Vsrc=0.0, Vabs=0.0, Vevap=0.0, Vinf=0.0, Vqi=0.0, Vqo=0.0, Vstore0=0.0, Vstore=0.0;
    bool store0_leaf = false;   // mesh refinement: the storage at the start over the leaf cells of all grids (sflow_amr_hydro.cpp)
    void step_rain(double q) {qrain = q;}

private:
    struct source
    {
        int ID = 0;
        double x = 0.0, y = 0.0, Q = 0.0;
        std::vector<double> ht, hq;     // hydrograph_src_<ID>.dat
        double q = 0.0;                 // discharge of the step
    };
    struct cell
    {
        int n = 0;                      // the source
        int ci = 0, cj = 0;
        double area = 0.0;
    };

    void stop(lexer*, ghostcell*, const std::string&) const;
    double xc(lexer*, int, int) const;
    double yc(lexer*, int, int) const;
    void load_record(lexer*, ghostcell*, int n, hydro_raster &dst);
public:
    bool locate(lexer*, double x, double y, int &ci, int &cj);      // a point (model coordinates) in a fluid cell of this grid
private:
    static double ga_step(double F0, double P, double Ksdt);
    static bool read_series(const std::string &file, std::vector<double> &t, std::vector<double> &v);
    static double ipol(const std::vector<double> &t, const std::vector<double> &v, double tt, int mode);
    static double rate(const hydro_raster &r, double x, double y);

    const hydro_sflow *in;                  // the inputs: level 0 (this on level 0)
    std::vector<cell> cells;                // the point sources in cells of this grid (one grid of one rank on level 0)

    // inputs (level 0)
    std::vector<source> src;
    std::vector<double> rt, rv;             // R 10 2
    std::vector<double> st;                 // R 10 3: times and files
    std::vector<std::string> sf;
    hydro_raster ra, rb;                    // R 10 3: records k, k+1
    int rk = -2;                            // record in ra (-2: none)
    double rw = 0.0;                        // R 10 3: time weight of rb, R 10 1/2: the rate [m/s]
    bool rlin = false;
    hydro_raster rman, rinf, rks, rpsi, rdth, rseal;
    bool hman = false, hinf = false, hks = false, hpsi = false, hdth = false, hseal = false;

    double qrain = 0.0, qsrc = 0.0;         // global rates of the step [m3/s]
    double bq_in = 0.0, bq_out = 0.0;       // the discharges over the open boundaries of the step [m3/s]
    bool bq_set = false;
    double tnext = 0.0;
    std::ofstream bout;
};

#endif
