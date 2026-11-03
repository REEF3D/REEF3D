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

#ifndef HYDRO_CORE_H_
#define HYDRO_CORE_H_

#include"increment.h"
#include<fstream>
#include<string>
#include<vector>

class lexer;
class ghostcell;

using namespace std;

/*
REEF3D::HYDRO catchment model (R 100)

Grid:     hydro_grid.dat from DIVEMesh (W 10): hydrology cells in topological order (upstream first),
          each with area, elevation bands (equal area, glacier fraction), downstream cell, channel length and slope
Forcing:  hydro_forcing.dat: stations (x, y, z) with precipitation [mm/h], temperature [C] and optionally
          potential evaporation [mm/h], piecewise constant in time; inverse distance weights (power 2),
          precipitation and temperature corrected to the band elevation (PCALT, TCALT)
          PET without columns: Oudin et al. (2005) from temperature, latitude (R 104) and day of the year (R 105)
Vertical: HBV (Bergstroem 1976, formulation of HBV-light, Seibert & Vis 2012) per cell and band:
          degree-day snow with liquid water and refreezing, glacier melt on snow-free glacier, soil box (FC, LP, BETA),
          response per cell: upper zone (K0 above UZL, K1, percolation PERC), lower zone (K2)
          recession coefficients per day are converted to the time step: k = 1 - (1-K)^(dt/1 day)
Routing:  MAXBAS triangular weights per cell, then the cell network: R 106 0 the outflow of a cell reaches the
          downstream cell in the same step, R 106 1 kinematic wave on the channel of every cell
          (implicit, Chow et al. 1988, A = alpha Q^0.6, width W = CH_WA (upstream area/km2)^CH_WB)
Parallel: cells are split into contiguous ranges over the ranks (vertical processes are local); the
          lateral inflows are gathered and the network is routed on every rank (identical on all ranks)
Output:   REEF3D_HYDRO/REEF3D-HYDRO-Gauges.dat, -Budget.dat, -Means.dat, REEF3D-HYDRO-<n>.vtk (R 103)
Coupled:  R 100 2 with SFLOW (hydro_couple.h): the share of the runoff of the cells under the SFLOW grid and the
          discharge of the cells whose flow path enters it go to SFLOW, the outflow of SFLOW comes back into
          the cell below the SFLOW grid, the infiltration of SFLOW into the upper zone
Hot start: R 107 1 writes REEF3D_HYDRO/hydro_state.dat at the end, R 108 1 starts from hydro_state.dat
Reservoirs: R 110 1, hydro_reservoir.dat: a level pool with a weir and a release series in place of the channel of a
          cell (hydro_core_reservoir.cpp), REEF3D_HYDRO/REEF3D-HYDRO-Reservoirs.dat
*/

class hydro_core : public increment
{
public:
    hydro_core(lexer*, ghostcell*);
    ~hydro_core();

    static bool active(lexer*);

    void ini(lexer*, ghostcell*);
    void step(lexer*, ghostcell*);              // advance the model from t to t+dt
    void run(lexer*, ghostcell*);               // R 100 1: catchment model only
    void begin(lexer*, ghostcell*);             // outputs at the start (after ini)
    void end(lexer*, ghostcell*);               // outputs at the end, the state (R 107 1)
    void regression_fields(vector<string>&, vector<vector<double>>&) const;   // owned state (regression dump)
    vector<double> regression_sums() const;

    enum par {TT,CFMAX,SFCF,PCORR,CFR,CWH,CFGLAC,FC,LP,BETA,PERC,UZL,K0,K1,K2,MAXBAS,PCALT,TCALT,
              CH_N,CH_WA,CH_WB,INI_SWE,INI_SM,INI_SUZ,INI_SLZ,NPAR};
    static const char *parname[NPAR];

    double t, dt;                               // time and time step [s]
    int N, nb;                                  // cells and bands

    // hot start (R 107 1 writes, R 108 1 reads REEF3D_HYDRO/hydro_state.dat / hydro_state.dat)
    void write_state(lexer*, ghostcell*);
    void read_state(lexer*, ghostcell*);
    void gather_cells(lexer*, ghostcell*, vector<double>&);     // rank 0: the state of all cells
    void scatter_cells(const double*);                          // the owned cells from the state of all cells

    // checkpoints of a coupled run (P 48, I 45): everything on rank 0 (collective), restored on all ranks
    void ckpt_state(lexer*, ghostcell*, vector<double>&);
    bool ckpt_restore(lexer*, const vector<double>&);
    void flush();
    vector<double> Q, qgen;                     // outflow of every cell [m3/s], its runoff into the channel [m3/s]
    vector<double> len, alpha;                  // channel length [m] and A = alpha Q^0.6 of every cell

    // coupling with SFLOW (R 100 2, hydro_couple): set before step(), read after it
    vector<double> share;                       // fraction of the runoff of a cell that goes to SFLOW
    vector<char> enter;                         // the outflow of the cell enters SFLOW instead of the downstream cell
    vector<double> qret;                        // outflow of SFLOW into every cell over the step [m3/s]
    vector<double> addv;                        // water from SFLOW into the upper zone of every cell over the step [m3] (all ranks)
    vector<double> qsf;                         // runoff of a cell to SFLOW over the step [m3/s]
    vector<double> Qent;                        // discharge of a cell into SFLOW over the step [m3/s]
    double Vto, Vfrom;                          // cumulative volumes to and from SFLOW [m3] (the same on all ranks)
    double total_storage(ghostcell*) const;     // collective

    // grid (hydro_grid.dat)
    int ni, nj;
    double x0, y0, D;
    vector<int> I, J, down, sub, zone;
    vector<double> xc, yc, area, zmean, slope, aup, lakef;
    vector<double> bz, bf, bg;                  // per cell and band
    vector<int> gid, gcell;
    vector<int> poff;                           // version 2: the fine flow path of cell n, points poff[n] .. poff[n+1]-1,
    vector<double> px, py;                      // from its outlet pixel to the outlet pixel of the downstream cell

    // reservoirs (R 110 1, hydro_reservoir.dat): a level pool in place of the channel of a cell
    struct reservoir
    {
        int id, cell;
        double area, zb, crest, width, cw;      // area [m2], bottom, crest [m], weir width [m], weir coefficient [m^0.5/s]
        vector<double> rt, rq;                  // release series hydro_release_<id>.dat (t [s], Q [m3/s], linear)
        double h, qin, qrel, qspill;            // level [m], inflow, release, spill of the last step [m3/s]
    };
    vector<reservoir> resv;
    vector<int> rescell;                        // reservoir of every cell, -1 none
    double reservoir_storage() const;

    // volumes [m3], cumulative (Vp, Vea, Vgm over the owned cells, Vout the same on all ranks)
    double Vp, Vea, Vgm, Vout, S0;

private:
    void stop(lexer*, ghostcell*, const string&) const;
    void read_grid(lexer*, ghostcell*);
    void read_param(lexer*, ghostcell*);
    void read_forcing(lexer*, ghostcell*);
    void forcing_mean(lexer*, double t0, double t1);
    double pet_oudin(lexer*, double T, double tday) const;
    void vertical(lexer*);
    void route(lexer*, ghostcell*);
    double kinwave(double c, double alpha, double R) const;
    void read_reservoirs(lexer*, ghostcell*);
    double reservoir_step(lexer*, reservoir&, double qin);     // the outflow of the step [m3/s]
    void reservoir_output(lexer*);
    ofstream fr;
    void gather(ghostcell*, vector<double>&) const;
    void budget(lexer*, ghostcell*, bool force);
    void print_vtk(lexer*, ghostcell*);
    double owned_storage() const;
    double P(int n, int k) const {return par[n*NPAR+k];}

    // partition: rank r owns cells [displ[r], displ[r]+count[r])
    int n0, n1;
    vector<int> cnt, displ;

    vector<double> par;                         // per cell

    // forcing
    vector<double> sx, sy, sz;                  // stations
    vector<double> ft;                          // times of the rows
    vector<vector<double>> fP, fT, fE;          // per row, per station
    bool fhasE;
    vector<double> mP, mT, mE;                  // mean over the step, per station
    vector<double> wgt;                         // inverse distance weights per owned cell and station

    // states, per owned cell (and band)
    vector<double> SP, WC, SM;                  // per band: dry snow, liquid water in snow, soil moisture [mm]
    vector<double> SUZ, SLZ;                    // per cell [mm]
    vector<double> buf;                         // MAXBAS buffer per cell [mm], nmb values
    vector<double> mbw;                         // MAXBAS weights per cell
    int nmb;

    // fluxes of the step per owned cell
    vector<double> cP, cT, cE, cEa, cRun;       // cell means for output [mm/h], [C]

    // routing (all cells, every rank)
    vector<double> Qold, Qin;

    double tout, tvtk, tlast;               // next output, next grid output, last output
    int nvtk;
    bool kw;                                    // kinematic wave (R 106 1)
    ofstream fg, fb, fm;
    bool first, opened;
};

#endif
