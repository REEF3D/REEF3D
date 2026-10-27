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

#ifndef IOFLOW_AMR_BLOCK_H_
#define IOFLOW_AMR_BLOCK_H_

#include"ioflow_void.h"
#include"bc_zone.h"
#include"background_state.h"
#include<string>
#include<vector>

// in- and outflow of a block 1.. of a multi-block grid (G 50 1, its level-0 grid of REEFAMR) through its open
// faces (G 54), as ioflow_f on the native grid with B 60 1: the discharge Q through the -i face (the velocity
// Q / wet area along the face normal in the ghost cells; NHFLOW: B 61 1 plain, else the log profile scaled
// to Q) and the water level at the +i face (start ramp F 64 from F 60).  G 54 types 3 and 4 take Q and the
// level from a time series file (discharge-<id>.dat, waterlevel-<id>.dat: lines "t value", linear in time,
// held beyond the ends), as the native hydrograph (B 60 3, 4).  A block lies on one rank: the sums over the
// face are local.  NHFLOW: an edge zone on the -i or +i face (B 520 methods 3-6 with B 522 b: Riemann, Flather,
// clamped level or discharge edge with its background B 523, B 510-515) sets the ghost cells as the edges of
// iowave on the native grid (iowave_nhflow_tide.cpp; a Riemann edge with the waves of its B 524 sources, the
// non-hydrostatic pressure of the incoming waves at the -i edge).  A relaxation zone (method 1) or a beach zone
// (method 2) on the -i or +i face reaches its B 521 width into the block, measured from the face along the grid
// lines, and relaxes as the zones of iowave (iowave_nhflow_relax.cpp): a relaxation zone towards its background
// and, with B 90 1 and B 98 2 on level 0, the waves of its B 524 sources (all sources without B 524); a beach
// zone towards its background (still water without B 523).  Everything else as ioflow_v.
class wave_interface;

class ioflow_amr_block : public ioflow_v
{
public:
    ioflow_amr_block(lexer*,ghostcell*,patchBC_interface*);
    virtual ~ioflow_amr_block();

    int in = 0;             // the -i face is open: discharge Q (m3/s)
    double Q = 0.0;
    int out = 0;            // the +i face is open: water level wsf (m)
    double wsf = 0.0;

    // the open faces of the G 54 lines nin (-i face) and nout (+i face), -1 none: types 1, 2 the constant
    // values, 3, 4 the time series (files checked by reefamr::blocks_open)
    void open_faces(lexer*,int,int);

    // the edge zones of the -i and +i faces (B 522), nullptr none, and the backgrounds (NHFLOW)
    void open_edges(lexer*,const bc_zone*,const bc_zone*,const background_state&);

    // the relaxation or beach zones of the -i and +i faces (B 522 methods 1, 2), nullptr none, the backgrounds and
    // the waves of the relaxation zones (B 90 1 and B 98 2 on level 0); NHFLOW
    // (the cells of the block: i0 .. i0+ni-1, j0 .. j0+nj-1 of its grid)
    void open_zones(lexer*,const bc_zone*,const bc_zone*,const background_state&,bool,int,int,int,int);
    bool has_zones() const { return !rzone.empty(); }
    void zones_ini(lexer*,fdm_nhf*,ghostcell*);

    // NHFLOW relaxation zones (B 95: the targets at the stage time t, the weights to the power dt / dt_ref)
    void wavegen_stage_nhflow(lexer*,fdm_nhf*,ghostcell*,double) override;
    void U_relax(lexer*,ghostcell*,double*,double*) override;
    void V_relax(lexer*,ghostcell*,double*,double*) override;
    void W_relax(lexer*,ghostcell*,double*,double*) override;
    void P_relax(lexer*,ghostcell*,double*) override;
    void WL_relax(lexer*,ghostcell*,slice&,slice&) override;

    // G 54 types 3, 4: the file of a time series and its value at time t
    static std::string series_file(int,int);
    static bool series_read(const std::string&,std::vector<double>&,std::vector<double>&);

    // SFLOW
    void discharge2D(lexer*,fdm2D*,ghostcell*) override;
    void waterlevel2D(lexer*,fdm2D*,ghostcell*,slice&) override;
    void inflow2D(lexer*,fdm2D*,ghostcell*,slice&,slice&,slice&,slice&) override;

    // NHFLOW
    void discharge_nhflow(lexer*,fdm_nhf*,ghostcell*) override;
    void inflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&) override;
    void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*,slice&) override;
    void rkinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*) override;
    void fsfinflow_nhflow(lexer*,fdm_nhf*,ghostcell*,slice&) override;

private:
    double level(lexer*) const;
    void series_update(lexer*);
    static double series_value(const std::vector<double>&,const std::vector<double>&,double);

    std::vector<double> sq_t, sq_v;     // discharge time series (G 54 type 3)
    std::vector<double> sw_t, sw_v;     // water level time series (G 54 type 4)

    // edge zones (B 522): ez[0] the -i face, ez[1] the +i face
    void edges_nhflow(lexer*,fdm_nhf*,ghostcell*,double*,double*,double*,double*,double*,double*);
    void edges_rk(lexer*,fdm_nhf*,double*,double*,double*,double*,double*,double*);
    void edges_wl(lexer*,fdm_nhf*,slice&);
    double col_ubar(lexer*,double*);
    double bg_prof(lexer*,int,int,double);

    std::vector<bc_zone> ezone;         // the zones of the faces
    int ez[2] = {-1,-1};                // their index per face, -1 none
    background_state bgs;
    std::vector<double> edge_h;         // the ghost depth of the Riemann / clamped level edges of the step
    void qin_nhflow(lexer*,fdm_nhf*,double&,double&);
    void inflow_q(lexer*,fdm_nhf*,double*,double*,double*,double*,double*,double*);
    void curv_dir(lexer*,double*,double*,double*,double*);

    double Ui = 0.0;        // inflow velocity normal to the -i face (m/s)
    double Ai = 0.0;        // wet area of the -i face (m2)

    // the waves of the relaxation zones and of the Riemann edges with sources (B 92, B 500-504).  The wave
    // members of the lexer are RUN members of the refinement (level 0 sets them, lexer_amr.cpp): the block keeps
    // the values its waves set (wown) and swaps them in only while it evaluates the waves (wave_in, wave_out)
    wave_interface *pwave = nullptr;
    bool wave_on(lexer*,ghostcell*);
    static const int NWAVE = 22;
    struct wave_run { double v[NWAVE]; int n; };
    wave_run wown;
    static void wave_get(const lexer*,wave_run&);
    static void wave_set(lexer*,const wave_run&);
    void wave_in(lexer*,wave_run&) const;
    void wave_out(lexer*,const wave_run&) const;
    double wramp(lexer*) const;
    void wave_frame(lexer*,double,double,double&,double&) const;
    std::vector<double> ew_eta, ew_u, ew_v;    // Riemann edge with waves: per column eta and the column means of u, v
    void edge_pressure(lexer*,fdm_nhf*,ghostcell*,double,double);

    // relaxation and beach zones (B 522 methods 1, 2): rz[0] the -i face, rz[1] the +i face
    std::vector<bc_zone> rzone;
    int rz[2] = {-1,-1};
    bool gen_waves = false;                     // B 90 1 and B 98 2 on level 0: the relaxation zones generate waves
    std::vector<int> zcol;                      // the columns (slice index) in a zone
    std::vector<double> zgen, zbch;             // their weights r (1: no relaxation), with B 95 to the power dt / dt_ref
    std::vector<double> zgen0, zbch0;           // ... as built
    std::vector<int> zgz, zbz;                  // their zone (index in rzone), -1 none
    double zfac = 1.0;
    int zi0 = 0, zni = 0, zj0 = 0, znj = 0;     // the cells of the block in its grid
    double tstage = -1.0e20;                    // B 95: the time of the targets
    bool zbuilt = false;
    bool zgen_only = false;                     // zones_ini: the water level of the relaxation zones only
    bool ztold = false;                         // a zone longer than the block reported
    void zones_build(lexer*);
    double zone_time(lexer*) const;
    void zone_stage(lexer*);
    template<class F> void zone_columns(lexer*,F&&);
};

#endif
