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

#ifndef SEASTATE_F_H_
#define SEASTATE_F_H_

#include"seastate.h"
#include"increment.h"
#include<fstream>
#include<vector>

class fdm_seastate;
class seastate_vtp;
class seastate_exchange;
class seastate_implicit;
class seastate_store;
class regression_dump;

using namespace std;

/*--------------------------------------------------------------------
REEF3D::SEASTATE - spectral (phase-averaged) wave model, A 10 7

2D horizontal grid of SFLOW (DIVEMesh), spectral grid (A 701-703),
block-sparse float32 action density (A 704), active cells from the
bathymetry (A 705), parametric initial spectrum (A 710), integrated wave
parameters, VTP output (P 181/P 182), integral log, regression dump.

Phase 1: implicit transport in x, y, sigma, theta (seastate_implicit),
nonstationary (A 700 1, time step A 706) or stationary (A 700 2),
iterations A 707, stationary convergence A 708, boundary spectrum
A 711 on the sides A 712, refraction A 713, frequency shift A 714,
halo exchange of the spectra (seastate_exchange). Source terms follow
in Phase 2.

  start   stand-alone run: ini, then the time loop calling step
  ini     set-up (environment, storage, initial spectrum)
  step    one time step of length A 706, so that a host model
          (SFLOW, NHFLOW) can drive the module later
--------------------------------------------------------------------*/

class seastate_f final : public seastate, public increment
{
public:
    seastate_f(lexer*, ghostcell*);
    virtual ~seastate_f();

    void start(lexer*, ghostcell*) override;

    void ini(lexer*, ghostcell*);
    void step(lexer*, ghostcell*);

    fdm_seastate *e;

private:
    void check_keys(lexer*, ghostcell*);
    void environment(lexer*, ghostcell*);
    void storage(lexer*, ghostcell*);
    void initial(lexer*, ghostcell*);
    void initial_parametric(lexer*, ghostcell*);
    void parametric_spectrum(lexer*, ghostcell*, std::vector<float>&, const char*);
    void boundary(lexer*, ghostcell*);
    void kinematics(lexer*, ghostcell*, double);
    void transport(lexer*, ghostcell*);
    void parameters(lexer*, ghostcell*);

    void log_ini(lexer*);
    void log_step(lexer*);

    seastate_vtp *pprint;
    regression_dump *preg;
    seastate_exchange *pex;
    seastate_implicit *psolv;
    seastate_store *N0;

    std::vector<float> Nb;      // boundary spectrum
    int side[4];                // sides with the boundary spectrum: x-, x+, y-, y+
    int iter_max, iter_done;
    double conv;
    std::vector<double> hs_old;

    ofstream integral;
    double etot, hsmax, hsmean, nmin, cells_active;
    double starttime, endtime;
};

#endif
