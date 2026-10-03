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

#ifndef WAVE_FIELD_H_
#define WAVE_FIELD_H_

#include"wave_source.h"
#include<vector>

class lexer;
class ghostcell;
class wave_lib;

/*--------------------------------------------------------------------
wave_field: superposition of wave sources (iowave redesign, phase 1).

Source 1 is the legacy wave (B 90, B 92, B 93, ...). It stays owned and
evaluated by wave_interface exactly as before, so a case without the
inputs below runs the unchanged code path. wave_field holds the
additional sources and returns their sum, which wave_interface adds.

Input (ctrl.txt, repeatable, read here on rank 0 and broadcast):
  B 500 id type H T                 source id (>= 2), B 92 type, H (Hs), T (Tp)
  B 501 id dir phase ts te t_ramp   direction and phase [deg], time window and ramp [s]
  B 502 id x0 y0                    origin in the legacy generation frame [m]
  B 504 id seed                     random phases of an irregular source
Spectrum shape, spreading and depth are taken from the B 8x / B 13x / B 94
input of the legacy source.

Rules (checked at setup): at most one nonlinear source in total; no HDC
or wavemaker (paddle) source beyond the legacy one; no decomposed
precalc (B 89 1) with more than one source until phase 2.
--------------------------------------------------------------------*/

class wave_field
{
public:
    wave_field(lexer*, ghostcell*);
    ~wave_field();

    int size() const {return (int)src.size();}

    double eta(lexer*, double, double);
    double u(lexer*, double, double, double);
    double v(lexer*, double, double, double);
    double w(lexer*, double, double, double);
    double fi(lexer*, double, double, double);

    void cache_points(lexer*, const std::vector<double>&, const std::vector<double>&);
    double eta_c(lexer*, int);
    double fi_c(lexer*, int, double);
    void uvw_c(lexer*, int, double, double&, double&, double&);

    void prestep(lexer*, ghostcell*);

    static bool nonlinear(int);
    bool exists(int) const;
    
    // sources used by the cached evaluation (eta_c, fi_c, uvw_c); nullptr: all
    const std::vector<int> *filter = nullptr;

private:
    void read(lexer*, ghostcell*);
    void check(lexer*);
    void build(lexer*, ghostcell*);
    void log(lexer*);

    // swaps a source's lexer context in for the scope of one evaluation
    struct scope
    {
        scope(lexer*, wave_source&);
        ~scope();
        lexer *p;
        wave_source &s;
        wave_lexer_context keep;
        double wavetime;
    };

    bool use(const wave_source*) const;
    
    std::vector<wave_source*> src;
};

#endif
