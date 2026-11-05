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

#ifndef MET_BLEND_H_
#define MET_BLEND_H_

#include"met_source.h"
#include<string>

/*--------------------------------------------------------------------
REEF3D atmospheric forcing W 200 4: composite of a background met file,
a moving nest and a parametric vortex. No lexer or MPI dependency (unit
test Regression/unit/met_test.cpp against the reference composite of
tools/met_forcing.py blend).

Layers, from the bottom:
  background  met-forcing.r3m, e.g. ERA5, BARRA-C2, NORA3 (met_forcing.py grid)
  nest        met-nest.r3m (W 220 1 x), e.g. a storm-following HAFS nest
              (met_forcing.py nest, a moving grid); weight
                wn = taper(d/W 221 width) taper_t(t, first, last record, W 221 T)
              d: distance of the point inside the edge of the nest grid at
              time t (its origin linear between the records); without pa in
              the nest file the pressure stays that of the background
  vortex      met-vortex.trk (W 220 x 1), the model W 215 with W 216-219; weight
                wv = taper((r2 - r)/(r2 - r1)) taper_t(t, first, last record, W 224)
              r: distance from the centre on the sphere, 1 inside r1, 0 beyond r2;
              W 222 basis r1 r2: 1 multiples of the rmw of the track, 2 multiples
              of the mean 34 kt radius, 3 metres
  Each layer: f = f_below + w (f_layer - f_below), u10 and v10 along the model
  axes and pa; w = 0 and w = 1 give the layers bitwise.

taper (W 225): 1 0.5 (1 - cos(pi s)) (continuous slope at both ends),
2 s (linear, as ADCIRC NWS 29/30); taper_t: in time, the same over T s
after the first and before the last record, 0 outside (T = 0: a step).

The ring: 16 points on the circle of radius r2 about the centre, at the
background record times (read directly from the file, the same on every
rank), linear in time between them; there the background (pb, ub) and the
vortex (pv, uv) are compared.

Vortex pressure W 223 1: the deficit of the vortex is rescaled so that the
vortex merges into the background on the ring (the ring means agree),
  pv' = pc + (pv - pc) max((<pb> - pc)/(<pv> - pc), 0)
The wind of the vortex is not changed. W 223 0: pv as fitted (pn of the
track).

Nest in time W 227 1: storm-relative, the two nest records around t are
shifted with the storm centre c (from the track) before the linear
interpolation,
  fn(t,x) = (1-w) f_k(x - (c(t) - c_k)) + w f_k+1(x - (c(t) - c_k+1))
so that the core of a moving storm is carried, not cross-faded between
the records (a storm at 20 km/h moves half its Rmax between hourly
records); the nest window of the rank grows by the largest move between
two records. Needs met-vortex.trk, also with the vortex layer off.
W 227 0: linear in time at fixed points (as met_file).

Environment W 226 1: the vortex carries the large-scale flow of the
background, i.e. the background minus the vortex on the ring: the mean
<ub - uv> is added to the wind of the vortex, the gradient of pb - pv
(least-squares plane through the ring) to its pressure,
  uv' = uv + <ub - uv>,  pv' = pv + grad(pb - pv) . (x - x_centre)
Subtracting the vortex keeps its own asymmetry (translation, quadrants)
from being counted twice. W 226 0: the vortex alone inside r1 (as ADCIRC
NWS 29/30).
--------------------------------------------------------------------*/

// the background on the ring r2 at a background record time
struct met_blend_ring
{
    double p = 0.0, pv = 0.0;           // mean pressure of the background and of the vortex
    double u = 0.0, v = 0.0;            // mean wind of the background minus the vortex (W 226 1)
    double gx = 0.0, gy = 0.0;          // gradient of the background pressure minus the vortex pressure
};

struct met_blend_par
{
    int nest = 0, vortex = 1;                           // W 220
    double nest_width = 50000.0, nest_time = 3600.0;    // W 221
    int basis = 1;                                      // W 222
    double r1 = 3.0, r2 = 6.0;
    int pressure = 1;                                   // W 223
    double vortex_time = 21600.0;                       // W 224
    int shape = 1;                                      // W 225
    int env = 0;                                        // W 226
    int nest_relative = 0;                              // W 227
};

class met_source_blend final : public met_source
{
public:
    met_source_blend(const met_frame &fr, const met_blend_par &p) : frame(fr), par(p) {}

    // background file, nest file and track (the names of unused layers are not opened)
    bool open(const std::string &bgfile, const std::string &nestfile, const std::string &trkfile, std::string &error);

    bool advance(double t, std::string &error) override;
    void at(double t, double x, double y, double &u, double &v, double &pa) const override;
    void window(double xs, double xe, double ys, double ye) override;
    double first_time() const override {return bg.time(0);}

    static double taper(double s, int shape);
    static double taper_time(double t, double t0, double t1, double T, int shape);

    // blend radii [m] at time t
    void radii(double t, double &r1m, double &r2m) const;

    // the background on the ring of radius r2 about the centre at the time of background record n (16 points)
    bool ring(int n, met_blend_ring &rg, std::string &error);

    met_frame frame;
    met_blend_par par;
    met_file bg, nest;
    met_vortex vx;

    // state of the last advance
    double tnow = 0.0;
    double wn_t = 0.0, wv_t = 0.0;                      // weights in time
    double nest_x0 = 0.0, nest_y0 = 0.0, nest_lx = 0.0, nest_ly = 0.0;
    double nest_w = 0.0, nest_d[2][2] = {{0.0,0.0},{0.0,0.0}}, nest_margin = 0.0;   // W 227 1: weight, shifts, window margin
    double lonc = 0.0, latc = 0.0, xcm = 0.0, ycm = 0.0, pc = 0.0, pn = 0.0, rmw = 0.0, r34 = 0.0, r1m = 0.0, r2m = 0.0;
    met_blend_ring env;                                 // the ring at time t (p = pv = pn without W 223 1, W 226 1)
    long long ring_samples = 0;                         // background records sampled on the ring (log, unit test)

private:
    int ring_n[2] = {-1, -1};
    met_blend_ring ring_s[2];
};

#endif
