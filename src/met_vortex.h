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

#ifndef MET_VORTEX_H_
#define MET_VORTEX_H_

#include<string>
#include<vector>

/*--------------------------------------------------------------------
REEF3D parametric tropical cyclone vortex (W 200 3). No lexer or MPI
dependency (unit test Regression/unit/met_test.cpp against the tables
of tools/met_track.py, which fits the vortex from a best track and
holds the reference implementation of every formula here).

Track file met-vortex.trk (ASCII, SI units; format in met_track.py):
  REEF3D-MET-TRACK 1 nrec rho_a, '#' comment lines, then per record
  YYYYMMDD.HHMMSS lon lat ut vt pc pn vmaxg rmw B80 vms bs xn rn Bg0 phi0
  and for the quadrants NE SE SW NW and the isotachs 34 50 64 kt:
  R Rm Bg phi (R = 0: no isotach).

Models (W 215):
  1 Holland (1980), symmetric:
      Vg = sqrt(B dp/rho_a (Rm/r)^B e^-(Rm/r)^B + (r f/2)^2) - r f/2
  2 GAHM (Gao 2018, ADCIRC NWS 20), asymmetric: Rm, Bg, phi of the
      quadrant, linear in r between its isotachs (outside the outermost
      the outermost, inside the innermost the innermost), inverse-square
      in azimuth between the two quadrants (the nearest one within 1 deg)
      Vg = sqrt(Vmax^2 (1+1/Ro) e^(phi (1-(Rm/r)^Bg)) (Rm/r)^Bg + (r f/2)^2) - r f/2
      Ro = Vmax/(|f| Rm), p = pc + dp e^(-phi (Rm/r)^Bg)
  3 Holland, Belanger and Fritz (2010): surface wind
      V = Vms [(Rm/r)^bs e^(1-(Rm/r)^bs)]^x, x = 0.5 inside Rm, linear
      to xn at rn outside
f = 2 Omega sin(lat) at the centre; in the southern hemisphere the
vortex turns clockwise. Inside 1 km: no wind, pc.

Surface wind (east, north):
  speed  W 216 Wrf Vg (models 1, 2), V (model 3)
  angle  (90 - bearing) + s (90 + beta), s = sign(f), beta the inflow
         angle of ADCIRC (W 218 1): 10 r/Rm inside Rm,
         10 + 75 (r/Rm - 1) up to 1.2 Rm, 25 deg outside
  plus   the translation (ut, vt) scaled by Vg/Vmax (W 219 1, ADCIRC),
         in full (W 219 2) or not (0)
  times  the 1-min to 10-min factor W 217 (0.8928)

In time: the centre and the translation linear between the records,
the fields of the two records around the time (both about that centre)
linear in time; before the first and after the last record the
first and the last record.
--------------------------------------------------------------------*/

struct met_vortex_rec
{
    double t = 0.0;
    double lon = 0.0, lat = 0.0, ut = 0.0, vt = 0.0, pc = 0.0, pn = 0.0;
    double vmaxg = 0.0, rmw = 0.0, B80 = 0.0, vms = 0.0, bs = 0.0, xn = 0.5, rn = 0.0, Bg0 = 0.0, phi0 = 1.0;
    double R[4][3], Rm[4][3], Bg[4][3], phi[4][3];
};

struct met_vortex_par
{
    int model = 2;              // W 215
    double wrf = 0.9;           // W 216
    double one2ten = 0.8928;    // W 217
    int inflow = 1;             // W 218
    int trans = 1;              // W 219
};

class met_vortex
{
public:
    bool open(const std::string &file, std::string &error);

    // wind (east, north) [m/s] and pressure [Pa] at (lon, lat) [deg] and time t [s since 1970]
    void at(double t, double lon, double lat, double &ue, double &vn, double &pa) const;

    // one record about the centre (lonc, latc) with the translation (ut, vt)
    void eval(const met_vortex_rec &r, double lonc, double latc, double ut, double vt,
              double lon, double lat, double &ue, double &vn, double &pa) const;

    // centre at time t (lon, lat [deg])
    void centre(double t, double &lon, double &lat) const;

    // centre, central and outer pressure [Pa], rmw and the mean 34 kt radius [m] at time t, linear between
    // the records (blend W 200 4); r34 of a record without 34 kt isotachs: 3 rmw
    void state(double t, double &lon, double &lat, double &pc, double &pn, double &rmw, double &r34) const;
    static double r34_mean(const met_vortex_rec &r);

    // point at distance d [m] and bearing brg [deg clockwise from north] from (lon, lat) on the sphere
    static void destination(double lon, double lat, double brg, double d, double &lon2, double &lat2);

    // first and last time of the track
    double t_first() const {return rec.front().t;}
    double t_last() const {return rec.back().t;}

    // profiles (public for the unit test)
    static void dist_bearing(double lon1, double lat1, double lon2, double lat2, double &d, double &brg);
    static double fang(double r, double rm);
    static double holland80(double r, double rm, double B, double dp, double rho, double f);
    static double gahm(double r, double vmaxg, double rm, double Bg, double phi, double f);
    static double holland10(double r, double vms, double rm, double bs, double xn, double rn);

    std::vector<met_vortex_rec> rec;
    double rho = 1.15;
    met_vortex_par par;

private:
    void quad_params(const met_vortex_rec &r, int q, double rr, double &rm, double &bg, double &ph) const;
    void azimuth_params(const met_vortex_rec &r, double theta, double rr, double &rm, double &bg, double &ph) const;
    int bracket(double t, double &w) const;
};

#endif
