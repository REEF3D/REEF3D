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

#ifndef MET_FILE_H_
#define MET_FILE_H_

#include<fstream>
#include<string>
#include<vector>

/*--------------------------------------------------------------------
REEF3D atmospheric forcing (W 200 2): reader of the binary met file
met-forcing.r3m. No lexer or MPI dependency (unit test
Regression/unit/met_test.cpp); written by tools/met_forcing.py.

File format (version 1, little-endian, no padding)
  header, 144 bytes
    char[8]  "REEF3DMT"
    int32    version (1), nx, ny, nvar (2: u10 v10, 3: u10 v10 pa),
             moving (0: one grid; 1: every record has its own origin,
             e.g. a storm-following nest), nrec, frame (0 none,
             1 transverse Mercator W 205, 2 projected CRS), epsg
    float64  x0 y0 dx dy          node (0,0) and spacing in model
                                  coordinates [m], axes along the model
                                  x and y axes
             lon0 lat0 rot X0 Y0  frame the file was written in (frame 1:
                                  W 205, W 206, checked against ctrl.txt;
                                  frame 2: epsg, rot and the offset X0 Y0
                                  subtracted from the projected coordinates)
             avg                  averaging period of u10 [s], 0: unknown
             reserved[3]
  nrec records, 32 + nvar*nx*ny*4 bytes each
    float64  time [s since 1970-01-01 00:00:00 UTC], xr0, yr0 (origin of
             the record's grid; = x0, y0 unless moving), reserved
    float32  u10[ny][nx], v10[ny][nx], (pa[ny][nx])
             wind 10 m above the sea [m/s] along the model x and y axes,
             air pressure at mean sea level [Pa]; i fastest
  Records in increasing time. nx = ny = 1: uniform time series.

Reading: every rank reads only the part of each record its grids need
(window, nodes around the box, two in every direction), and only the
two records around the model time are kept, so files of long hindcasts
on fine grids need little memory and no communication.

  open      header and the record times
  window    box in model coordinates the rank needs (grows the window;
            the next advance reads the records again)
  advance   records k, k+1 with time(k) <= t < time(k+1), clamped at
            the first and the last record
  at        bilinear in space (clamped to the grid), linear in time
            (clamped to the first and the last record)
  origin    origin of the grid at time t (moving grids: linear between the
            two loaded records)
  sample    one variable of record n at a point, read directly from the
            file (four nodes, no window): the same value on every rank

met_datetime(v)   YYYYMMDD.HHMMSS -> s since 1970-01-01 00:00:00 UTC
met_date(t)       s since 1970 -> "YYYY-MM-DD HH:MM:SS"
--------------------------------------------------------------------*/

double met_datetime(double yyyymmdd_hhmmss);
std::string met_date(double t);

class met_file
{
public:
    bool open(const std::string &file, std::string &error);
    void window(double xs, double xe, double ys, double ye);
    bool advance(double t, std::string &error);
    void at(double t, double x, double y, double &u, double &v, double &pa) const;

    void origin(double t, double &ox, double &oy) const;
    bool sample(int n, int var, double x, double y, double &val, std::string &error);

    // the two loaded records (q = 0, 1): time, origin, values at a point (bilinear, clamped), time weight
    // (storm-relative time interpolation of the blend W 227 1)
    double rec_time(int q) const {return rec[q].t;}
    void rec_origin(int q, double &ox, double &oy) const {ox = rec[q].xr0;  oy = rec[q].yr0;}
    void at_record(int q, double x, double y, double &u, double &v, double &pa) const {space(rec[q],x,y,u,v,pa);}
    double time_weight(double t) const {return weight(t);}

    double time(int n) const {return times[n];}

    int nx = 0, ny = 0, nvar = 0, moving = 0, nrec = 0, frame = 0, epsg = 0;
    double x0 = 0.0, y0 = 0.0, dx = 1.0, dy = 1.0;
    double lon0 = 0.0, lat0 = 0.0, rot = 0.0, X0 = 0.0, Y0 = 0.0, avg = 0.0;

    long long values_read = 0;      // floats read so far (log, unit test)

    static const int header_bytes = 144;
    static const int record_head_bytes = 32;

private:
    struct record
    {
        int n = -1;                 // record number, -1: not loaded
        double t = 0.0, xr0 = 0.0, yr0 = 0.0;
        int i0 = 0, i1 = -1, j0 = 0, j1 = -1;   // window of nodes, inclusive
        std::vector<float> u, v, pa;
    };

    bool read_record(int n, record &r, std::string &error);
    void space(const record &r, double x, double y, double &u, double &v, double &pa) const;
    double weight(double t) const;
    long long offset(int n) const;

    std::ifstream in;
    std::string file;
    std::vector<double> times, xr, yr;
    record rec[2];
    bool haswin = false;
    double wxs = 0.0, wxe = 0.0, wys = 0.0, wye = 0.0;
    bool swap = false;
};

#endif
