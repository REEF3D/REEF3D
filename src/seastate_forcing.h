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


#ifndef SEASTATE_FORCING_H_
#define SEASTATE_FORCING_H_

#include<fstream>
#include<string>
#include<vector>

class seastate_grid;

/*--------------------------------------------------------------------
REEF3D::SEASTATE external forcing - file readers (Phase 5a)

No lexer or MPI dependency (unit test Regression/unit/seastate_test.cpp).
Both readers stream through their file: only the two records around
the model time are kept, so long hindcasts need little memory. Times
must increase from record to record.

seastate_datetime(v)   SWAN date-time YYYYMMDD.HHMMSS -> seconds since
                       1970-01-01 00:00:00 (proleptic Gregorian, UTC)

seastate_spc_series    nonstationary SWAN 2D spectral file
                       (seastate-boundary.spc, A 711 3): many locations
                       (LOCATIONS: x y in model coordinates; LONLAT is
                       rejected, convert with tools/seastate_forcing.py),
                       AFREQ/RFREQ, CDIR/NDIR, QUANT VaDens/EnDens, then
                       for every time (TIME, line YYYYMMDD.HHMMSS) one
                       FACTOR/ZERO/NODATA block per location; also a
                       stationary file (no TIME: one record at time 0).
                       Spectra are stored on the model's spectral grid
                       (action density, seastate_swan_spc::to_grid).
  open      header; the first two records
  advance   read records until rec[1] is the first record at or after t
  N(k,loc)  spectrum of record k = 0, 1, location loc on the model grid
  weight(t) linear time weight of rec[1] (clamped to [0,1])

The wind field file seastate-wind.dat (A 730 2) is retired (storm surge
P6): the wind of SEASTATE is that of the REEF3D met forcing W 200
(met_provider; met-forcing.r3m, tools/met_forcing.py, which converts an
old seastate-wind.dat with 'met_forcing.py seastate').
--------------------------------------------------------------------*/

double seastate_datetime(double yyyymmdd_hhmmss);

class seastate_spc_series
{
public:
    bool open(const std::string &file, const seastate_grid &g, std::string &error);
    bool advance(double t, std::string &error);       // t: seconds, same clock as time()
    double weight(double t) const;

    const std::vector<float> &N(int k, int loc) const {return rec[k][loc];}
    double time(int k) const {return rtime[k];}
    bool stationary() const {return !timed;}

    std::vector<double> xs, ys;                      // locations
    int nloc = 0, records = 0;
    std::vector<double> f, dir;                      // file frequencies [Hz], Cartesian directions [deg]

private:
    bool read_record(std::vector<std::vector<float>> &out, double &t, std::string &error);

    const seastate_grid *g = nullptr;
    std::ifstream in;
    std::string file;
    bool timed = false, nautical = false, energy = false, eof = false;
    double exception = -99.0;
    std::vector<size_t> order;                        // sorted direction -> file direction
    std::vector<std::vector<float>> rec[2];
    double rtime[2] = {0.0, 0.0};
    std::string pending;                              // first key of the next record
};

#endif
