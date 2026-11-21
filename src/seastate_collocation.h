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


#ifndef SEASTATE_COLLOCATION_H_
#define SEASTATE_COLLOCATION_H_

#include<string>
#include<vector>

/*--------------------------------------------------------------------
Online collocation of the model with observations (A 787 1).

seastate_collocation holds the observations of seastate-collocation.csv
and the running statistics; no lexer or MPI dependency (unit test
Regression/unit/seastate_test.cpp). seastate_f samples the model on the
rank that owns the cell of an observation (seastate_f_collocation.cpp).

The file: comma separated, one header line with the column names, then
one observation per line ('#' starts a comment):
  time    YYYYMMDD.HHMMSS or YYYY-MM-DD HH:MM:SS[.sss] (also with T),
          the file time of the forcing (A 780); or
  t       the model time [s] instead of time
  x, y    position in the model frame [m] (tools/seastate_collocation.py
          converts longitude and latitude)
  hs      observed significant wave height [m]
  source  optional: the kind (alti, buoy, ...), statistics per source
  id      optional: mission or station
  tp, tm01, tm10
          optional: observed peak period, mean period Tm01 and energy
          period Tm-1,0 [s]; empty or nan where not observed
Other columns are ignored.

Matchups: the model Hs, Tm01, Tp and mean direction, bilinear between
the four cell centres around the point (all four wet, else the matchup
is flagged dry), linear in time between the ends of the two time steps
around the observation (nonstationary, A 700 1); in a stationary run
(A 700 2) all observations take the final field.

Statistics per source and over all (matchups counted from model time
A 788 on, in a stationary run all; wet only): N, mean observation, bias, RMSE, scatter index,
correlation.

The regimes (Phase 25), Hs per source and over all:
  - by the observed Hs: 0-1, 1-2, 2-4, 4-6, 6-8, >= 8 m;
  - by the model depth at the point: 0-10, 10-20, 20-50, 50-100,
    100-200, >= 200 m;
  - by the wave age cp/U10 of the model (cp at the model Tp and depth,
    U10 of the model wind): < 0.8 young wind sea, 0.8-1.2 developed,
    1.2-2 old, >= 2 (or U10 < 1 m/s) swell;
the periods where observed (tp, tm01, tm10 of the file): the same
statistics of the model Tp, Tm01, Tm-1,0; the quantiles of the observed
and the model Hs (quantile-quantile, each sorted on its own) at 50, 90,
95, 99 and 99.9 % and the largest, per source and over all.
--------------------------------------------------------------------*/

struct seastate_obs
{
    double t;           // model time [s]
    double x, y;        // model frame [m]
    double hs;          // observed Hs [m]
    double tp, tm01, tm10;  // observed periods [s], nan if not observed
    int src, id;        // index into seastate_collocation::sources, ids
};

// a matchup counted in the statistics: observation k, the model at the point
struct seastate_matchup
{
    int src = 0;
    double hs_obs = 0.0, tp_obs = 0.0, tm01_obs = 0.0, tm10_obs = 0.0;      // the observation (periods nan if none)
    double hs = 0.0, tm01 = 0.0, tp = 0.0, tm10 = 0.0, depth = 0.0, u10 = 0.0;  // the model
};

struct seastate_colloc_stats
{
    long n = 0;
    double so = 0.0, sm = 0.0, so2 = 0.0, sm2 = 0.0, som = 0.0;

    void add(double o, double m) {++n; so += o; sm += m; so2 += o*o; sm2 += m*m; som += o*m;}
    void merge(const seastate_colloc_stats &b) {n += b.n; so += b.so; sm += b.sm; so2 += b.so2; sm2 += b.sm2; som += b.som;}
    double mean_obs() const {return n>0 ? so/n : 0.0;}
    double bias() const {return n>0 ? (sm-so)/n : 0.0;}
    double rmse() const;
    double si() const;          // standard deviation of the error over the mean observation
    double r() const;
};

class seastate_collocation
{
public:
    // read the file; times of the time column minus tref (file time of model time 0); observations sorted by time
    bool read(const std::string &file, double tref, std::string &error);

    // YYYYMMDD.HHMMSS or YYYY-MM-DD[ T]HH:MM:SS[.sss] -> seconds since 1970 (UTC); false if neither
    static bool parse_time(const std::string &s, double &t);
    // seconds since 1970 (UTC) -> YYYYMMDD.HHMMSS (whole seconds)
    static std::string stamp(double t);

    std::vector<seastate_obs> obs;
    std::vector<std::string> sources, ids;
    bool model_time = false;    // the file gives t (model seconds)

    bool has_tp = false, has_tm01 = false, has_tm10 = false;   // the file has the column

    // statistics per source (index) and over all
    std::vector<seastate_colloc_stats> per_source;
    seastate_colloc_stats total;
    void count(int src, double o, double m);
    // a matchup: Hs as count(), the regimes, the periods, the quantiles
    void count(const seastate_matchup &m);

    // the summary table (one line per source and the total)
    std::string summary() const;

    // the regimes: the classes of each kind, the statistics [kind][class][source], the last source index all
    static const int nkind = 3;
    static const char *kind_name(int k);
    static int nclass(int k);
    static std::string class_name(int k, int c);
    static int classify(int k, const seastate_matchup &m);
    const seastate_colloc_stats &regime(int k, int c, int src) const;   // src -1: all
    // the periods: 0 Tp, 1 Tm01, 2 Tm-1,0 [period][source], the last all
    const seastate_colloc_stats &period(int q, int src) const;
    // the quantile f (0..1, 1: the largest) of the observed (model false) or model Hs of the counted matchups
    double quantile(int src, double f, bool model) const;

    // the tables of REEF3D_SEASTATE_Collocation (header lines '#', numeric rows)
    std::string regimes_table() const;
    std::string periods_table() const;
    std::string quantiles_table() const;

    // the phase velocity [m/s] at period T [s] and depth d [m] (linear dispersion, deep water for d <= 0)
    static double phase_speed(double T, double d);

private:
    std::vector<seastate_colloc_stats> reg;      // [kind offset + class][source+1]
    std::vector<seastate_colloc_stats> per;      // [period][source+1]
    std::vector<std::vector<float>> qo, qm;      // the counted Hs per source+1: observed, model
    int reg_off(int k) const;
};

#endif
