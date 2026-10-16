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
// REEF3D::CORAL: hydrodynamic database of linear potential flow for N bodies (6N modes):
// hydrostatic restoring C, added mass A(w), radiation damping B(w), infinite-frequency added mass,
// wave excitation X(w, beta) per unit wave amplitude, about one reference point per body, in the
// global frame at equilibrium, time convention exp(i w t) (F = Re{a X exp(i (w t + phase))}).
//
// Sources: CORAL's own text format (read/write, below), WAMIT output (.1, .3, .hst), and later the
// CORAL panel solver and REEF3D forced-oscillation runs.
//
// CORAL format (tokens separated by white space, # starts a comment):
//   coral_hydro_db 1
//   rho <rho>  g <g>  depth <d>                       (d <= 0: infinite depth)
//   bodies <N>
//   body <i> name <name> ref <x y z> [volume <V>] [cob <x y z>]      (i = 1..N)
//   hydrostatics <6N x 6N>
//   frequencies <n> <w_1 .. w_n>                       (rad/s, ascending)
//   headings <m> <beta_1 .. beta_m>                    (degrees)
//   [added_mass_inf <6N x 6N>]
//   frequency <k>  added_mass <6N x 6N>  damping <6N x 6N>
//                  excitation <j> <6N pairs re im>  (one per heading)       (k = 1..n)
// Interpolation in w is linear; outside the frequency range A is kept at the end values and B and X
// are zero.

#ifndef CORAL_HYDRO_DB_H_
#define CORAL_HYDRO_DB_H_

#include"coral/core/coral_types.h"
#include<string>
#include<vector>

namespace coral
{
    using VecXc = Eigen::VectorXcd;

    class hydro_db
    {
    public:
        struct body_info
        {
            std::string name;
            Vec3 ref = Vec3::Zero();     // reference point of the coefficients (global, equilibrium)
            double volume = -1.0;        // displaced volume; < 0: unknown (buoyancy = weight)
            Vec3 cob = Vec3::Zero();     // centre of buoyancy (global, equilibrium)
        };

        int nbody() const {return (int)body.size();}
        int ndof() const {return 6*nbody();}
        int nfreq() const {return (int)w.size();}
        int nheading() const {return (int)beta.size();}

        // set up empty arrays for N bodies, n frequencies, m headings
        void resize(int N, int n, int m);

        MatX added_mass(double omega) const;
        MatX damping(double omega) const;
        VecXc excitation(double omega, double heading) const;   // heading [rad]

        void read(const std::string &file);
        void write(const std::string &file) const;

        // WAMIT output root.1, root.3 (optional), root.hst (optional), nondimensional with the length
        // ulen; periods in the first column (PER = -1: infinite frequency, PER = 0: zero frequency)
        void read_wamit(const std::string &root, int N, double ulen, double rho, double g, double depth=-1.0);

        // consistency of the sizes (throws)
        void check() const;

        double rho = 1025.0;
        double g = 9.81;
        double depth = -1.0;
        std::vector<body_info> body;
        MatX C;                          // hydrostatic restoring
        std::vector<double> w;           // rad/s, ascending
        std::vector<double> beta;        // rad
        std::vector<MatX> A, B;          // per frequency
        MatX Ainf;
        bool has_ainf = false;
        std::vector<std::vector<VecXc>> X;    // [frequency][heading]
    };
}

#endif
