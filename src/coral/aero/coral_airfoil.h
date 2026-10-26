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
// REEF3D::CORAL: airfoil polars: lift, drag and moment coefficients against the angle of attack, one or more tables
// per airfoil (Reynolds numbers), linear in the angle of attack (periodic when a table spans -180 .. 180 deg, held
// at its ends otherwise) and between the tables in log(Re) (held outside them).
// Readers: AeroDyn's airfoil file (AirfoilInfo v1, as in OpenFAST: NumTabs, per table Re, the unsteady
// aerodynamics parameters (kept by name for the dynamic stall model) and NumAlf rows alpha cl cd [cm] in deg),
// or a plain table of rows alpha cl cd [cm] in deg ('#' and '!' start comments).

#ifndef CORAL_AIRFOIL_H_
#define CORAL_AIRFOIL_H_

#include<map>
#include<string>
#include<vector>

namespace coral
{
    struct polar_table
    {
        double re = 0.0;                         // Reynolds number [-] (AeroDyn files: in millions, converted)
        std::vector<double> alpha, cl, cd, cm;   // alpha in rad, increasing
        std::map<std::string, double> params;    // the unsteady aerodynamics parameters of the file (alpha0, Cn1, ...)
    };

    struct airfoil_coef {double cl = 0.0, cd = 0.0, cm = 0.0;};

    class airfoil
    {
    public:
        std::string name;
        std::vector<polar_table> tables;         // increasing Re

        // coefficients at alpha [rad] and the Reynolds number (ignored with one table)
        airfoil_coef coefficients(double alpha, double re=0.0) const;
        airfoil_coef table_coefficients(int i, double alpha) const;

        // AeroDyn airfoil file or plain table; false with a message
        static bool read(const std::string &file, airfoil&, std::string &msg);
        static bool read_string(const std::string &text, airfoil&, std::string &msg);
    };

    // a in (-pi, pi]
    double wrap_angle(double a);
}

#endif
