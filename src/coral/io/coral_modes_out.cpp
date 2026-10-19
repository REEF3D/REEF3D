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
#include"coral/io/coral_modes_out.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_version.h"
#include"coral/linear/coral_modes.h"
#include<cmath>
#include<cstdio>
#include<fstream>

namespace coral
{
    std::array<double,6> energy_split(const model &m, const mode_result &r, int q)
    {
        std::array<double,6> e = {0, 0, 0, 0, 0, 0};
        const VecX &s = r.shape.col(q);
        for(int b=0; b<m.nbody(); ++b)
        {
            const rigid_body &bd = m.bodies[b];
            const Vec3 v = s.segment<3>(6*b);
            const Vec3 wg = bd.R*s.segment<3>(6*b+3);
            const Vec3 Iw = bd.R*bd.inertia*bd.R.transpose()*wg;
            for(int k=0; k<3; ++k)
            {
                e[k] += bd.mass*v(k)*v(k);
                e[3+k] += wg(k)*Iw(k);
            }
        }
        double sum = 0.0;
        for(double x : e) sum += std::fabs(x);
        if(sum > 0.0)
            for(double &x : e) x = std::fabs(x)/sum;
        return e;
    }

    namespace
    {
        int largest(const std::array<double,6> &e)
        {
            int k = 0;
            for(int i=1; i<6; ++i) if(e[i] > e[k]) k = i;
            return k;
        }
    }

    std::string dominant(const std::array<double,6> &e)
    {
        static const char *nm[6] = {"x", "y", "z", "rx", "ry", "rz"};
        return nm[largest(e)];
    }

    void write_modes(const std::string &dir, const model &m, const mode_result &r)
    {
        std::ofstream mf(dir + "modes.dat");
        mf<<"# "<<banner()<<"\n# kinetic energy split of each mode over the global directions\n"
          <<"# mode frequency[Hz] period[s] damping_ratio Ex Ey Ez Erx Ery Erz dominant\n";
        char line[256];
        for(int q=0; q<r.frequency.size(); ++q)
        {
            const double f = r.frequency(q);
            const std::array<double,6> e = energy_split(m, r, q);
            std::snprintf(line, sizeof(line), "%d %.9g %.9g %.6g %.4f %.4f %.4f %.4f %.4f %.4f %s\n", q+1, f, f > 0.0 ? 1.0/f : 0.0,
                          r.damping(q), e[0], e[1], e[2], e[3], e[4], e[5], dominant(e).c_str());
            mf<<line;
        }
        std::ofstream sf(dir + "mode_shapes.dat");
        sf<<"# "<<banner()<<"\n# per mode: body x y z ux uy uz rx ry rz (translation and rotation global, largest entry 1)\n";
        for(int q=0; q<r.frequency.size(); ++q)
        {
            std::snprintf(line, sizeof(line), "# mode %d %.9g Hz\n", q+1, r.frequency(q));
            sf<<line;
            for(int b=0; b<m.nbody(); ++b)
            {
                const rigid_body &bd = m.bodies[b];
                const Vec3 u = r.shape.col(q).segment<3>(6*b);
                const Vec3 w = bd.R*r.shape.col(q).segment<3>(6*b+3);
                std::snprintf(line, sizeof(line), "%s %.6f %.6f %.6f % .6e % .6e % .6e % .6e % .6e % .6e\n",
                              bd.name.empty() ? std::to_string(b).c_str() : bd.name.c_str(), bd.x(0), bd.x(1), bd.x(2),
                              u(0), u(1), u(2), w(0), w(1), w(2));
                sf<<line;
            }
            sf<<"\n\n";
        }
    }

    void print_modes(const model &m, const mode_result &r, bool damped)
    {
        std::printf("linear analysis (%s):\n", damped ? "damped" : "undamped");
        for(int q=0; q<r.frequency.size(); ++q)
        {
            const double f = r.frequency(q);
            const std::array<double,6> e = energy_split(m, r, q);
            std::printf("  mode %3d  %12.6f Hz  %12.6f s  zeta %.4f  %-2s (%3.0f%% of the kinetic energy)\n", q+1, f, f > 0.0 ? 1.0/f : 0.0,
                        r.damping(q), dominant(e).c_str(), 100.0*e[largest(e)]);
        }
    }
}
