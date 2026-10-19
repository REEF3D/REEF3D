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

#include"coral/io/coral_output.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include"coral/core/coral_version.h"
#include<cmath>
#include<cstdio>
#include<stdexcept>

namespace coral
{

void time_series::open(const std::string &file, const model &m, bool append)
{
    // a restart continues an existing file (no new header)
    if(append)
    {
        std::ifstream test(file);
        append = test.good();
    }
    if(append)
        out_.open(file, std::ios::app);
    else
        out_.open(file);
    if(!out_)
        throw std::runtime_error("CORAL: cannot write " + file);
    cols_.clear();
    cols_.push_back("t");
    for(const auto &b : m.bodies)
        for(const char *c : {"x", "y", "z", "roll", "pitch", "yaw", "vx", "vy", "vz", "wx", "wy", "wz"})
            cols_.push_back(b.name + "_" + c);
    for(const auto &l : m.load_list)
        l->channel_names(cols_);
    if(append)
        return;
    out_<<"# "<<banner()<<"\n#";
    for(const auto &c : cols_)
        out_<<" "<<c;
    out_<<"\n";
}

void time_series::write(const model &m, double t)
{
    std::vector<double> v;
    v.push_back(t);
    for(const auto &b : m.bodies)
    {
        const Vec3 e = euler_from_matrix(b.R)*180.0/M_PI;
        const Vec3 w = b.omega_global();
        for(int i=0; i<3; ++i) v.push_back(b.x(i));
        for(int i=0; i<3; ++i) v.push_back(e(i));
        for(int i=0; i<3; ++i) v.push_back(b.v(i));
        for(int i=0; i<3; ++i) v.push_back(w(i));
    }
    for(const auto &l : m.load_list)
        l->channel_values(m, t, v);
    char buf[32];
    for(size_t i=0; i<v.size(); ++i)
    {
        std::snprintf(buf, sizeof(buf), i ? " %.10e" : "%.10e", v[i]);
        out_<<buf;
    }
    // flushed per line: the file is complete when a coupled run ends without closing it
    out_<<"\n"<<std::flush;
}

}
