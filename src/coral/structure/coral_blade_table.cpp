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

#include"coral/structure/coral_blade_table.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>

namespace coral
{

double blade_structure::mass_at(double f) const
{
    if(frac.empty()) return 0.0;
    if(f <= frac.front()) return mass.front();
    if(f >= frac.back()) return mass.back();
    const size_t k = std::upper_bound(frac.begin(), frac.end(), f) - frac.begin() - 1;
    const double s = (f - frac[k])/(frac[k+1] - frac[k]);
    return mass[k] + s*(mass[k+1] - mass[k]);
}

double blade_structure::total_mass(double L) const
{
    double m = 0.0;
    for(size_t k=0; k+1<frac.size(); ++k)
        m += 0.5*(mass[k] + mass[k+1])*(frac[k+1] - frac[k])*L;
    return m;
}

bool blade_structure::read(const std::string &file, blade_structure &b, std::string &msg)
{
    std::ifstream in(file);
    if(!in)
    {
        msg = "cannot open '" + file + "'";
        return false;
    }
    std::stringstream ss;
    ss << in.rdbuf();
    if(!read_string(ss.str(), b, msg))
    {
        msg = file + ": " + msg;
        return false;
    }
    return true;
}

bool blade_structure::read_string(const std::string &text, blade_structure &b, std::string &msg)
{
    b = blade_structure();
    std::istringstream in(text);
    std::string line;
    int n = -1;
    double adj_m = 1.0, adj_f = 1.0, adj_e = 1.0;
    bool elastodyn = text.find("NBlInpSt") != std::string::npos;
    auto num = [](const std::string &w, double &x) {try {size_t k; x = std::stod(w, &k); return k==w.size();} catch(...) {return false;}};
    while(std::getline(in, line))
    {
        const size_t h = line.find('#');
        if(h != std::string::npos) line.erase(h);
        std::istringstream ls(line);
        std::vector<std::string> w;
        std::string x;
        while(ls >> x) w.push_back(x);
        if(w.empty()) continue;
        double v;
        if(elastodyn && w.size() >= 2 && num(w[0], v))
        {
            if(w[1]=="NBlInpSt") {n = (int)v; continue;}
            if(w[1]=="AdjBlMs") {adj_m = v; continue;}
            if(w[1]=="AdjFlSt") {adj_f = v; continue;}
            if(w[1]=="AdjEdSt") {adj_e = v; continue;}
        }
        if(w.size() < 6) continue;
        double r[6];
        bool ok = true;
        for(int q=0; q<6 && ok; ++q) ok = num(w[q], r[q]);
        if(!ok) continue;
        if(elastodyn && (n < 0 || (int)b.frac.size()==n)) continue;   // the mode shape tables after the stations
        if(!b.frac.empty() && r[0] <= b.frac.back())
        {
            msg = "the blade fractions must increase";
            return false;
        }
        b.frac.push_back(r[0]); b.pitch_axis.push_back(r[1]); b.twist.push_back(r[2]*M_PI/180.0);
        b.mass.push_back(r[3]); b.EI_flap.push_back(r[4]); b.EI_edge.push_back(r[5]);
    }
    if(elastodyn && (int)b.frac.size() != n)
    {
        msg = std::to_string(b.frac.size()) + " of " + std::to_string(n) + " blade stations found";
        return false;
    }
    if(b.frac.size() < 2)
    {
        msg = "at least two blade stations needed";
        return false;
    }
    for(size_t k=0; k<b.frac.size(); ++k) {b.mass[k] *= adj_m; b.EI_flap[k] *= adj_f; b.EI_edge[k] *= adj_e;}
    return true;
}

}
