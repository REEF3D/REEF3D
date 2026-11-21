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

#include"coral/pto/coral_generator.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>

namespace coral
{

generator_power generator_model::convert(double F, double qd) const
{
    const generator_settings &s = settings;
    generator_power p;
    p.P_mech = -F*qd;
    const double g = (p.P_mech >= 0.0) ? s.eta_gear : 1.0/s.eta_gear;
    p.w = s.ratio*qd;
    p.T = -F/s.ratio*g;
    p.P_shaft = p.T*p.w;
    switch(s.kind)
    {
    case generator_kind::ideal:
        p.P_el = (p.P_shaft >= 0.0) ? s.eta*p.P_shaft : p.P_shaft/s.eta;
        break;
    case generator_kind::map:
    {
        const double e = map_eta(p.w, p.T);
        p.P_el = (p.P_shaft >= 0.0) ? e*p.P_shaft : p.P_shaft/e;
        break;
    }
    case generator_kind::pmsg:
        p.I = std::fabs(p.T)/s.kt;
        p.P_el = p.P_shaft - 1.5*s.R*p.I*p.I - s.k_fe*p.w*p.w;
        break;
    }
    p.P_grid = (p.P_el >= 0.0) ? s.eta_conv*p.P_el : p.P_el/s.eta_conv;
    return p;
}

namespace
{
    // index n and weight r with x between v[n] and v[n+1], the end values outside
    void bracket(const std::vector<double> &v, double x, size_t &n, double &r)
    {
        if(v.size() < 2 || x <= v.front()) {n = 0; r = 0.0; return;}
        if(x >= v.back()) {n = v.size() - 2; r = 1.0; return;}
        n = size_t(std::upper_bound(v.begin(), v.end(), x) - v.begin()) - 1;
        r = (x - v[n])/(v[n+1] - v[n]);
    }
}

double generator_model::map_eta(double w, double T) const
{
    const generator_settings &s = settings;
    const size_t nT = s.map_T.size();
    if(s.map_w.empty() || nT == 0)
        return 1.0;
    size_t i = 0, j = 0;
    double ri = 0.0, rj = 0.0;
    bracket(s.map_w, std::fabs(w), i, ri);
    bracket(s.map_T, std::fabs(T), j, rj);
    auto e = [&](size_t a, size_t b) {return s.map_eta[std::min(a, s.map_w.size()-1)*nT + std::min(b, nT-1)];};
    const double e0 = (1.0 - rj)*e(i, j) + rj*e(i, j+1);
    const double e1 = (1.0 - rj)*e(i+1, j) + rj*e(i+1, j+1);
    return std::max(1.0e-3, (1.0 - ri)*e0 + ri*e1);
}

bool generator_model::read_map(const std::string &file, generator_settings &s, std::string &msg)
{
    std::ifstream in(file);
    if(!in)
    {
        msg = "cannot open '" + file + "'";
        return false;
    }
    s.map_w.clear(); s.map_T.clear(); s.map_eta.clear();
    std::string line;
    int n = 0;
    bool head = true;
    while(std::getline(in, line))
    {
        ++n;
        const size_t h = line.find('#');
        if(h != std::string::npos) line.erase(h);
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream ls(line);
        std::vector<double> v;
        double x;
        while(ls >> x) v.push_back(x);
        if(v.empty())
            continue;
        if(head)
        {
            if(v.size() < 2)
            {
                msg = file + " line " + std::to_string(n) + ": the first line is 0 T_1 .. T_m";
                return false;
            }
            s.map_T.assign(v.begin()+1, v.end());
            head = false;
            continue;
        }
        if(v.size() != s.map_T.size() + 1)
        {
            msg = file + " line " + std::to_string(n) + ": w and " + std::to_string(s.map_T.size()) + " efficiencies expected";
            return false;
        }
        s.map_w.push_back(v[0]);
        for(size_t k=1; k<v.size(); ++k)
        {
            if(v[k] <= 0.0 || v[k] > 1.0)
            {
                msg = file + " line " + std::to_string(n) + ": efficiencies must be in (0, 1]";
                return false;
            }
            s.map_eta.push_back(v[k]);
        }
    }
    auto increasing = [](const std::vector<double> &v)
    {
        for(size_t k=1; k<v.size(); ++k) if(v[k] <= v[k-1]) return false;
        return true;
    };
    if(s.map_w.empty() || !increasing(s.map_w) || !increasing(s.map_T))
    {
        msg = file + ": speeds and torques must increase, at least one speed line";
        return false;
    }
    s.kind = generator_kind::map;
    return true;
}

std::string generator_model::check() const
{
    const generator_settings &s = settings;
    if(s.ratio == 0.0) return "the gear ratio must not be 0";
    if(s.eta_gear <= 0.0 || s.eta_gear > 1.0) return "eta_gear must be in (0, 1]";
    if(s.eta_conv <= 0.0 || s.eta_conv > 1.0) return "eta_conv must be in (0, 1]";
    if(s.kind == generator_kind::ideal && (s.eta <= 0.0 || s.eta > 1.0)) return "eta must be in (0, 1]";
    if(s.kind == generator_kind::pmsg && (s.kt <= 0.0 || s.R < 0.0 || s.k_fe < 0.0)) return "pmsg needs kt > 0, R >= 0, k_fe >= 0";
    if(s.kind == generator_kind::map && s.map_w.empty()) return "the efficiency map is empty";
    return "";
}

}
