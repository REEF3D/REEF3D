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

#include"coral/aero/coral_airfoil.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>

namespace coral
{

double wrap_angle(double a)
{
    a = std::fmod(a + M_PI, 2.0*M_PI);
    if(a <= 0.0) a += 2.0*M_PI;
    return a - M_PI;
}

airfoil_coef airfoil::table_coefficients(int i, double alpha) const
{
    const polar_table &p = tables[i];
    airfoil_coef c;
    const int n = (int)p.alpha.size();
    if(n==0)
        return c;
    if(n==1)
        return airfoil_coef{p.cl[0], p.cd[0], p.cm.empty() ? 0.0 : p.cm[0]};
    const bool full = p.alpha.front() <= -M_PI + 1e-6 && p.alpha.back() >= M_PI - 1e-6;
    double a = full ? wrap_angle(alpha) : alpha;
    a = std::clamp(a, p.alpha.front(), p.alpha.back());
    int k = (int)(std::upper_bound(p.alpha.begin(), p.alpha.end(), a) - p.alpha.begin()) - 1;
    k = std::clamp(k, 0, n - 2);
    const double s = (a - p.alpha[k])/(p.alpha[k+1] - p.alpha[k]);
    c.cl = p.cl[k] + s*(p.cl[k+1] - p.cl[k]);
    c.cd = p.cd[k] + s*(p.cd[k+1] - p.cd[k]);
    if(!p.cm.empty())
        c.cm = p.cm[k] + s*(p.cm[k+1] - p.cm[k]);
    return c;
}

airfoil_coef airfoil::coefficients(double alpha, double re) const
{
    if(tables.empty())
        return airfoil_coef{};
    if(tables.size()==1 || re <= tables.front().re)
        return table_coefficients(0, alpha);
    if(re >= tables.back().re)
        return table_coefficients((int)tables.size() - 1, alpha);
    int k = 0;
    while(k + 2 < (int)tables.size() && tables[k+1].re < re) ++k;
    const double s = (std::log(re) - std::log(tables[k].re))/(std::log(tables[k+1].re) - std::log(tables[k].re));
    const airfoil_coef a = table_coefficients(k, alpha), b = table_coefficients(k + 1, alpha);
    return airfoil_coef{a.cl + s*(b.cl - a.cl), a.cd + s*(b.cd - a.cd), a.cm + s*(b.cm - a.cm)};
}

bool airfoil::read(const std::string &file, airfoil &af, std::string &msg)
{
    std::ifstream in(file);
    if(!in)
    {
        msg = "cannot open '" + file + "'";
        return false;
    }
    std::stringstream ss;
    ss << in.rdbuf();
    if(!read_string(ss.str(), af, msg))
    {
        msg = file + ": " + msg;
        return false;
    }
    if(af.name.empty())
    {
        const size_t s = file.find_last_of('/');
        std::string n = s==std::string::npos ? file : file.substr(s + 1);
        const size_t d = n.find_last_of('.');
        af.name = d==std::string::npos ? n : n.substr(0, d);
    }
    return true;
}

namespace
{
    bool number(const std::string &w, double &x)
    {
        try {size_t n; x = std::stod(w, &n); return n==w.size();}
        catch(...) {return false;}
    }
}

bool airfoil::read_string(const std::string &text, airfoil &af, std::string &msg)
{
    af.tables.clear();
    std::istringstream in(text);
    std::vector<std::vector<std::string>> lines;
    std::string raw;
    while(std::getline(in, raw))
    {
        for(const char c : {'!', '#'})
        {
            const size_t h = raw.find(c);
            if(h != std::string::npos) raw.erase(h);
        }
        std::istringstream ls(raw);
        std::vector<std::string> w;
        std::string x;
        while(ls >> x) w.push_back(x);
        if(!w.empty()) lines.push_back(w);
    }
    auto rows = [&](size_t &i, int n, polar_table &p) -> bool
    {
        for(int r=0; r<n; ++r, ++i)
        {
            if(i >= lines.size())
            {
                msg = "the table ends after " + std::to_string(r) + " of " + std::to_string(n) + " rows";
                return false;
            }
            double v[4] = {0, 0, 0, 0};
            const size_t m = std::min<size_t>(lines[i].size(), 4);
            for(size_t q=0; q<m; ++q)
                if(!number(lines[i][q], v[q]))
                {
                    msg = "row " + std::to_string(r + 1) + " of the table is not numeric";
                    return false;
                }
            if(m < 3)
            {
                msg = "row " + std::to_string(r + 1) + ": alpha cl cd expected";
                return false;
            }
            const double a = v[0]*M_PI/180.0;
            if(!p.alpha.empty() && a <= p.alpha.back())
            {
                msg = "the angles of attack must increase (row " + std::to_string(r + 1) + ")";
                return false;
            }
            p.alpha.push_back(a); p.cl.push_back(v[1]); p.cd.push_back(v[2]);
            if(m==4) p.cm.push_back(v[3]);
            else if(!p.cm.empty()) {msg = "cm missing in row " + std::to_string(r + 1); return false;}
        }
        return true;
    };
    bool aerodyn = false;
    for(const auto &l : lines)
        if(l.size() >= 2 && l[1]=="NumAlf") aerodyn = true;
    if(aerodyn)
    {
        polar_table *cur = nullptr;
        for(size_t i=0; i<lines.size(); )
        {
            const auto &l = lines[i];
            if(l.size() >= 2 && l[1]=="Re")
            {
                af.tables.emplace_back();
                cur = &af.tables.back();
                double re;
                if(number(l[0], re)) cur->re = re*1.0e6;
                ++i;
            }
            else if(l.size() >= 2 && l[1]=="NumAlf")
            {
                if(!cur) {af.tables.emplace_back(); cur = &af.tables.back();}
                double n;
                if(!number(l[0], n) || n < 1) {msg = "NumAlf is not a positive number"; return false;}
                ++i;
                if(!rows(i, (int)n, *cur)) return false;
                cur = nullptr;
            }
            else
            {
                double v;
                if(cur && l.size() >= 2 && number(l[0], v))
                    cur->params[l[1]] = v;
                ++i;
            }
        }
    }
    else
    {
        af.tables.emplace_back();
        size_t i = 0;
        int n = 0;
        while(i + n < lines.size()) ++n;
        if(!rows(i, n, af.tables.back())) return false;
    }
    af.tables.erase(std::remove_if(af.tables.begin(), af.tables.end(), [](const polar_table &p) {return p.alpha.empty();}), af.tables.end());
    if(af.tables.empty())
    {
        msg = "no polar table";
        return false;
    }
    std::stable_sort(af.tables.begin(), af.tables.end(), [](const polar_table &a, const polar_table &b) {return a.re < b.re;});
    return true;
}

}
