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


#include"seastate_bathy.h"
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>

namespace
{
    // next non-empty line without the comment
    bool next_line(std::ifstream &in, std::string &s)
    {
        while(std::getline(in,s))
        {
        const size_t c = s.find('$');
        if(c!=std::string::npos)
        s.erase(c);

        if(s.find_first_not_of(" \t\r")!=std::string::npos)
        return true;
        }

        return false;
    }
}

bool seastate_bathy::read(const std::string &file, std::string &err)
{
    std::ifstream in(file.c_str());

    if(!in)
    {
    err = "cannot open " + file;
    return false;
    }

    std::string s;

    // title (may be empty or a comment)
    if(!std::getline(in,s))
    {
    err = file + ": empty file";
    return false;
    }

    if(!next_line(in,s))
    {
    err = file + ": missing line 'nx ny'";
    return false;
    }

        {
        std::istringstream is(s);
        if(!(is>>nx>>ny) || nx<2 || ny<2)
        {
        err = file + ": 'nx ny' must be at least 2 2";
        return false;
        }
        }

    if(!next_line(in,s))
    {
    err = file + ": missing line 'x0 y0 dx dy'";
    return false;
    }

        {
        std::istringstream is(s);
        if(!(is>>x0>>y0>>dx>>dy) || !(dx>0.0) || !(dy>0.0))
        {
        err = file + ": 'x0 y0 dx dy' with positive spacing expected";
        return false;
        }
        }

    z.assign(size_t(nx)*ny,0.0);
    size_t n = 0;

    while(n<z.size() && next_line(in,s))
    {
    std::istringstream is(s);
    double v;

        while(n<z.size() && is>>v)
        z[n++] = v;
    }

    if(n<z.size())
    {
    std::ostringstream o;
    o<<file<<": "<<n<<" of "<<z.size()<<" bed levels read";
    err = o.str();
    return false;
    }

    return true;
}

double seastate_bathy::at(double x, double y) const
{
    const double fx = std::min(std::max((x-x0)/dx,0.0),double(nx-1));
    const double fy = std::min(std::max((y-y0)/dy,0.0),double(ny-1));

    const int i = std::min(int(fx),nx-2);
    const int j = std::min(int(fy),ny-2);
    const double a = fx-i, b = fy-j;

    return (1.0-a)*(1.0-b)*node(i,j) + a*(1.0-b)*node(i+1,j) + (1.0-a)*b*node(i,j+1) + a*b*node(i+1,j+1);
}

double seastate_bathy::cell(double xa, double xb, double ya, double yb, double zdry) const
{
    // nodes with xa <= x < xb, ya <= y < yb
    const int i0 = std::max(int(std::ceil((xa-x0)/dx - 1.0e-9)),0);
    const int i1 = std::min(int(std::ceil((xb-x0)/dx - 1.0e-9))-1,nx-1);
    const int j0 = std::max(int(std::ceil((ya-y0)/dy - 1.0e-9)),0);
    const int j1 = std::min(int(std::ceil((yb-y0)/dy - 1.0e-9))-1,ny-1);

    if(i1<i0 || j1<j0)
    return at(0.5*(xa+xb),0.5*(ya+yb));

    double sw=0.0, sd=0.0;
    int nw=0, nd=0;

    for(int j=j0; j<=j1; ++j)
    for(int i=i0; i<=i1; ++i)
    {
        const double z = node(i,j);
        if(z>zdry)
        {
        sd += z;
        ++nd;
        }
        else
        {
        sw += z;
        ++nw;
        }
    }

    // majority of the nodes (ties: wet); without zdry all nodes are wet
    if(nd>nw)
    return sd/double(nd);

    return sw/double(nw);
}
