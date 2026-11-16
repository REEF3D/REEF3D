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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ks_map.h"
#include"lexer.h"
#include<fstream>
#include<sstream>
#include<iostream>
#include<cmath>

using namespace std;

const ks_map& ks_map::get(lexer *p)
{
    static ks_map m;
    static bool read = false;
    if(read)
    return m;
    read = true;

    ifstream in("ks-map.dat");
    if(!in.is_open())
    {
        if(p->mpirank==0)
        cout<<"!!! B 58 1: ./ks-map.dat not found, ks = B 50 everywhere"<<endl;
        return m;
    }
    // the header and the values, comment lines (#) skipped
    stringstream all;
    string line;
    while(getline(in,line))
    {
        const size_t c = line.find('#');
        if(c!=string::npos)
        line.erase(c);
        all<<line<<' ';
    }
    all>>m.x0>>m.y0>>m.dx>>m.dy>>m.nx>>m.ny;
    if(!all || m.nx<1 || m.ny<1 || m.dx<=0.0 || m.dy<=0.0)
    {
        if(p->mpirank==0)
        cout<<"!!! B 58 1: bad header in ./ks-map.dat (x0 y0 dx dy nx ny), ks = B 50 everywhere"<<endl;
        return m;
    }
    m.v.resize((size_t)m.nx*m.ny);
    for(size_t n=0; n<m.v.size(); ++n)
    if(!(all>>m.v[n]))
    {
        if(p->mpirank==0)
        cout<<"!!! B 58 1: ./ks-map.dat has "<<n<<" of "<<m.v.size()<<" values, ks = B 50 everywhere"<<endl;
        m.v.clear();
        return m;
    }
    m.ok = true;
    if(p->mpirank==0)
    {
        double vmin=1e30, vmax=-1e30;
        for(double a : m.v)
        if(a>0.0)
        {
            vmin = min(vmin,a);
            vmax = max(vmax,a);
        }
        cout<<"B 58 1: bed roughness map ./ks-map.dat, "<<m.nx<<" x "<<m.ny<<" points, ks "<<vmin<<" .. "<<vmax<<" m"<<endl;
    }
    return m;
}

double ks_map::at(double x, double y, double def) const
{
    if(!ok)
    return def;
    const double fx = (x-x0)/dx, fy = (y-y0)/dy;
    if(fx<0.0 || fy<0.0 || fx>nx-1 || fy>ny-1)
    return def;
    int i = (int)floor(fx), j = (int)floor(fy);
    i = min(i,max(nx-2,0));
    j = min(j,max(ny-2,0));
    const double wx = (nx>1) ? fx-i : 0.0, wy = (ny>1) ? fy-j : 0.0;
    const int i1 = min(i+1,nx-1), j1 = min(j+1,ny-1);
    const double a = v[(size_t)j*nx+i], b = v[(size_t)j*nx+i1], c = v[(size_t)j1*nx+i], d = v[(size_t)j1*nx+i1];
    if(a<=0.0 || b<=0.0 || c<=0.0 || d<=0.0)
    {
        // a hole in the map: the nearest point, B 50 if that one is a hole as well
        const double q = v[(size_t)(wy<0.5 ? j : j1)*nx+(wx<0.5 ? i : i1)];
        return q>0.0 ? q : def;
    }
    return (1.0-wx)*(1.0-wy)*a + wx*(1.0-wy)*b + (1.0-wx)*wy*c + wx*wy*d;
}
