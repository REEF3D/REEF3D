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

#include"sediment_inflow.h"
#include"lexer.h"
#include<fstream>
#include<iostream>
#include<cstdlib>
#include<cstdio>
#include<string>

sediment_inflow::sediment_inflow(lexer *p) : mode(p->S129)
{
    if(mode!=2)
    return;
    
    std::ifstream f("sediment-inflow.dat", std::ios_base::in);
    
    if(!f)
    {
        if(p->mpirank==0)
        std::cout<<"!!! S 129 2: no 'sediment-inflow.dat' file found !!!"<<std::endl;
        
        exit(1);
    }
    
    // t qb C per line; lines starting with # are comments
    std::string line;
    while(std::getline(f,line))
    {
        size_t s = line.find_first_not_of(" \t");
        if(s==std::string::npos || line[s]=='#')
        continue;
        
        double a,b,d;
        if(sscanf(line.c_str(),"%lf %lf %lf",&a,&b,&d)==3)
        {
        t.push_back(a);
        q.push_back(b);
        c.push_back(d);
        }
    }
    
    if(t.empty())
    {
        if(p->mpirank==0)
        std::cout<<"!!! S 129 2: 'sediment-inflow.dat' has no lines 't qb C' !!!"<<std::endl;
        
        exit(1);
    }
    
    if(p->mpirank==0)
    std::cout<<"S 129 2: sediment inflow time series, "<<t.size()<<" lines, t = "<<t.front()<<" - "<<t.back()<<" s"<<std::endl;
}

double sediment_inflow::ipol(lexer *p, int col) const
{
    const std::vector<double> &v = (col==0) ? q : c;
    const double tt = p->simtime;
    
    if(tt<=t.front())
    return v.front();
    
    if(tt>=t.back())
    return v.back();
    
    size_t n=0;
    while(n+1<t.size() && !(tt<t[n+1]))
    ++n;
    
    const double dt = t[n+1]-t[n];
    
    return (dt>0.0) ? v[n] + (v[n+1]-v[n])*(tt-t[n])/dt : v[n+1];
}

double sediment_inflow::qb(lexer *p) const
{
    if(mode==1)
    return p->S130;
    
    if(mode==2)
    return ipol(p,0);
    
    return -1.0;
}

double sediment_inflow::conc(lexer *p) const
{
    if(mode==2)
    return ipol(p,1);
    
    return p->S124;
}
