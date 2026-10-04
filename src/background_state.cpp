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

#include"background_state.h"
#include"lexer.h"
#include<cmath>
#include<cstdlib>
#include<iostream>
#include<string>

void background_state::read(lexer *p)
{
    auto fail = [&](const std::string &msg)
    {
        if(p->mpirank==0)
        std::cout<<std::endl<<"!!! background: "<<msg<<" !!!"<<std::endl<<std::endl;
        std::exit(1);
    };
    
    g = fabs(p->W22)>0.0 ? fabs(p->W22) : 9.81;
    
    for(int n=0; n<p->B510; ++n)
    {
        if(p->B510_id[n]<1)
        fail("B 510: background ids start at 1");
        if(index(p->B510_id[n])>=0)
        fail("background "+std::to_string(p->B510_id[n])+" defined twice in B 510");
        if(p->B510_mode[n]!=1 && p->B510_mode[n]!=3)
        fail("B 510 mode is 1 (harmonic) or 3 (constant); time series (2) are not available yet");
        
        item b;
        b.id = p->B510_id[n];
        b.mode = p->B510_mode[n];
        b.dir = p->B510_dir[n]*(3.14159265358979323846/180.0);
        b.tramp = p->B510_tramp[n];
        bg.push_back(b);
    }
    
    for(int n=0; n<p->B511; ++n)
    {
        int b = index(p->B511_id[n]);
        if(b<0)
        fail("B 511 refers to background "+std::to_string(p->B511_id[n])+", which has no B 510");
        if(bg[b].mode!=1)
        fail("B 511: background "+std::to_string(p->B511_id[n])+" is not harmonic (B 510 mode 1)");
        if(p->B511_T[n]<=0.0)
        fail("B 511: the period must be positive");
        
        bg[b].c.push_back({p->B511_a[n], p->B511_T[n], p->B511_phase[n]*(3.14159265358979323846/180.0)});
    }
    
    for(int n=0; n<p->B514; ++n)
    {
        int b = index(p->B514_id[n]);
        if(b<0)
        fail("B 514 refers to background "+std::to_string(p->B514_id[n])+", which has no B 510");
        
        bg[b].eta0 = p->B514_eta0[n];
        bg[b].U = p->B514_U[n];
        bg[b].V = p->B514_V[n];
    }
    
    if(p->mpirank==0)
    for(const item &b : bg)
    {
        std::cout<<"background "<<b.id<<": "<<(b.mode==1 ? "harmonic, " : "constant, ")<<b.c.size()<<" constituents, eta0 "<<b.eta0
                 <<" U "<<b.U<<" V "<<b.V<<", t_ramp "<<b.tramp<<std::endl;
    }
}

int background_state::index(int id) const
{
    for(int b=0; b<int(bg.size()); ++b)
    if(bg[b].id==id)
    return b;
    
    return -1;
}

void background_state::update(lexer *p, double t)
{
    const double pi = 3.14159265358979323846;
    
    for(item &b : bg)
    {
        b.r = 1.0;
        if(b.tramp>0.0 && t<b.tramp)
        b.r = 0.5*(1.0 - cos(pi*fmax(t,0.0)/b.tramp));
        
        double s = 0.0;
        for(const constituent &c : b.c)
        s += c.a*cos(2.0*pi*t/c.T - c.phase);
        
        b.eta_h = b.r*s;
        b.eta = b.r*b.eta0 + b.eta_h;
    }
}

void background_state::vel(int k, double h, double &u, double &v) const
{
    const item &b = bg[k];
    const double c = h>1.0e-6 ? sqrt(g/h) : 0.0;
    
    u = b.r*b.U + b.eta_h*c*cos(b.dir);
    v = b.r*b.V + b.eta_h*c*sin(b.dir);
}
