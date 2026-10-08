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
#include<algorithm>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<iostream>
#include<sstream>
#include<string>

void background_state::read(lexer *p)
{
    const double pi = 3.14159265358979323846;
    
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
        if(p->B510_mode[n]<1 || p->B510_mode[n]>3)
        fail("B 510 mode is 1 (harmonic), 2 (time series) or 3 (constant)");
        
        item b;
        b.id = p->B510_id[n];
        b.mode = p->B510_mode[n];
        b.dir = p->B510_dir[n]*(pi/180.0);
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
        
        constituent c;
        c.a = p->B511_a[n];
        c.T = p->B511_T[n];
        c.phase = p->B511_phase[n]*(pi/180.0);
        c.k = 0.0;
        bg[b].c.push_back(c);
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
    
    for(int n=0; n<p->B513; ++n)
    {
        int b = index(p->B513_id[n]);
        if(b<0)
        fail("B 513 refers to background "+std::to_string(p->B513_id[n])+", which has no B 510");
        if(p->B513_profile[n]<0 || p->B513_profile[n]>2)
        fail("B 513 profile is 0 (uniform), 1 (log-law) or 2 (power law)");
        if(p->B513_profile[n]==1 && p->B513_par[n]<=0.0)
        fail("B 513: the log-law needs a roughness length z0 > 0");
        
        bg[b].prof = p->B513_profile[n];
        bg[b].ppar = p->B513_profile[n]==2 && p->B513_par[n]<=0.0 ? 7.0 : p->B513_par[n];
    }
    
    for(int n=0; n<p->B515; ++n)
    {
        int b = index(p->B515_id[n]);
        if(b<0)
        fail("B 515 refers to background "+std::to_string(p->B515_id[n])+", which has no B 510");
        if(bg[b].mode==3)
        fail("B 515: background "+std::to_string(p->B515_id[n])+" is constant (B 510 mode 3), it does not travel");
        
        item &it = bg[b];
        it.prog = true;
        it.x0 = p->B515_x0[n];
        it.y0 = p->B515_y0[n];
        it.href = p->B515_href[n]>0.0 ? p->B515_href[n] : p->wd;
        if(it.href<=0.0)
        fail("B 515: the reference depth must be positive (h_ref or F 60)");
        it.cel = sqrt(g*it.href);
    }
    
    for(item &b : bg)
    {
        // long-wave wave numbers of the constituents
        for(constituent &c : b.c)
        c.k = b.prog ? 2.0*pi/(c.T*b.cel) : 0.0;
        
        // time series: background-<id>.dat, "t eta" or "t eta U V" per line, '#' comments
        if(b.mode==2)
        {
            const std::string name = "background-"+std::to_string(b.id)+".dat";
            std::ifstream f(name);
            if(!f)
            fail("background "+std::to_string(b.id)+" (B 510 mode 2) needs the file "+name);
            
            std::string line;
            int cols=-1;
            while(std::getline(f,line))
            {
                const size_t h = line.find('#');
                if(h!=std::string::npos)
                line.erase(h);
                
                std::istringstream is(line);
                std::vector<double> v;
                double x;
                while(is>>x)
                v.push_back(x);
                
                if(v.empty())
                continue;
                if(v.size()!=2 && v.size()!=4)
                fail(name+": each line needs 2 (t eta) or 4 (t eta U V) values");
                if(cols>=0 && int(v.size())!=cols)
                fail(name+": all lines need the same number of values");
                if(!b.ft.empty() && v[0]<=b.ft.back())
                fail(name+": the times must increase");
                
                cols = int(v.size());
                b.ft.push_back(v[0]);
                b.feta.push_back(v[1]);
                if(cols==4)
                {
                b.fu.push_back(v[2]);
                b.fv.push_back(v[3]);
                }
            }
            
            if(b.ft.size()<2)
            fail(name+": at least two lines are needed");
            
            b.file_uv = (cols==4);
        }
    }
    
    if(p->mpirank==0)
    for(const item &b : bg)
    {
        std::cout<<"background "<<b.id<<": "<<(b.mode==1 ? "harmonic, " : b.mode==2 ? "time series, " : "constant, ");
        if(b.mode==1)
        std::cout<<b.c.size()<<" constituents, ";
        if(b.mode==2)
        std::cout<<b.ft.size()<<" times "<<b.ft.front()<<" - "<<b.ft.back()<<" s"<<(b.file_uv ? " with U, V, " : ", ");
        std::cout<<"eta0 "<<b.eta0<<" U "<<b.U<<" V "<<b.V<<", t_ramp "<<b.tramp;
        if(b.prog)
        std::cout<<", progressive from ("<<b.x0<<", "<<b.y0<<"), h_ref "<<b.href;
        std::cout<<std::endl;
    }
}

int background_state::index(int id) const
{
    for(int b=0; b<int(bg.size()); ++b)
    if(bg[b].id==id)
    return b;
    
    return -1;
}

double background_state::interp(const std::vector<double> &t, const std::vector<double> &f, double tt)
{
    if(tt<=t.front())
    return f.front();
    
    if(tt>=t.back())
    return f.back();
    
    const size_t n = size_t(std::upper_bound(t.begin(),t.end(),tt) - t.begin());
    const double w = (tt-t[n-1])/(t[n]-t[n-1]);
    
    return (1.0-w)*f[n-1] + w*f[n];
}

void background_state::update(lexer *p, double t)
{
    const double pi = 3.14159265358979323846;
    time = t;
    
    for(item &b : bg)
    {
        b.r = 1.0;
        if(b.tramp>0.0 && t<b.tramp)
        b.r = 0.5*(1.0 - cos(pi*fmax(t,0.0)/b.tramp));
        
        for(constituent &c : b.c)
        {
            c.C = cos(2.0*pi*t/c.T - c.phase);
            c.S = sin(2.0*pi*t/c.T - c.phase);
        }
    }
}

double background_state::tide(const item &b, double x, double y, double &uf, double &vf) const
{
    const double s = b.prog ? (x-b.x0)*cos(b.dir) + (y-b.y0)*sin(b.dir) : 0.0;
    double e = 0.0;
    uf = vf = 0.0;
    
    // cos(w t - phase - k s) = C cos(k s) + S sin(k s)
    if(b.mode==1)
    for(const constituent &c : b.c)
    e += c.a*(c.C*cos(c.k*s) + c.S*sin(c.k*s));
    
    if(b.mode==2)
    {
        const double tt = b.prog ? time - s/b.cel : time;
        e = interp(b.ft,b.feta,tt);
        
        if(b.file_uv)
        {
        uf = b.r*interp(b.ft,b.fu,tt);
        vf = b.r*interp(b.ft,b.fv,tt);
        }
    }
    
    return b.r*e;
}

double background_state::eta(int k, double x, double y) const
{
    const item &b = bg[k];
    double uf, vf;
    
    return b.r*b.eta0 + tide(b,x,y,uf,vf);
}

void background_state::vel(int k, double h, double x, double y, double &u, double &v) const
{
    const item &b = bg[k];
    double uf, vf;
    const double et = tide(b,x,y,uf,vf);
    
    if(b.file_uv)
    {
    u = b.r*b.U + uf;
    v = b.r*b.V + vf;
    return;
    }
    
    const double c = h>1.0e-6 ? sqrt(g/h) : 0.0;
    
    u = b.r*b.U + et*c*cos(b.dir);
    v = b.r*b.V + et*c*sin(b.dir);
}

bool background_state::carries_current(int b) const
{
    if(b<0)
    return false;
    
    const item &it = bg[b];
    
    // a tide or a time series always drives a current (long wave or file U, V)
    return it.mode==1 || it.mode==2 || it.U!=0.0 || it.V!=0.0;
}

bool background_state::carries_level(int b) const
{
    if(b<0)
    return false;
    
    const item &it = bg[b];
    
    return it.mode==1 || it.mode==2 || it.eta0!=0.0;
}

double background_state::shape(int k, double zeta, double h) const
{
    const item &b = bg[k];
    zeta = fmin(fmax(zeta,0.0),1.0);
    
    if(b.prof==1)
    {
        // ln(1 + z/z0), divided by its depth average ((h+z0) ln(1+h/z0) - h)/h
        const double z0 = b.ppar;
        const double m = ((h+z0)*log(1.0+h/z0) - h)/h;
        return m>0.0 ? log(1.0 + zeta*h/z0)/m : 1.0;
    }
    
    if(b.prof==2)
    {
        const double n = b.ppar;
        return (n+1.0)/n*pow(zeta,1.0/n);
    }
    
    return 1.0;
}
