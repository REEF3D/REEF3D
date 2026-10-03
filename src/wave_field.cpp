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

#include"wave_field.h"
#include"wave_lib.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>
#include<cstdlib>
#include<iostream>
#include<string>

wave_lib* wave_lib_create(lexer*, ghostcell*, int);

wave_field::scope::scope(lexer *pp, wave_source &ss) : p(pp), s(ss)
{
    keep.save(p);
    s.ctx.load(p);
    wavetime = p->wavetime;
    p->wavetime = wavetime - s.tshift;
}

wave_field::scope::~scope()
{
    p->wavetime = wavetime;
    s.ctx.save(p);
    keep.load(p);
}

wave_field::wave_field(lexer *p, ghostcell *pgc)
{
    read(p,pgc);

    if(src.empty())
    return;

    check(p);
    build(p,pgc);
    log(p);
}

wave_field::~wave_field()
{
    for(wave_source *s : src)
    delete s;
}

bool wave_field::nonlinear(int t)
{
    // linear theories: shallow, linear, deep, 1st-order irregular
    return !(t==1 || t==2 || t==3 || t==31 || t==41 || t==51);
}

bool wave_field::exists(int k) const
{
    for(const wave_source *s : src)
    if(s->id==k)
    return true;
    
    return false;
}

bool wave_field::use(const wave_source *s) const
{
    if(filter==nullptr)
    return true;
    
    for(int k : *filter)
    if(k==s->id)
    return true;
    
    return false;
}

// ---------------------------------------------------------------------
// input
// ---------------------------------------------------------------------

void wave_field::read(lexer *p, ghostcell *pgc)
{
    // B 500-504 (control.h, read by read_control and sent to all ranks)
    auto fail = [&](const std::string &msg)
    {
        if(p->mpirank==0)
        std::cout<<std::endl<<"!!! wave_field: "<<msg<<" !!!"<<std::endl<<std::endl;
        std::exit(1);
    };
    
    auto find = [&](int k)->wave_source*
    {
        for(wave_source *s : src)
        if(s->id==k)
        return s;
        return nullptr;
    };
    
    for(int n=0; n<p->B500; ++n)
    {
        const int k = p->B500_id[n];
        
        if(k<2)
        fail("source ids start at 2 (id 1 is the B 92 wave)");
        if(find(k)!=nullptr)
        fail("source id "+std::to_string(k)+" defined twice in B 500");
        
        wave_source *s = new wave_source(k,p->B500_type[n]);
        s->H = p->B500_H[n];
        s->T = p->B500_T[n];
        src.push_back(s);
    }
    
    for(int n=0; n<p->B501; ++n)
    {
        wave_source *s = find(p->B501_id[n]);
        if(s==nullptr)
        fail("B 501 refers to source "+std::to_string(p->B501_id[n])+", which has no B 500");
        
        s->rot    = p->B501_dir[n];
        s->phase  = p->B501_phase[n];
        s->ts     = p->B501_ts[n];
        s->te     = p->B501_te[n];
        s->t_ramp = p->B501_tramp[n];
    }
    
    for(int n=0; n<p->B502; ++n)
    {
        wave_source *s = find(p->B502_id[n]);
        if(s==nullptr)
        fail("B 502 refers to source "+std::to_string(p->B502_id[n])+", which has no B 500");
        
        s->x0 = p->B502_x[n];
        s->y0 = p->B502_y[n];
    }
    
    for(int n=0; n<p->B504; ++n)
    {
        wave_source *s = find(p->B504_id[n]);
        if(s==nullptr)
        fail("B 504 refers to source "+std::to_string(p->B504_id[n])+", which has no B 500");
        
        s->seed = p->B504_seed[n];
    }
    
    for(wave_source *s : src)
    {
        const double r = s->rot*(3.14159265358979323846/180.0);
        s->cr = cos(r);
        s->sr = sin(r);
        s->tshift = s->phase/360.0*s->T;
    }
}

void wave_field::check(lexer *p)
{
    std::string err;
    int nnonlin = nonlinear(p->B92) ? 1 : 0;

    for(wave_source *s : src)
    {
        const int t = s->type;

        if(t==61)
        err = "source "+std::to_string(s->id)+": HDC coupling (61) is only available as the B 92 wave";

        else
        if(t>=20 && t<30)
        err = "source "+std::to_string(s->id)+": wavemaker theories (20-24) are only available as the B 92 wave";

        else
        if(t==0)
        err = "source "+std::to_string(s->id)+": type 0 is no wave";

        if(nonlinear(t))
        ++nnonlin;

        if(s->T<=0.0 || s->H<0.0)
        err = "source "+std::to_string(s->id)+": needs H >= 0 and T > 0";

        if(s->te<=s->ts)
        err = "source "+std::to_string(s->id)+": te must be larger than ts";
    }

    if(nnonlin>1)
    err = "at most one nonlinear wave theory may be superposed with linear ones ("+std::to_string(nnonlin)+" given)";

    if(p->B89==1)
    err = "decomposed precalc (B 89 1) supports one wave source only, until the precalc engine of phase 2";

    if(!err.empty())
    {
        if(p->mpirank==0)
        std::cout<<std::endl<<"!!! wave_field: "<<err<<" !!!"<<std::endl<<std::endl;
        std::exit(1);
    }
}

// ---------------------------------------------------------------------
// construction: every source starts from the legacy lexer state, gets its
// own wave input, constructs its wave_lib and keeps the resulting state
// ---------------------------------------------------------------------

void wave_field::build(lexer *p, ghostcell *pgc)
{
    wave_lexer_context legacy;
    legacy.save(p);

    for(wave_source *s : src)
    {
        legacy.load(p);

        p->B92   = s->type;
        p->B91   = 0;
        p->B93   = 1;
        p->B93_1 = s->H;
        p->B93_2 = s->T;
        p->B105_1 = legacy.B105_1 + s->rot;

        if(s->seed>0)
        {
            p->B139   = s->seed;
            p->B138   = 1;
            p->B138_1 = s->seed;
            p->B138_2 = s->seed+7;
        }

        p->wts = 0.0;
        p->wte = 1.0e20;

        if(p->mpirank==0)
        std::cout<<"wave_field: source "<<s->id<<std::endl;

        s->lib = wave_lib_create(p,pgc,s->type);
        s->ctx.save(p);
    }

    legacy.load(p);
}

void wave_field::log(lexer *p)
{
    if(p->mpirank!=0)
    return;

    std::cout<<"wave_field: "<<src.size()+1<<" sources, equivalent input:"<<std::endl;

    const bool irr = p->B92>30 && p->B92!=70;
    std::cout<<"  B 500 1 "<<p->B92<<" "<<(irr?p->wHs:p->wH)<<" "<<(irr?p->wTp:p->wT)
             <<"      (B 92 wave)"<<std::endl;

    for(wave_source *s : src)
    {
        std::cout<<"  B 500 "<<s->id<<" "<<s->type<<" "<<s->H<<" "<<s->T<<std::endl;
        std::cout<<"  B 501 "<<s->id<<" "<<s->rot<<" "<<s->phase<<" "<<s->ts<<" "<<s->te<<" "<<s->t_ramp<<std::endl;
        if(s->x0!=0.0 || s->y0!=0.0)
        std::cout<<"  B 502 "<<s->id<<" "<<s->x0<<" "<<s->y0<<std::endl;
        if(s->seed>0)
        std::cout<<"  B 504 "<<s->id<<" "<<s->seed<<std::endl;
    }
}

// ---------------------------------------------------------------------
// evaluation: the sum of the additional sources
// ---------------------------------------------------------------------

double wave_field::eta(lexer *p, double x, double y)
{
    double val=0.0, xs, ys;

    for(wave_source *s : src)
    if(s->active(p))
    {
        s->local(x,y,xs,ys);
        scope sc(p,*s);
        val += s->ramp(p)*s->lib->wave_eta(p,xs,ys);
    }

    return val;
}

double wave_field::u(lexer *p, double x, double y, double z)
{
    double val=0.0, xs, ys;

    for(wave_source *s : src)
    if(s->active(p))
    {
        s->local(x,y,xs,ys);
        scope sc(p,*s);
        val += s->ramp(p)*s->lib->wave_u(p,xs,ys,z);
    }

    return val;
}

double wave_field::v(lexer *p, double x, double y, double z)
{
    double val=0.0, xs, ys;

    for(wave_source *s : src)
    if(s->active(p))
    {
        s->local(x,y,xs,ys);
        scope sc(p,*s);
        val += s->ramp(p)*s->lib->wave_v(p,xs,ys,z);
    }

    return val;
}

double wave_field::w(lexer *p, double x, double y, double z)
{
    double val=0.0, xs, ys;

    for(wave_source *s : src)
    if(s->active(p))
    {
        s->local(x,y,xs,ys);
        scope sc(p,*s);
        val += s->ramp(p)*s->lib->wave_w(p,xs,ys,z);
    }

    return val;
}

double wave_field::fi(lexer *p, double x, double y, double z)
{
    double val=0.0, xs, ys;

    for(wave_source *s : src)
    if(s->active(p))
    {
        s->local(x,y,xs,ys);
        scope sc(p,*s);
        val += s->ramp(p)*s->lib->wave_fi(p,xs,ys,z);
    }

    return val;
}

void wave_field::cache_points(lexer *p, const std::vector<double> &x, const std::vector<double> &y)
{
    std::vector<double> xs(x.size()), ys(y.size());

    for(wave_source *s : src)
    {
        for(size_t q=0; q<x.size(); ++q)
        s->local(x[q],y[q],xs[q],ys[q]);

        scope sc(p,*s);
        s->lib->wave_cache_points(p,xs,ys);
    }
}

double wave_field::eta_c(lexer *p, int q)
{
    double val=0.0;

    for(wave_source *s : src)
    if(use(s) && s->active(p))
    {
        scope sc(p,*s);
        val += s->ramp(p)*s->lib->wave_eta_c(p,q);
    }

    return val;
}

double wave_field::fi_c(lexer *p, int q, double z)
{
    double val=0.0;

    for(wave_source *s : src)
    if(use(s) && s->active(p))
    {
        scope sc(p,*s);
        val += s->ramp(p)*s->lib->wave_fi_c(p,q,z);
    }

    return val;
}

void wave_field::uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
{
    double us,vs,ws;

    for(wave_source *s : src)
    if(use(s) && s->active(p))
    {
        scope sc(p,*s);
        const double r = s->ramp(p);
        us=vs=ws=0.0;
        s->lib->wave_uvw_c(p,q,z,us,vs,ws);
        u += r*us;
        v += r*vs;
        w += r*ws;
    }
}

void wave_field::prestep(lexer *p, ghostcell *pgc)
{
    for(wave_source *s : src)
    {
        scope sc(p,*s);
        s->lib->wave_prestep(p,pgc);
    }
}
