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

//  iowave in checkpoints (P 48, I 45): the state of the wave boundaries that carries over to the next
//  step, local to the rank (low-pass filters and edge values per column, the waves on the background,
//  the edge mass balance, the cached evaluation of the irregular waves, the samples of SIRW absorption),
//  restored on the same decomposition.  The geometry of the zones is built again from the input.

#include"iowave.h"
#include"lexer.h"
#include"wave_lib.h"
#include<iomanip>

namespace
{
struct ckw
{
    std::vector<double> &v;
    void d(double x) { v.push_back(x); }
    void vd(const std::vector<double> &a) { v.push_back(double(a.size())); v.insert(v.end(),a.begin(),a.end()); }
    void vi(const std::vector<int> &a) { v.push_back(double(a.size())); for(int x : a) v.push_back(double(x)); }
};

struct ckr
{
    const std::vector<double> &v;
    size_t k = 0;
    bool ok = true;
    double d() { if(k>=v.size()) { ok=false; return 0.0; } return v[k++]; }
    void vd(std::vector<double> &a) { const size_t n = size_t(d()); if(!ok || k+n>v.size()) { ok=false; return; } a.assign(v.begin()+k,v.begin()+k+n); k+=n; }
    void vi(std::vector<int> &a) { const size_t n = size_t(d()); if(!ok || k+n>v.size()) { ok=false; return; } a.resize(n); for(size_t m=0; m<n; ++m) a[m] = int(v[k+m]); k+=n; }
};
}

// the state of the wave libraries of all sources (cached evaluation), n = 0 the B 92 wave
static void wave_libs_state(const wave_interface &w, std::vector<double> &v)
{
    const int n = w.wave_nsources();
    v.push_back(double(n));
    for(int q=0; q<n; ++q)
    {
        int id,type; double rot;
        wave_lib *lib = w.wave_source_lib(q,id,type,rot);
        std::vector<double> s;
        if(lib!=nullptr)
        lib->ckpt_state(s);
        v.push_back(double(s.size()));
        v.insert(v.end(),s.begin(),s.end());
    }
}

static void wave_libs_restore(const wave_interface &w, const std::vector<double> &v, size_t &k)
{
    const int nv = (k<v.size()) ? int(v[k++]) : 0;
    const int n = w.wave_nsources();
    for(int q=0; q<nv; ++q)
    {
        const size_t len = (k<v.size()) ? size_t(v[k++]) : 0;
        if(k+len>v.size())
        {
            k = v.size();
            return;
        }
        if(q<n && len>0)
        {
            int id,type; double rot;
            wave_lib *lib = w.wave_source_lib(q,id,type,rot);
            std::vector<double> s(v.begin()+k,v.begin()+k+len);
            size_t j=0;
            if(lib!=nullptr)
            lib->ckpt_restore(s,j);
        }
        k += len;
    }
}

void iowave::ckpt_state(lexer *p, std::vector<double> &v)
{
    ckw w{v};
    w.d(2.0);   // version (2: SIRW)

    for(lowpass *l : {&lp_wl,&lp_u,&lp_uh,&lp_v,&lp_vh,&lp_w,&lp_wh,&lp_fifsf,
                      &lpg_eta,&lpg_fifsf,&lpg_wl,&lpg_u,&lpg_uh,&lpg_v,&lpg_vh,&lpg_w,&lpg_wh})
    {
        w.vd(l->f);
        w.vi(l->c);
    }

    for(std::vector<double> *a : {&edge_eta,&edge_us,&edge_uin,&edge_h,&edge_etaw,&edge_uw,&edge_vw,&ua_edge,
                                  &ab2_lin,&ab2_lout,&ab2_sin})
    w.vd(*a);
    w.vi(ab2_cin);
    w.vi(ab2_cout);

    w.d(double(wbg_count));
    w.d(double(wbg.size()));
    for(const wave_bg_state &s : wbg)
    {
        w.vd(s.k0); w.vd(s.k1); w.vd(s.u0); w.vd(s.u1); w.vd(s.a0); w.vd(s.a1);
        w.d(s.h0); w.d(s.h1); w.d(double(s.c0)); w.d(double(s.nblock)); w.d(s.on ? 1.0 : 0.0);
    }

    for(double x : mb_qs) w.d(x);
    for(double x : mb_int) w.d(x);
    w.d(double(mb_scount)); w.d(double(mb_stage)); w.d(mb_t); w.d(mb_t0); w.d(mb_V0); w.d(double(mb_n));

    // SIRW active absorption (B 98 8, B 99 8): the gauge samples and filter outputs, its low-pass
    sirw_in.ckpt_state(v);
    sirw_out.ckpt_state(v);
    w.vd(sw_l1);
    w.vi(sw_c1);

    // the wave libraries (cached evaluation of irregular waves)
    wave_libs_state(*this,v);
}

void iowave::ckpt_restore(lexer *p, const std::vector<double> &v)
{
    if(v.empty())
    return;

    ckr r{v};
    const int version = int(r.d());
    if(version!=1 && version!=2)
    return;

    for(lowpass *l : {&lp_wl,&lp_u,&lp_uh,&lp_v,&lp_vh,&lp_w,&lp_wh,&lp_fifsf,
                      &lpg_eta,&lpg_fifsf,&lpg_wl,&lpg_u,&lpg_uh,&lpg_v,&lpg_vh,&lpg_w,&lpg_wh})
    {
        r.vd(l->f);
        r.vi(l->c);
    }

    for(std::vector<double> *a : {&edge_eta,&edge_us,&edge_uin,&edge_h,&edge_etaw,&edge_uw,&edge_vw,&ua_edge,
                                  &ab2_lin,&ab2_lout,&ab2_sin})
    r.vd(*a);
    r.vi(ab2_cin);
    r.vi(ab2_cout);

    wbg_count = int(r.d());
    wbg.assign(size_t(r.d()),wave_bg_state());
    for(wave_bg_state &s : wbg)
    {
        r.vd(s.k0); r.vd(s.k1); r.vd(s.u0); r.vd(s.u1); r.vd(s.a0); r.vd(s.a1);
        s.h0 = r.d(); s.h1 = r.d(); s.c0 = int(r.d()); s.nblock = int(r.d()); s.on = r.d()>0.5;
    }

    for(double &x : mb_qs) x = r.d();
    for(double &x : mb_int) x = r.d();
    mb_scount = int(r.d()); mb_stage = int(r.d()); mb_t = r.d(); mb_t0 = r.d(); mb_V0 = r.d(); mb_n = int(r.d());

    if(version>=2 && r.ok)
    {
        sirw_in.ckpt_restore(v,r.k);
        sirw_out.ckpt_restore(v,r.k);
        r.vd(sw_l1);
        r.vi(sw_c1);
    }

    if(r.ok && r.k<v.size())
    wave_libs_restore(*this,v,r.k);

    // the edge mass balance (B 529) continues its file
    if(mb_t>=0.0 && p->mpirank==0 && !mb_out.is_open())
    mb_out.open(p->A10==2 ? "./REEF3D_SFLOW_Log/REEF3D-SFLOW-iowave-mass-balance.dat"
                          : "./REEF3D_Log/REEF3D-NHFLOW-iowave-mass-balance.dat",std::ios::app);

    if(!r.ok && p->mpirank==0)
    std::cout<<"I 45: the state of iowave in the checkpoint is incomplete"<<std::endl;
}
