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

#include"wave_lib_irregular_2nd_a.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_irregular_2nd_a::wave_lib_irregular_2nd_a(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    if(p->B85!=4 && p->B85!=5 && p->B85!=6 && p->B92!=52)
	{
        irregular_parameters(p);
        parameters(p,pgc);
        
        if(p->B92==32)
        {
        amplitudes_irregular(p);
        phases_irregular(p);
        pgc->bcast_double(ei,p->wN,0);
        }
        
        if(p->B92==42)
        {
        amplitudes_focused(p);
        phases_focused(p);
        }
	}
	
    if(p->B92==52)
    {
    recon_read(p,pgc);
    recon_parameters(p,pgc);
    parameters(p,pgc);
    }
    
	if(p->B85==4 || p->B85==5 || p->B85==6)
	{
	wavepackets_parameters(p);
	parameters(p,pgc);
	}
    
    print_components(p);
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: 2nd-order irregular waves A"<<endl;
    
    cout<<"Hs: "<<p->wHs<<" Tp: "<<p->wTp<<" wp: "<<p->wwp<<" cp: "<<p->wC<<endl;
    if(p->B92>40 && p->B92<50)
    cout<<"Focused Wave   xF: "<< p->B81_1 << " yF: " << p->B81_3 <<" tF: "<<p->B81_2<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
    
    // evaluation without the lexer (iowave redesign, step 2): number of components and
    // spreading switch of this wave, as constructed
    Nw = p->wN;
    B130v = p->B130;
}

wave_lib_irregular_2nd_a::~wave_lib_irregular_2nd_a()
{
}

void wave_lib_irregular_2nd_a::parameters(lexer *p, ghostcell *pgc)
{
}

void wave_lib_irregular_2nd_a::wave_prestep(lexer *p, ghostcell *pgc)
{
}

// ---------------------------------------------------------------------
// second-order theory (wave_lib_irregular_2nd_cache.h), built at the first
// evaluation, when amplitudes and phases are set
// ---------------------------------------------------------------------

void wave_lib_irregular_2nd_a::terms(lexer *p)
{
    if(tt.on)
    return;
    
    tt.build(Nw,Ai,wi,ki,cosbeta,sinbeta,wdt,9.81);
    dc.M = Nw;
}

void wave_lib_irregular_2nd_a::direct(lexer *p, double x, double y)
{
    terms(p);
    dc.phases_at(x,y,p->wavetime,ki,cosbeta,sinbeta,wi,ei);
}

void wave_lib_irregular_2nd_a::rotate(lexer *p, double &u, double &v)
{
    if(B130v==0)
    {
        u*=cosgamma;
        v*=singamma;
    }
}

double wave_lib_irregular_2nd_a::wave_eta(lexer *p, double x, double y)
{
    direct(p,x,y);
    return tt.eta(dc.C.data(),dc.S.data());
}

double wave_lib_irregular_2nd_a::wave_fi(lexer *p, double x, double y, double z)
{
    double u,v,w,f;
    direct(p,x,y);
    tt.pair_phases(dc.C.data(),dc.S.data(),dc.cp,dc.sp);
    tt.kin(dc.C.data(),dc.S.data(),dc.cp.data(),dc.sp.data(),z,4,u,v,w,f,dc.P,dc.iP,dc.Q);
    return f;
}

double wave_lib_irregular_2nd_a::wave_u(lexer *p, double x, double y, double z)
{
    double u,v,w,f;
    direct(p,x,y);
    tt.pair_phases(dc.C.data(),dc.S.data(),dc.cp,dc.sp);
    tt.kin(dc.C.data(),dc.S.data(),dc.cp.data(),dc.sp.data(),z,1,u,v,w,f,dc.P,dc.iP,dc.Q);
    rotate(p,u,v);
    return u;
}

double wave_lib_irregular_2nd_a::wave_v(lexer *p, double x, double y, double z)
{
    double u,v,w,f;
    direct(p,x,y);
    tt.pair_phases(dc.C.data(),dc.S.data(),dc.cp,dc.sp);
    tt.kin(dc.C.data(),dc.S.data(),dc.cp.data(),dc.sp.data(),z,2,u,v,w,f,dc.P,dc.iP,dc.Q);
    rotate(p,u,v);
    return v;
}

double wave_lib_irregular_2nd_a::wave_w(lexer *p, double x, double y, double z)
{
    double u,v,w,f;
    direct(p,x,y);
    tt.pair_phases(dc.C.data(),dc.S.data(),dc.cp,dc.sp);
    tt.kin(dc.C.data(),dc.S.data(),dc.cp.data(),dc.sp.data(),z,1,u,v,w,f,dc.P,dc.iP,dc.Q);
    return w;
}

// cached-point evaluation
void wave_lib_irregular_2nd_a::wave_cache_points(lexer *p, const std::vector<double> &x, const std::vector<double> &y)
{
    cache_x=x;
    cache_y=y;
    terms(p);
    cc.points(x,y,Nw,ki,cosbeta,sinbeta);
}

double wave_lib_irregular_2nd_a::wave_eta_c(lexer *p, int q)
{
    cc.time(p->wavetime,wi,ei);
    cc.phases(q);
    return tt.eta(cc.C.data(),cc.S.data());
}

double wave_lib_irregular_2nd_a::wave_fi_c(lexer *p, int q, double z)
{
    double u,v,w,f;
    cc.time(p->wavetime,wi,ei);
    cc.phases_pairs(q,tt);
    tt.kin(cc.C.data(),cc.S.data(),cc.cp.data(),cc.sp.data(),z,4,u,v,w,f,cc.P,cc.iP,cc.Q);
    return f;
}

void wave_lib_irregular_2nd_a::wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
{
    double f;
    cc.time(p->wavetime,wi,ei);
    cc.phases_pairs(q,tt);
    tt.kin(cc.C.data(),cc.S.data(),cc.cp.data(),cc.sp.data(),z,p->j_dir==1 ? 3 : 1,u,v,w,f,cc.P,cc.iP,cc.Q);
    rotate(p,u,v);
}

double wave_lib_irregular_2nd_a::wave_u_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return u;
}

double wave_lib_irregular_2nd_a::wave_v_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return v;
}

double wave_lib_irregular_2nd_a::wave_w_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return w;
}

// eta at one point for the times tv (iowave::timeseries)
void wave_lib_irregular_2nd_a::wave_eta_series(lexer *p, double x, double y, const std::vector<double> &tv, std::vector<double> &ev)
{
    terms(p);
    wave_lib_irregular_2nd_cache tc;
    tc.points(std::vector<double>(1,x),std::vector<double>(1,y),Nw,ki,cosbeta,sinbeta);
    ev.assign(tv.size(),0.0);
    
    for(size_t i=0; i<tv.size(); ++i)
    {
        tc.time(tv[i],wi,ei);
        tc.phases(0);
        ev[i] = tt.eta(tc.C.data(),tc.S.data());
    }
}
