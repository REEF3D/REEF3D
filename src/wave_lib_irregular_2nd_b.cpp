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

#include"wave_lib_irregular_2nd_b.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_irregular_2nd_b::wave_lib_irregular_2nd_b(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    if(p->B85!=4 && p->B85!=5 && p->B85!=6 && p->B92!=53)
	{
        irregular_parameters(p);
        parameters(p,pgc);
        
        if(p->B92==33)
        {
        amplitudes_irregular(p);
        phases_irregular(p);
        pgc->bcast_double(ei,p->wN,0);
        }
        
        if(p->B92==43)
        {
        amplitudes_focused(p);
        phases_focused(p);
        }
	}
    
    if(p->B92==53)
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
    cout<<"Wave_Lib: 2nd-order irregular waves B"<<endl;
    cout<<"Hs: "<<p->wHs<<" Tp: "<<p->wTp<<" wp: "<<p->wwp<<" cp: "<<p->wC<<endl;
    if(p->B92>40 && p->B92<50)
    cout<<"Focused Wave   xF: "<<p->B81_1 << " yF: " << p->B81_3 <<" tF: "<<p->B81_2<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));

    
    // evaluation without the lexer (iowave redesign, step 2): number of components and
    // spreading switch of this wave, as constructed
    Nw = p->wN;
    B130v = p->B130;
    
    // reconstructed waves with directions: evaluated as a directional sea
    if(recon_dir)
    B130v = 1;
}

wave_lib_irregular_2nd_b::~wave_lib_irregular_2nd_b()
{
}

void wave_lib_irregular_2nd_b::parameters(lexer *p, ghostcell *pgc)
{
}

void wave_lib_irregular_2nd_b::wave_prestep(lexer *p, ghostcell *pgc)
{
}

// ---------------------------------------------------------------------
// second-order theory (wave_lib_irregular_2nd_cache.h), built at the first
// evaluation, when amplitudes and phases are set
// ---------------------------------------------------------------------

void wave_lib_irregular_2nd_b::terms(lexer *p)
{
    if(tt.on)
    return;
    
    tt.build(Nw,Ai,wi,ki,cosbeta,sinbeta,wdt,9.81);
    dc.M = Nw;
}

void wave_lib_irregular_2nd_b::direct(lexer *p, double x, double y)
{
    terms(p);
    dc.phases_at(x,y,p->wavetime,ki,cosbeta,sinbeta,wi,ei);
}

void wave_lib_irregular_2nd_b::rotate(lexer *p, double &u, double &v)
{
    if(B130v==0)
    {
        u*=cosgamma;
        v*=singamma;
    }
}

double wave_lib_irregular_2nd_b::wave_eta(lexer *p, double x, double y)
{
    direct(p,x,y);
    return tt.eta(dc);
}

double wave_lib_irregular_2nd_b::wave_fi(lexer *p, double x, double y, double z)
{
    double u,v,w,f;
    direct(p,x,y);
    tt.pair_phases(dc);
    tt.kin(dc,z,4,u,v,w,f);
    return f;
}

double wave_lib_irregular_2nd_b::wave_u(lexer *p, double x, double y, double z)
{
    double u,v,w,f;
    direct(p,x,y);
    tt.pair_phases(dc);
    tt.kin(dc,z,1,u,v,w,f);
    rotate(p,u,v);
    return u;
}

double wave_lib_irregular_2nd_b::wave_v(lexer *p, double x, double y, double z)
{
    double u,v,w,f;
    direct(p,x,y);
    tt.pair_phases(dc);
    tt.kin(dc,z,2,u,v,w,f);
    rotate(p,u,v);
    return v;
}

double wave_lib_irregular_2nd_b::wave_w(lexer *p, double x, double y, double z)
{
    double u,v,w,f;
    direct(p,x,y);
    tt.pair_phases(dc);
    tt.kin(dc,z,1,u,v,w,f);
    return w;
}

// cached-point evaluation
void wave_lib_irregular_2nd_b::wave_cache_points(lexer *p, const std::vector<double> &x, const std::vector<double> &y)
{
    cache_x=x;
    cache_y=y;
    terms(p);
    cc.points(x,y,Nw,ki,cosbeta,sinbeta);
    
    // points with the same phases (the cells across the crests of a unidirectional sea); the profiles
    // exp(+-K z) of the terms have K <= 2 kmax
    zg.build(x,y,Nw,cosbeta,sinbeta,2.0*tt.kmax,wdt);
}

// cached points: eta and the potential once per group of points with the same phases, u, v, w
// from the tables in z of the group where it asked for enough values (wave_lib_zgroup.h)
double wave_lib_irregular_2nd_b::wave_eta_c(lexer *p, int q)
{
    cc.time(p->wavetime,wi,ei);
    
    return zg.eta(q,cc.t,[&]()
    {
        cc.phases(q);
        return tt.eta(cc);
    });
}

double wave_lib_irregular_2nd_b::wave_fi_c(lexer *p, int q, double z)
{
    cc.time(p->wavetime,wi,ei);
    
    return zg.fi(q,cc.t,z,[&](double zz)
    {
        double u,v,w,f;
        cc.phases_pairs(q,tt);
        tt.kin(cc,zz,4,u,v,w,f);
        return f;
    });
}

void wave_lib_irregular_2nd_b::wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
{
    cc.time(p->wavetime,wi,ei);
    
    const int mode = p->j_dir==1 ? 3 : 1;
    auto ev = [&](double zz, double &uu, double &vv, double &ww)
    {
        double f;
        cc.phases_pairs(q,tt);
        tt.kin(cc,zz,mode,uu,vv,ww,f);
    };
    
    if(!zg.uvw(q,cc.t,z,mode,ev,u,v,w))
    ev(z,u,v,w);
    
    rotate(p,u,v);
}

double wave_lib_irregular_2nd_b::wave_u_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return u;
}

double wave_lib_irregular_2nd_b::wave_v_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return v;
}

double wave_lib_irregular_2nd_b::wave_w_c(lexer *p, int q, double z)
{
    double u,v,w;
    wave_uvw_c(p,q,z,u,v,w);
    return w;
}

// eta at one point for the times tv (iowave::timeseries)
void wave_lib_irregular_2nd_b::wave_eta_series(lexer *p, double x, double y, const std::vector<double> &tv, std::vector<double> &ev)
{
    terms(p);
    
    std::vector<double> phi(Nw);
    for(int n=0; n<Nw; ++n)
    phi[n] = ki[n]*(cosbeta[n]*x + sinbeta[n]*y) - ei[n];
    
    tt.series(phi.data(),tv,ev);
}
