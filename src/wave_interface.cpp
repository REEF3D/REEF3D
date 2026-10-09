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

#include"wave_interface.h"
#include"wave_lib_void.h"
#include"wave_lib_shallow.h"
#include"wave_lib_deep.h"
#include"wave_lib_linear.h"
#include"wave_lib_flap.h"
#include"wave_lib_flap_double.h"
#include"wave_lib_piston.h"
#include"wave_lib_piston_eta.h"
#include"wave_lib_flap_eta.h"
#include"wave_lib_Stokes_2nd.h"
#include"wave_lib_Stokes_5th.h"
#include"wave_lib_Stokes_5th_SH.h"
#include"wave_lib_cnoidal_shallow.h"
#include"wave_lib_cnoidal_1st.h"
#include"wave_lib_cnoidal_5th.h"
#include"wave_lib_solitary_1st.h"
#include"wave_lib_solitary_3rd.h"
#include"wave_lib_irregular_1st.h"
#include"wave_lib_irregular_2nd_a.h"
#include"wave_lib_irregular_2nd_b.h"
#include"wave_lib_hdc.h"
#include"wave_lib_ssgw.h"
#include"wave_field.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

// wave theory by B 92 code (also used for the additional sources of wave_field)
wave_lib* wave_lib_create(lexer *p, ghostcell *pgc, int wtype)
{
    wave_lib *pwave = nullptr;

    if(wtype==0)
    pwave = new wave_lib_void(p,pgc);
	
    if(wtype==1)
    pwave = new wave_lib_shallow(p,pgc);
    
    if(wtype==2)
    pwave = new wave_lib_linear(p,pgc);
    
    if(wtype==3)
    pwave = new wave_lib_deep(p,pgc);
    
    if(wtype==4)
    pwave = new wave_lib_Stokes_2nd(p,pgc);
    
    if(wtype==5)
    pwave = new wave_lib_Stokes_5th(p,pgc);
    
	if(wtype==6)
    pwave = new wave_lib_cnoidal_shallow(p,pgc);
    
	if(wtype==7)
    pwave = new wave_lib_cnoidal_1st(p,pgc);
    
	if(wtype==8)
    pwave = new wave_lib_cnoidal_5th(p,pgc);
    
	if(wtype==9)
    pwave = new wave_lib_solitary_1st(p,pgc);
    
	if(wtype==10)
    pwave = new wave_lib_solitary_3rd(p,pgc);
    
    if(wtype==11)
    pwave = new wave_lib_Stokes_5th_SH(p,pgc);
    
    if(wtype==20)
    pwave = new wave_lib_piston_eta(p,pgc);
    
	if(wtype==21)
    pwave = new wave_lib_piston(p,pgc);
    
	if(wtype==22)
    pwave = new wave_lib_flap(p,pgc);
    
    if(wtype==23)
    pwave = new wave_lib_flap_double(p,pgc);
    
    if(wtype==24)
    pwave = new wave_lib_flap_eta(p,pgc);
    
	if(wtype==31)
    pwave = new wave_lib_irregular_1st(p,pgc);
    
    if(wtype==32)
    pwave = new wave_lib_irregular_2nd_a(p,pgc);
    
	if(wtype==33)
    pwave = new wave_lib_irregular_2nd_b(p,pgc);
    
	if(wtype==41)
    pwave = new wave_lib_irregular_1st(p,pgc);
    
    if(wtype==42)
    pwave = new wave_lib_irregular_2nd_a(p,pgc);
    
	if(wtype==43)
    pwave = new wave_lib_irregular_2nd_b(p,pgc);
    
    if(wtype==51)
    pwave = new wave_lib_irregular_1st(p,pgc);
    
    if(wtype==52)
    pwave = new wave_lib_irregular_2nd_a(p,pgc);
    
	if(wtype==53)
    pwave = new wave_lib_irregular_2nd_b(p,pgc);
    
    if(wtype==61)
    pwave = new wave_lib_hdc(p,pgc);
    
    if(wtype==70)
    pwave = new wave_lib_ssgw(p,pgc);

    return pwave;
}

wave_interface::wave_interface(lexer *p, ghostcell *pgc) 
{ 
    p->wts=0.0;
    p->wte=1.0e20;
    
    wtype=p->B92;
    
    if(p->B94==0)
	 wD=p->phimean;
	
	if(p->B94==1)
	wD=p->B94_wdt;
    
    pwave = wave_lib_create(p,pgc,wtype);

    pfield = new wave_field(p,pgc);
    
    legacy_on = true;
}

wave_interface::~wave_interface()
{
    delete pfield;
}

double wave_interface::wave_u(lexer *p, ghostcell *pgc, double x, double y, double z)
{

    double uvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte && legacy_on)
    uvel = pwave->wave_u(p,x,y,z);

    if(pfield->size()>0)
    uvel += pfield->u(p,x,y,z);
	
    return uvel;
}

double wave_interface::wave_v(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    double vvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte && legacy_on)
    vvel = pwave->wave_v(p,x,y,z);

    if(pfield->size()>0)
    vvel += pfield->v(p,x,y,z);

    return vvel;
}

double wave_interface::wave_w(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    double wvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte && legacy_on)
    wvel = pwave->wave_w(p,x,y,z);

    if(pfield->size()>0)
    wvel += pfield->w(p,x,y,z);

    return wvel;
}

double wave_interface::wave_h(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    double lsv=p->phimean;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte && legacy_on)
    lsv=p->phimean + pwave->wave_eta(p,x,y);

    if(pfield->size()>0)
    lsv += pfield->eta(p,x,y);

    return lsv;
}

double wave_interface::wave_fi(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    double pval=0.0;
    
    z = MAX(z,-wD);
    
    if(legacy_on)
    pval = pwave->wave_fi(p,x,y,z);

    if(pfield->size()>0)
    pval += pfield->fi(p,x,y,z);
	
    return pval;
}

void wave_interface::wave_cache_points(lexer *p, ghostcell *pgc, const std::vector<double> &x, const std::vector<double> &y)
{
    pwave->wave_cache_points(p,x,y);
    pfield->cache_points(p,x,y);
}

double wave_interface::wave_eta_c(lexer *p, ghostcell *pgc, int q)
{
    double eta=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte && legacy_on)
    eta = pwave->wave_eta_c(p,q);

    if(pfield->size()>0)
    eta += pfield->eta_c(p,q);
	
    return eta;
}

double wave_interface::wave_fi_c(lexer *p, ghostcell *pgc, int q, double z)
{
    z = MAX(z,-wD);
    
    double val = legacy_on ? pwave->wave_fi_c(p,q,z) : 0.0;

    if(pfield->size()>0)
    val += pfield->fi_c(p,q,z);

    return val;
}

void wave_interface::wave_uvw_c(lexer *p, ghostcell *pgc, int q, double z, double &u, double &v, double &w)
{
    u=v=w=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte && legacy_on)
    pwave->wave_uvw_c(p,q,z,u,v,w);

    if(pfield->size()>0)
    pfield->uvw_c(p,q,z,u,v,w);
}

double wave_interface::wave_eta(lexer *p, ghostcell *pgc, double x, double y)
{
    double eta=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte && legacy_on)
    eta = pwave->wave_eta(p,x,y);

    if(pfield->size()>0)
    eta += pfield->eta(p,x,y);
	
    return eta;
}

double wave_interface::wave_um(lexer *p, ghostcell *pgc, double x, double y)
{
    double um;
    
    return um;
}

double wave_interface::wave_vm(lexer *p, ghostcell *pgc, double x, double y)
{
    double vm;
    
    return vm;
}

void wave_interface::wave_prestep(lexer *p, ghostcell *pgc)
{
    pwave->wave_prestep(p,pgc);
    pfield->prestep(p,pgc);
}

int wave_interface::printcheck=0;

double wave_interface::wave_paddle_Q(lexer *p, ghostcell *pgc, double z)
{
    double val=0.0;

    if(p->simtime>=p->wts && p->simtime<=p->wte)
    val = pwave->wave_paddle_Q(p,z);

    return val;
}

void wave_interface::select_sources(const std::vector<int> *ids)
{
    legacy_on = true;
    
    if(ids!=nullptr)
    {
        legacy_on = false;
        
        for(int k : *ids)
        if(k==1)
        legacy_on = true;
    }
    
    pfield->filter = ids;
}

bool wave_interface::source_exists(int k) const
{
    return k==1 || pfield->exists(k);
}

int wave_interface::wave_nsources() const
{
    return 1 + pfield->size();
}

wave_lib *wave_interface::wave_source_lib(int n, int &id, int &type, double &rot) const
{
    if(n==0)
    {
        id = 1;
        type = wtype;
        rot = 0.0;
        return pwave;
    }
    
    const wave_source *s = pfield->source(n-1);
    id = s->id;
    type = s->type;
    rot = s->rot;
    return s->lib;
}

// default: wave_eta at each time
void wave_lib::wave_eta_series(lexer *p, double x, double y, const std::vector<double> &tv, std::vector<double> &ev)
{
    const double wt = p->wavetime;
    ev.assign(tv.size(),0.0);
    
    for(size_t i=0; i<tv.size(); ++i)
    {
        p->wavetime = tv[i];
        ev[i] = wave_eta(p,x,y);
    }
    
    p->wavetime = wt;
}

// eta of the B 92 wave at one point for the times tv; false with additional sources
// (iowave::timeseries then evaluates wave_eta itself)
bool wave_interface::wave_eta_series(lexer *p, ghostcell *pgc, double x, double y, const std::vector<double> &tv, std::vector<double> &ev)
{
    if(pfield->size()>0)
    return false;
    
    ev.assign(tv.size(),0.0);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte && legacy_on)
    pwave->wave_eta_series(p,x,y,tv,ev);
    
    return true;
}

// default: u, v, w at the stored coordinates; v only in 3D (unused in 2D)
void wave_lib::wave_uvw_c(lexer *p, int q, double z, double &u, double &v, double &w)
{
    u=wave_u_c(p,q,z);
    v=p->j_dir==1 ? wave_v_c(p,q,z) : 0.0;
    w=wave_w_c(p,q,z);
}
