/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"wave_lib_reconstruct.h"
#include"wave_lib_hdc.h"
#include"wave_lib_ssgw.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_interface::wave_interface(lexer *p, ghostcell *pgc) 
{ 
    p->wts=0.0;
    p->wte=1.0e20;
    
    wtype=p->B92;
    
    if(p->B94==0)
	 wD=p->phimean;
	
	if(p->B94==1)
	wD=p->B94_wdt;
    
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
}

wave_interface::~wave_interface()
{
}

double wave_interface::wave_u(lexer *p, ghostcell *pgc, double x, double y, double z)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_u(p,x,y,z);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_v(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    starttime=pgc->timer();
	
    double vvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    vvel = pwave->wave_v(p,x,y,z);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return vvel;
}

double wave_interface::wave_w(lexer *p, ghostcell *pgc, double x, double y, double z)
{
	starttime=pgc->timer();
	
    double wvel=0.0;
    
    z = MAX(z,-wD);
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    wvel = pwave->wave_w(p,x,y,z);

	p->wavetime+=pgc->timer()-starttime;

    return wvel;
}

double wave_interface::wave_h(lexer *p, ghostcell *pgc, double x, double y, double z)
{
	starttime=pgc->timer();
	
    double lsv=p->phimean;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    lsv=p->phimean + pwave->wave_eta(p,x,y);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return lsv;
}

double wave_interface::wave_fi(lexer *p, ghostcell *pgc, double x, double y, double z)
{
	starttime=pgc->timer();
	
    double pval=0.0;
    
    z = MAX(z,-wD);
    
    pval = pwave->wave_fi(p,x,y,z);
	
	p->wavetime+=pgc->timer()-starttime;

    return pval;
}

double wave_interface::wave_eta(lexer *p, ghostcell *pgc, double x, double y)
{
    starttime=pgc->timer();
	
    double eta=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    eta = pwave->wave_eta(p,x,y);
	
	p->wavetime+=pgc->timer()-starttime;
	
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
}

int wave_interface::printcheck=0;
