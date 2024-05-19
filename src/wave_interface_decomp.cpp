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
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"wave_lib.h"

// U
double wave_interface::wave_u_space_sin(lexer *p, ghostcell *pgc, double x, double y, double z, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_u_space_sin(p,x,y,z,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_u_space_cos(lexer *p, ghostcell *pgc, double x, double y, double z, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_u_space_cos(p,x,y,z,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}


double wave_interface::wave_u_time_sin(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_u_time_sin(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_u_time_cos(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_u_time_cos(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

// V
double wave_interface::wave_v_space_sin(lexer *p, ghostcell *pgc, double x, double y, double z, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_v_space_sin(p,x,y,z,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_v_space_cos(lexer *p, ghostcell *pgc, double x, double y, double z, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_v_space_cos(p,x,y,z,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}


double wave_interface::wave_v_time_sin(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_v_time_sin(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_v_time_cos(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_v_time_cos(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

// W
double wave_interface::wave_w_space_sin(lexer *p, ghostcell *pgc, double x, double y, double z, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_w_space_sin(p,x,y,z,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_w_space_cos(lexer *p, ghostcell *pgc, double x, double y, double z, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_w_space_cos(p,x,y,z,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}


double wave_interface::wave_w_time_sin(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_w_time_sin(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_w_time_cos(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_w_time_cos(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

// ETA
double wave_interface::wave_eta_space_sin(lexer *p, ghostcell *pgc, double x, double y, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_eta_space_sin(p,x,y,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_eta_space_cos(lexer *p, ghostcell *pgc, double x, double y, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_eta_space_cos(p,x,y,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}


double wave_interface::wave_eta_time_sin(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_eta_time_sin(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_eta_time_cos(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_eta_time_cos(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

// FI
double wave_interface::wave_fi_space_sin(lexer *p, ghostcell *pgc, double x, double y, double z, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_fi_space_sin(p,x,y,z,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_fi_space_cos(lexer *p, ghostcell *pgc, double x, double y, double z, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_fi_space_cos(p,x,y,z,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}


double wave_interface::wave_fi_time_sin(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_fi_time_sin(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}

double wave_interface::wave_fi_time_cos(lexer *p, ghostcell *pgc, int n)
{
	starttime=pgc->timer();
	
    double uvel=0.0;
    
    if(p->simtime>=p->wts && p->simtime<=p->wte)
    uvel = pwave->wave_fi_time_cos(p,n);
	
	p->wavetime+=pgc->timer()-starttime;
	
    return uvel;
}
