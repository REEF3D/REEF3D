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

#include"wave_lib_cnoidal_5th.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_cnoidal_5th::wave_lib_cnoidal_5th(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    parameters(p,pgc);
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_cnoidal_5th::~wave_lib_cnoidal_5th()
{
}

double wave_lib_cnoidal_5th::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);
	
    return cosgamma*vel;
}

double wave_lib_cnoidal_5th::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);
	
    return singamma*vel;
}

double wave_lib_cnoidal_5th::wave_horzvel(lexer *p, double x, double y, double z)
{
    return model.u(x,z,p->wavetime,pshift);
}

double wave_lib_cnoidal_5th::wave_w(lexer *p, double x, double y, double z)
{
    return model.w(x,z,p->wavetime,pshift);
}

double wave_lib_cnoidal_5th::wave_eta(lexer *p, double x, double y)
{
    return model.eta(x,p->wavetime,pshift);
}

double wave_lib_cnoidal_5th::wave_fi(lexer *p, double x, double y, double z)
{
    // the theory gives no potential: the antiderivative of the horizontal velocity (wave_lib_steady_fi.h)
    if(!sfi.built())
    sfi.build(wL, wL/wT, wdt, fabs(p->W22), p->wavetime,
              [&](double s, double zz) {return wave_horzvel(p,s,0.0,zz);},
              [&](double s, double zz) {return wave_w(p,s,0.0,zz);},
              [&](double s) {return wave_eta(p,s,0.0);});
    
    return sfi.fi(x,z,p->wavetime);
}

void wave_lib_cnoidal_5th::parameters(lexer *p, ghostcell *pgc)
{
    // self-consistent cnoidal wave: h from mean level = d, m from T (B 93) or L (B 91), zero mean current
    const bool period = (p->B93==1);
    
    if(p->mpirank==0)
    cout<<"Ursell number: "<<(wH*wL*wL)/pow(wdt,3.0)<<"  (linear L)"<<endl;
    
    try
    {
    model = cnoidal_consistent::Model::solve(wdt, wH, fabs(p->W22), period?wT:wL, period);
    }
    catch(const std::exception &e)
    {
        if(p->mpirank==0)
        cout<<"Wave_Lib cnoidal 5th: "<<e.what()<<endl;
        
    exit(1);
    }
    
    wL = model.L;
    wT = model.T;
    wC = model.c;
    wk = 2.0*PI/wL;
    ww = 2.0*PI/wT;
    wf = 1.0/wT;
    
    modulus = model.m;
    Km = model.K;
    Em = model.E;
    ell = Em/Km;
    
    p->wL = wL;
    p->wT = wT;
    p->wC = wC;
    p->wk = wk;
    p->ww = ww;
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: cnoidal waves, order "<<model.order<<endl;
    
    if(model.order==3)
    cout<<"Wave_Lib: m < 0.96, third-order cnoidal with full m dependence (fifth-order Iwagaki coefficients need m >= 0.96)"<<endl;
    
    cout<<"MODULUS: "<<modulus<<"  WAVE TROUGH: "<<model.h<<endl;
    cout<<"k: "<<wk<<" w: "<<ww<<" f: "<<wf<<" T: "<<wT<<" L: "<<wL<<" d: "<<wdt<<" kd: "<<wdt*wk<<" c: "<<wC<<endl;
    cout<<"d/gT^2: "<<wdt/(fabs(p->W22)*wT*wT)<<" H/gT^2: "<<wH/(fabs(p->W22)*wT*wT)<<endl;
    }
}

void wave_lib_cnoidal_5th::wave_prestep(lexer *p, ghostcell *pgc)
{
}
