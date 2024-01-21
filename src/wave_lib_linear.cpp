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

#include"wave_lib_linear.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_linear::wave_lib_linear(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: linear waves "<<endl;
    cout<<"k: "<<wk<<" w: "<<ww<<" f: "<<wf<<" T: "<<wT<<" L: "<<wL<<" d: "<<wdt<<" kd: "<<wdt*wk<<" c: "<<p->wC<<endl;
    cout<<"d/gT^2: "<<wdt/(fabs(p->W22)*wT*wT)<<" H/gT^2: "<<wH/(fabs(p->W22)*wT*wT)<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_linear::~wave_lib_linear()
{
}

double wave_lib_linear::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_linear::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_linear::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel=0.0;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    vel = ww*wa*( cosh(wk*(wdt+z))/sinh(wk*wdt) ) * cos(teta);

    return vel;
}

double wave_lib_linear::wave_w(lexer *p, double x, double y, double z)
{
    double vel=0.0;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    vel = ww*wa*( sinh(wk*(wdt+z))/sinh(wk*wdt) ) * sin(teta);

    return vel;
}

double wave_lib_linear::wave_eta(lexer *p, double x, double y)
{
    double eta=0.0;

	teta = wk*x-ww*(p->simtime) + pshift;

    eta =  wa * cos(teta);

    return eta;
}

double wave_lib_linear::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    teta = wk*x-ww*(p->simtime) + pshift;
    
    fi = ((ww*0.5*wH)/(wk))*( cosh(wk*(wdt+z))/sinh(wk*wdt) ) * sin(teta);
    
    vel = ww*wa*( cosh(wk*(wdt+z))/sinh(wk*wdt) ) * cos(teta);
    
    return fi;
}

void wave_lib_linear::parameters(lexer *p, ghostcell *pgc)
{
}

void wave_lib_linear::wave_prestep(lexer *p, ghostcell *pgc)
{
}
