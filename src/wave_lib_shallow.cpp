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

#include"wave_lib_shallow.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_shallow::wave_lib_shallow(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc)
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: shallow water waves "<<endl;
    cout<<"k: "<<wk<<" w: "<<ww<<" f: "<<wf<<" T: "<<wT<<" L: "<<wL<<" d: "<<wdt<<" kd: "<<wdt*wk<<" c: "<<p->wC<<endl;
    cout<<"d/gT^2: "<<wdt/(fabs(p->W22)*wT*wT)<<" H/gT^2: "<<wH/(fabs(p->W22)*wT*wT)<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_shallow::~wave_lib_shallow()
{
}

double wave_lib_shallow::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_shallow::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_shallow::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    vel = wa*sqrt(9.81/wdt) * cos(teta);

    return vel;
}

double wave_lib_shallow::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    vel = wa*sqrt(9.81/wdt)*((z+wdt)/wdt) * sin(teta);

    return vel;
}

double wave_lib_shallow::wave_eta(lexer *p, double x, double y)
{
    double eta;
	
	teta = wk*x-ww*(p->simtime) + pshift;

    eta =  wa * cos(teta);

    return eta;
}

double wave_lib_shallow::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    return fi;
}

void wave_lib_shallow::parameters(lexer *p, ghostcell *pgc)
{
    
    
}

void wave_lib_shallow::wave_prestep(lexer *p, ghostcell *pgc)
{
}
