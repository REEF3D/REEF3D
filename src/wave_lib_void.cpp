/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"wave_lib_void.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_void::wave_lib_void(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc)
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave Tank: no wave specified; "<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_void::~wave_lib_void()
{
}

// U -------------------------------------------------------------
double wave_lib_void::wave_u(lexer *p, double x, double y, double z)
{
    double vel=0.0;

    return cosgamma*vel;
}

// V -------------------------------------------------------------
double wave_lib_void::wave_v(lexer *p, double x, double y, double z)
{
    double vel=0.0;

    return singamma*vel;
}

// W -------------------------------------------------------------
double wave_lib_void::wave_w(lexer *p, double x, double y, double z)
{
    double vel=0.0;

    return vel;
}

// ETA -------------------------------------------------------------
double wave_lib_void::wave_eta(lexer *p, double x, double y)
{
    double eta=0.0;

    return eta;
}

// FI -------------------------------------------------------------
double wave_lib_void::wave_fi(lexer *p, double x, double y, double z)
{
    double fi=0.0;
    
    return fi;
}

//  -------------------------------------------------------------
void wave_lib_void::parameters(lexer *p, ghostcell *pgc)
{
}

void wave_lib_void::wave_prestep(lexer *p, ghostcell *pgc)
{
}
