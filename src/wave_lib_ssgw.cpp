/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"wave_lib_ssgw.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_ssgw::wave_lib_ssgw(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc) 
{
}

wave_lib_ssgw::~wave_lib_ssgw()
{
}

double wave_lib_ssgw::wave_eta(lexer *p, double x, double y)
{
    return eta;
}

double wave_lib_ssgw::wave_fi(lexer *p, double x, double y, double z)
{
    return fi;
}

double wave_lib_ssgw::wave_u(lexer *p, double x, double y, double z)
{
    double vel = 0.0;
    cout<<"Not implemented yet"<<endl;
    return cosgamma*vel;
}

double wave_lib_ssgw::wave_v(lexer *p, double x, double y, double z)
{
    double vel = 0.0;
    cout<<"Not implemented yet"<<endl;
    return singamma*vel;
}

double wave_lib_ssgw::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel = 0.0;
    cout<<"Not implemented yet"<<endl;
    return vel;
}

double wave_lib_ssgw::wave_w(lexer *p, double x, double y, double z)
{
    double vel = 0.0;
    cout<<"Not implemented yet"<<endl;
    return vel;
}

void wave_lib_ssgw::parameters(lexer *p, ghostcell *pgc){}

void wave_lib_ssgw::wave_prestep(lexer *p, ghostcell *pgc){}
