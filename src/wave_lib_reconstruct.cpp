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

#include"wave_lib_reconstruct.h"
#include"wave_lib_shallow.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"wave_lib_irregular_1st.h"
#include"wave_lib_irregular_2nd_a.h"
#include"wave_lib_irregular_2nd_b.h"

wave_lib_reconstruct::wave_lib_reconstruct(lexer *p, ghostcell *pgc)
{ 
    
    if(p->mpirank==0)
    cout<<"Wave_Lib: reconstruct water waves; "<<endl;

    
    if(p->B92==51)
    ppwave = new wave_lib_irregular_1st(p,pgc);
    
    if(p->B92==52)
    ppwave = new wave_lib_irregular_2nd_a(p,pgc);
    
	if(p->B92==53)
    ppwave = new wave_lib_irregular_2nd_b(p,pgc);
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_reconstruct::~wave_lib_reconstruct()
{
}

double wave_lib_reconstruct::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = ppwave->wave_u(p,x,y,z);

    return vel;
}

double wave_lib_reconstruct::wave_v(lexer *p, double x, double y, double z)
{
    double vel;
	
    vel = ppwave->wave_v(p,x,y,z);

    return vel;
}

double wave_lib_reconstruct::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
	
    vel = ppwave->wave_w(p,x,y,z);

    return vel;
}

double wave_lib_reconstruct::wave_eta(lexer *p, double x, double y)
{
    double vel;
	
    vel = ppwave->wave_eta(p,x,y);

    return vel;
}

double wave_lib_reconstruct::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    fi = ppwave->wave_fi(p,x,y,z);
    
    return fi;
}

void wave_lib_reconstruct::parameters(lexer *p, ghostcell *pgc)
{

}

void wave_lib_reconstruct::wave_prestep(lexer *p, ghostcell *pgc)
{
}

