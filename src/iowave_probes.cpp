/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

double iowave::wave_fsf(lexer *p, ghostcell *pgc, double x)
{
    double val=0.0;
    
    p->wavetime = p->simtime;
    
    val = wave_h(p,pgc,x,0.0,0.0);

    return val;
}

double iowave::wave_xvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    double val=0.0;
    
    p->wavetime = p->simtime;
    
    z -= p->phimean;
    
    val = wave_u(p,pgc,x,y,z);

    return val;
}

double iowave::wave_yvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    double val=0.0;
    
    p->wavetime = p->simtime;
    
    z -= p->phimean;
    
    val = wave_v(p,pgc,x,y,z);

    return val;
}

double iowave::wave_zvel(lexer *p, ghostcell *pgc, double x, double y, double z)
{
    double val=0.0;
    
    p->wavetime = p->simtime;
    
    z -= p->phimean;
    
    val = wave_w(p,pgc,x,y,z);

    return val;
}

