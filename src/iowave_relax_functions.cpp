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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

double iowave::rb1_ext(lexer *p, int var)
{
    double x0,y0;
    
    if(var==1)
    {
    x0 = p->pos1_x();
    y0 = p->pos1_y();
    }
    
    if(var==2)
    {
    x0 = p->pos2_x();
    y0 = p->pos2_y();
    }
    
    if(var==3||var==4)
    {
    x0 = p->pos_x();
    y0 = p->pos_y();
    }

    return zones.relax_weight(x0,y0);
}

int iowave::rb1_flag(lexer *p, int var)
{
    double x0,y0;
    
    if(var==1)
    {
    x0 = p->pos1_x();
    y0 = p->pos1_y();
    }
    
    if(var==2)
    {
    x0 = p->pos2_x();
    y0 = p->pos2_y();
    }
    
    if(var==3||var==4)
    {
    x0 = p->pos_x();
    y0 = p->pos_y();
    }

    return zones.relax_flag(x0,y0);
}

double iowave::rb3_ext(lexer *p, int var)
{
    double x0,y0;
    
    if(var==1)
    {
    x0 = p->pos1_x();
    y0 = p->pos1_y();
    }
    
    if(var==2)
    {
    x0 = p->pos2_x();
    y0 = p->pos2_y();
    }
    
    if(var==3||var==4)
    {
    x0 = p->pos_x();
    y0 = p->pos_y();
    }

    return zones.beach_weight(x0,y0);
}

double iowave::ramp(lexer *p)
{
    double f=1.0;

    if(p->B101==1 && p->simtime<p->B102*p->wT)
    {
    f = p->simtime/(p->B102*p->wT) - (1.0/PI)*sin(PI*(p->simtime/(p->B102*p->wT)));
    }
    
    if(p->B101==2 && p->simtime<p->B102)
    {
    f = p->simtime/(p->B102) - (1.0/PI)*sin(PI*(p->simtime/(p->B102)));
    }

    return f;
}
