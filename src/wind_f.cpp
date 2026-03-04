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

#include"wind_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"

wind_f::wind_f(lexer *p) 
{
    
    
    xs = -1.0e10;
    xe =  1.0e10;
    ys = -1.0e10;
    ye =  1.0e10;
    
    if(p->A10==3)
    {
        wind_forcing_drag_coeff_fnpf(p);
        
        cosa = cos(p->A371_dir*(PI/180.0));
        sina = sin(p->A371_dir*(PI/180.0));
        
        if(p->A372==1)
        {
        xs = p->A372_xs;
        xe = p->A372_xe;
        ys = p->A372_ys;
        ye = p->A372_ye;
        }
    }
    
    
    if(p->A10==5)
    {
        wind_forcing_drag_coeff_nhflow(p);
        
        cosa = cos(p->A571_dir*(PI/180.0));
        sina = sin(p->A571_dir*(PI/180.0));
        
        if(p->A572==1)
        {
        xs = p->A572_xs;
        xe = p->A572_xe;
        ys = p->A572_ys;
        ye = p->A572_ye;
        }
    }
    
    
    
    Uref = 31.5;

}

wind_f::~wind_f()
{
}

void wind_f::wind_forcing_ini(lexer *p, ghostcell *pgc)
{
}


