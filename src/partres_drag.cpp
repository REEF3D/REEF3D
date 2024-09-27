/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"partres.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"turbulence.h"

double partres::drag_model(lexer *p, double d50, double rhoS, double vel, double Ts)
{    
        Tf = 1.0-Ts;
        
        vel = fabs(vel);

        Rep = vel*d50/p->W2;

        Cd = (24.0/Rep)*(pow(Tf,-2.65) + (1.0/6.0)*pow(Rep,2.0/3.0)*pow(Tf,-1.78));
        
        Cd = MIN(Cd,10.0);
        Cd = MAX(Cd,0.0);
        
        Dp = Cd*(3.0/8.0)*(p->W1/rhoS)*(vel/(0.5*d50));
        
        return Dp;
    
}