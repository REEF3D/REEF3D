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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"partres.h"
#include"particles_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

double partres::settling_velocity(lexer *p, double d, double du, double dv, double dw, double thetas) const
{
        const double dU=sqrt(du*du+dv*dv+dw*dw);
        if(dU==0) // Saveguard
        return 0;

        const double Rep=dU*d*invKinVis;

        const double Cd=24.0/Rep+4.0/sqrt(Rep)+0.4;
        const double ws_single=sqrt(4.0/3.0*(p->S22-p->W1)/p->W1*d*fabs(p->W22)/Cd);
        double n;
        if(Rep<=0.2)
        n=4.65;
        else if(Rep<=1)
        n=4.35*pow(Rep,-0.03);
        else if(Rep<=500)
        n=4.45*pow(Rep,-0.1);
        else
        n=2.39;
        const double ws_swarm = ws_single*pow((1.0-thetas),n);
        return ws_swarm;
}
