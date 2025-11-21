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
Author: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cmath>

rheology_f::rheology_f(lexer *p) : strain(p), tau_x(p), tau_y(p), tau_z(p), epsi(p->F45*p->DXM), gravity(sqrt(p->W20*p->W20+p->W21*p->W21+p->W22*p->W22)), density_interstitial_fluid(1000.0)
{
    tanphi=0.0;
    if(p->W101>0)
    tanphi=tan(fabs(p->W102_phi)*(M_PI/180.0));
}

double rheology_f::viscosity(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    switch(p->W90)
    {
    case 1:
        val = Herschel_Bulkley(p,a,pgc,u,v,w);
        break;
    case 2:
        val = Mohr_Coulomb_and_Herschel_Bulkley(p,a,pgc);
        break;
    default:
        val=0.0;
        break;
    }

    return val;
}

double rheology_f::heaviside(int phival)
{
    double H;
    
    if(phival>epsi)
        H=1.0;
    else if(phival<-epsi)
        H=0.0;
    else
        H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));

    return H;
}

