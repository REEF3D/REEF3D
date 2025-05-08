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

#include"density_df.h"
#include"lexer.h"
#include"fdm.h"

density_df::density_df(lexer* p)
{
    if(p->j_dir==0)
        fbpsi = 0.5*(1.0/2.0)*(p->DRM+p->DTM);
    else
        fbpsi = 0.5*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
}

double density_df::roface(lexer *p, fdm *a, int aa, int bb, int cc)
{
    double factor = (p->X46==1 && (0.5*(a->fb(i,j,k) + a->fb(i+aa,j+bb,k+cc)) < -fbpsi)) ? 2.0 : 1.0;
    
    double phival = 0.5*(a->phi(i,j,k) + a->phi(i+aa,j+bb,k+cc));

    double H;
    if(phival>factor*p->psi)
        H=1.0;
    else if(phival<-factor*p->psi)
        H=0.0;
    else
        H=0.5*(1.0 + phival/(factor*p->psi) + (1.0/PI)*sin((PI*phival)/(factor*p->psi)));

    return p->W1*H + p->W3*(1.0-H);
}
