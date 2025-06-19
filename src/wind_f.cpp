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

#include"wind_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"

wind_f::wind_f(lexer *p) 
{
    wind_forcing_drag_coeff(p);
    
    cosa = cos(p->A571_dir*(PI/180.0));
    sina = sin(p->A571_dir*(PI/180.0));

}

wind_f::~wind_f()
{
}

void wind_f::wind_forcing_ini(lexer *p, ghostcell *pgc)
{
    if(p->A570==2)
    {
    if(p->A571_u<7.5)
    Cd = 1.2875e-3;
    
    if(p->A571_u>=7.5)
    Cd = (0.8 + 0.065*p->A571_u)*1.0e-3;
    }
}

void wind_f::wind_forcing_nhf_x(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *F, slice &WL, slice &eta)
{
    k=p->knoz-1;
    
    SLICELOOP4
    WETDRY
    {
    F[IJK] += WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*cosa;
    
    double drag = WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*cosa;
    
    //if(drag!=drag)
    //cout<<"DRAGL: "<<drag<<" "<<WL(i,j)<<endl;
    }
}

void wind_f::wind_forcing_nhf_y(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *G, slice &WL, slice &eta)
{
    k=p->knoz-1;
    
    SLICELOOP4
    WETDRY
    G[IJK] += WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*sina;
}

