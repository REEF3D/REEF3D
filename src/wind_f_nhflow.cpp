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

void wind_f::wind_forcing_nhf_x(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *F, slice &WL, slice &eta)
{
    k=p->knoz-1;
    double psi, xpos; // wind decay
    
    if(p->A573<3)
    SLICELOOP4
    WETDRY
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    if(p->A573==1 || eta(i,j)>0.0)
        {
            psi = cos(0.5*PI*(p->XP[IP]-xs)/(xe-xs));
            psi = psi*psi;
            F[IJK] += psi*WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*cosa;
        }
    
    if(p->A573==3 || p->A573==4)
    SLICELOOP4
    WETDRY
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    {
    Sx = d->Ex(i,j);
    xpos = p->XP[IP];
    
    psi = 1.0;
    
    if(p->A574==1)
    {
    psi = cos(0.5*PI*(xpos - xs)/(xe - xs));
    psi = psi*psi;
    }
    
    
    if(p->A573==3 || eta(i,j)>0.0)
    if(Sx*cosa>0.0)
    F[IJK] += psi*WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*cosa;
    
    
    d->test2D(i,j) = 0.0;
    
    if(Sx*cosa>0.0)
    if(p->A573==3 || eta(i,j)>0.0)
    d->test2D(i,j) = psi*WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*cosa;
    }
    
    LOOP
    d->test[IJK] = d->test2D(i,j);
}

void wind_f::wind_forcing_nhf_y(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *G, slice &WL, slice &eta)
{
    k=p->knoz-1;
    
    if(p->A573<3)
    SLICELOOP4
    WETDRY
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    if(p->A573==1 || eta(i,j)>0.0)
    G[IJK] += WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*sina;
    
    if(p->A573==3 || p->A573==3)
    SLICELOOP4
    WETDRY
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    {
    Sx = d->Ex(i,j);
    
    if(p->A573==3 || eta(i,j)>0.0)
    if(Sy*sina>0.0)
    G[IJK] += WL(i,j)*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*sina;
    }
}

