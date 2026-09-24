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
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"slice.h"

void wind_f::wind_forcing_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &K, slice &eta)
{
    double windforce,xc,yc,xref,yref;
    
    // A373 1/2: along-wind potential ramp
    // The ramp is measured from the start of the forcing area and the coordinates are
    // clamped to [xs,xe] x [ys,ye]: zero upstream of the area, constant downstream of it.
    // This keeps Fifsf continuous at the area edges (previously it jumped by Cd*U^2*(xs-xmin)
    // at xs and dropped to zero at xe). Without A372 the ramp is measured from the domain origin
    // as before.
    if(p->A373<3)
    {
    xref = (p->A372==1) ? xs : p->global_xmin;
    yref = (p->A372==1) ? ys : p->global_ymin;
    
    SLICELOOP4
    WETDRY
    if(p->A373==1 || eta(i,j)>0.0)
    {
    xc = MIN(MAX(p->XP[IP],xs),xe) - xref;
    yc = (p->j_dir==1) ? (MIN(MAX(p->YP[JP],ys),ye) - yref) : 0.0;
    
    windforce = (p->W3/p->W1)*Cd*p->A371_u*p->A371_u*(cosa*xc + sina*yc);
    
    K(i,j) -= windforce;
    }
    }
    
    
    // A373 3/4: Jeffreys sheltering, applied as an atmospheric surface pressure
    // p_a = rho_a * s * (U - c)^2 * d(eta)/dn  on wind-facing slopes (d(eta)/dn > 0),
    // n = wind direction, s = sheltering coefficient (A 375), c = wave phase speed.
    // Dynamic FSBC: dFi/dt = ... - p_a/rho_w
    // A373==4 restricts the forcing further to the crests (eta > 0).
    // A374==1 (with A372) applies the cos^2 downwind decay over [xs,xe] as in NHFLOW.
    if(p->A373==3 || p->A373==4)
    {
    double cph,dU,slope,psi,pa;
    
    cph = (p->A375_c>0.0) ? p->A375_c : p->wC;
    dU = p->A371_u - cph;
    
    if(dU>0.0)
    SLICELOOP4
    WETDRY
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    if(p->A373==3 || eta(i,j)>0.0)
    {
    slope = cosa*c->Ex(i,j);
    
    if(p->j_dir==1)
    slope += sina*c->Ey(i,j);
    
    if(slope>0.0)
    {
    psi = 1.0;
    
    if(p->A374==1 && p->A372==1)
    {
    psi = cos(0.5*PI*(p->XP[IP]-xs)/(xe-xs));
    psi = psi*psi;
    }
    
    pa = p->W3*p->A375_s*dU*dU*slope;
    
    K(i,j) -= psi*pa/p->W1;
    }
    }
    }
}
