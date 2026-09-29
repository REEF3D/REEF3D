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
#include"slice4.h"

// FNPF wind forcing, added to the dynamic free surface boundary condition:
//
//   dFi/dt = ... + K_wind,   K_wind = -p_a/rho_w  (pressure modes)
//
// A373 == 1 : wind setup from the surface stress tau = rho_a*Cd*U10^2, spread over a reference
//             depth h_ref (A377). Applied as the along-wind potential  +tau/(rho_w*h_ref)*s,
//             s = along-wind coordinate clamped to the forcing area, so Fifsf stays continuous at
//             the area edges. Equilibrium: g*deta/ds = tau/(rho_w*h_ref), setup downwind.
//             Exact for constant depth; with varying depth there is no exact surface potential
//             (curl(tau/h) != 0), h_ref is the mean still water depth of the forcing area unless given.
// A373 == 2 : deprecated (crest-masked ramp), treated as A373 == 1.
// A373 == 3 : Miles/Plant wind input, zero-mean surface pressure in phase with the slope
//             p_a = beta*rho_a*u*^2*deta/dn,  u*^2 = Cd*U10^2,  n = wind direction,  beta = A375.
//             Energy growth gamma/omega = beta*(rho_a/rho_w)*(u*/c)^2 for every wave component,
//             beta ~ 32 reproduces Plant (1982), gamma/omega = 0.04*(u*/c)^2.
// A373 == 4 : modified Jeffreys (Touboul et al. 2006, Kharif et al. 2008): 
//             p_a = s*rho_a*(U10-c)^2*deta/dn, only on steep faces |deta/dn| > critical slope,
//             switched on smoothly over 0.8..1.0 of the critical slope. s, slope_c, c from A376.
// A374 == 1 (with A372): cos^2 downwind decay of the pressure forcing over [xs,xe], as in NHFLOW.
// A378 : number of 1-2-1 low-pass passes on the forcing slope (default 4).

void wind_f::wind_forcing_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &K, slice &eta)
{
    // wind setup
    if(p->A373==1 || p->A373==2)
    {
    double xref,yref,xc,yc,f;
    
    // reference depth, computed once: mean still water depth of the wet forcing area
    if(href_set==0)
    {
        if(p->A377>0.0)
        href = p->A377;
        
        else
        {
        double hsum=0.0, hcount=0.0;
        
        SLICELOOP4
        WETDRY
        if( p->XP[IP]>xs && p->XP[IP]<xe)
        if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
        {
        hsum += MAX(p->wd - c->bed(i,j), 0.0);
        hcount += 1.0;
        }
        
        hsum = pgc->globalsum(hsum);
        hcount = pgc->globalsum(hcount);
        
        href = (hcount>0.0) ? hsum/hcount : 0.0;
        }
        
        href = MAX(href, 1.0e-3);
        href_set = 1;
        
        if(p->mpirank==0)
        cout<<"FNPF wind setup: h_ref = "<<href<<" m"<<endl;
    }
    
    f = (p->W3/p->W1)*Cd*p->A371_u*p->A371_u/href;
    
    xref = (p->A372==1) ? xs : p->global_xmin;
    yref = (p->A372==1) ? ys : p->global_ymin;
    
    SLICELOOP4
    WETDRY
    {
    xc = MIN(MAX(p->XP[IP],xs),xe) - xref;
    yc = (p->j_dir==1) ? (MIN(MAX(p->YP[JP],ys),ye) - yref) : 0.0;
    
    K(i,j) += f*(cosa*xc + sina*yc);
    }
    }
    
    
    // wind input via surface pressure, dFi/dt = ... - p_a/rho_w
    if(p->A373==3 || p->A373==4)
    {
    wind_pressure(p,pgc,eta,p->A373,p->A371_u,p->A375,p->A376_s,p->A376_sc,p->A376_c,p->A378,p->A374,p->A372);
    
    SLICELOOP4
    WETDRY
    K(i,j) -= (*Pa)(i,j)/p->W1;
    }
}
