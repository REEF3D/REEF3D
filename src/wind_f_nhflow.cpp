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
#include"slice4.h"

// NHFLOW wind forcing, momentum sources d(UH)/dt = F, d(VH)/dt = G
//
// A573 == 1 : surface shear stress tau = rho_a*Cd*U10^2 in the top sigma layer,
//             F += tau_x/(rho_w*dsigma_top), i.e. an acceleration tau_x/(rho_w*dz_top) in that layer
//             (same scaling as the bed shear stress in nhflow_bcmom::roughness_u).
//             Depth-integrated: d(hU)/dt = tau_x/rho_w, closed-basin setup g*h*deta/dx = tau_x/rho_w.
// A573 == 2 : deprecated (crest-masked stress), treated as A573 == 1.
// A573 == 3 : Miles/Plant wave growth via surface pressure,  A573 == 4 : modified Jeffreys,
//             see wind_f::wind_pressure (A 575, A 576, A 578). The atmospheric pressure is transmitted
//             through the water column, so it acts in every layer: F -= WL/rho_w * dp_a/dx.
// A574 == 1 (with A572): cos^2 downwind decay over [xs,xe], all modes.

void wind_f::wind_forcing_nhf_x(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *F, slice &WL, slice &eta)
{
    wind_forcing_nhf_dir(p,d,pgc,F,WL,eta,0);
}

void wind_f::wind_forcing_nhf_y(lexer *p, fdm_nhf *d, ghostcell *pgc, double *U, double *V, double *G, slice &WL, slice &eta)
{
    if(p->j_dir==1)
    wind_forcing_nhf_dir(p,d,pgc,G,WL,eta,1);
}

void wind_f::wind_forcing_nhf_dir(lexer *p, fdm_nhf *d, ghostcell *pgc, double *F, slice &WL, slice &eta, int dir)
{
    const double comp = (dir==0) ? cosa : sina;
    double psi,dpdx;
    
    // surface shear stress, top layer
    if(p->A573==1 || p->A573==2)
    {
    k = p->knoz-1;
    
    SLICELOOP4
    WETDRY
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    {
    psi = 1.0;
    
    if(p->A574==1 && p->A572==1)
    {
    psi = cos(0.5*PI*(p->XP[IP]-xs)/(xe-xs));
    psi = psi*psi;
    }
    
    F[IJK] += psi*(p->W3/p->W1)*Cd*p->A571_u*p->A571_u*comp/p->DZN[KP];
    }
    }
    
    // wave growth, surface pressure gradient over the whole column
    if(p->A573==3 || p->A573==4)
    {
    wind_pressure(p,pgc,eta,p->A573,p->A571_u,p->A575,p->A576_s,p->A576_sc,p->A576_c,p->A578,p->A574,p->A572);
    
    pgc->gcsl_start4(p,*Pa,1);
    
    LOOP
    WETDRY
    {
    if(dir==0)
    dpdx = ((*Pa)(i+1,j) - (*Pa)(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
    
    if(dir==1)
    dpdx = ((*Pa)(i,j+1) - (*Pa)(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
    
    F[IJK] -= WL(i,j)*dpdx/p->W1;
    }
    }
}
