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

#include"bedconc_VR.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedconc_VR::bedconc_VR(lexer *p)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    shields=p->S30;
    visc=p->W2;
    kappa=0.4;
    adist=3.0*d50;
    deltab=3.0*d50;
    Rstar=(rhosed-rhowat)/rhowat;
}

bedconc_VR::~bedconc_VR()
{
}

void bedconc_VR::start(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    double Ts,Tb;
    double ca,h,z1,P,R;
    
    SLICELOOP4
    s->cbn(i,j) = s->cbe(i,j);
    
    // cb* van Rijn
    k=0;
    SEDSLICELOOP
    {
    Ts = s->tau_crit(i,j);
    Tb = s->tau_i(i,j);                 // skin friction, not tau_eff

    Ti = MAX((Tb-Ts)/Ts, 0.0);

    Ds = d50*pow((Rstar*g)/(visc*visc),1.0/3.0);
    Ds = Ds>1.0e-10?Ds:1.0e10;

    h = MAX(s->waterlevel(i,j), 1.0e-6);

    // van Rijn reference height: physical, grid independent
    // evaluate van Rijn at a physical reference height, then transfer to the first cell centre.
    adist = MAX(s->ks_eff(i,j), 0.01*h);
    adist = MAX(adist, 2.0*d50);
    adist = MIN(adist, 0.5*h);

    // reference concentration at z = a
    ca = MIN( (0.015*d50*pow(Ti,1.5))/(pow(Ds,0.3)*adist), 0.05); // c_max = 0.05

    // height of the first cell centre above the bed
    if(p->A10==5)               // NHFLOW: ZP is sigma in [0,1]
    {
    k  = 0;
    z1 = p->ZP[KP]*h;
    }

    if(p->A10==6)               // CFD: ZP is physical z
    {
    k  = s->bedk(i,j);
    z1 = p->ZP[KP] - s->bedzh(i,j);
    }

    z1 = MAX(z1,adist);         // never extrapolate below the reference level
    z1 = MIN(z1,0.99*h);

    // Rouse transfer from z = a to z = z1
    P = s->ws/(kappa*MAX(s->shearvel_eff(i,j),1.0e-6));
    P = MIN(P,5.0);

    R = pow( (adist/(h-adist)) * ((h-z1)/z1), P );

    s->cbe(i,j) = ca*R;
    
    
    //cout<<"adist: "<<adist<<" ca: "<<ca<<" cbe: "<<s->cbe(i,j)<<endl;
    }
}




