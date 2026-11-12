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

    // NHFLOW (ZP is sigma in [0,1]): the equilibrium concentration of the bottom cell is the mean of the Rouse
    // profile over the cell, c = c_a below z = a. The cell value of the suspended solve is a cell mean, and the
    // profile is convex: its value at the cell centre is far below its mean in a thick bottom cell. rbe = mean/c_a
    // goes to the suspended solve, which exchanges ws (c_a - c_0/rbe) with the bed: pick-up ws c_a at the
    // reference level and deposition with the near-bed concentration c_0/rbe, on any vertical grid
    // (nhflow_suspended_IM1::suspsource).
    if(p->A10==5)
    {
    k  = 0;
    
    P = s->ws/(kappa*MAX(s->shearvel_eff(i,j),1.0e-6));
    P = MIN(P,5.0);
    
    const double zt = MIN(p->DZN[KP]*h, h);
    
    s->rbe(i,j) = rouse_mean(adist, zt, h, P);
    s->cbe(i,j) = ca*s->rbe(i,j);
    
    // the node of the bottom cell for the flux to the cell above: the cell centre, not below z = a
    z1 = MIN(MAX(p->ZP[KP]*h, adist), zt);
    
    s->rbz(i,j) = z1;
    s->rbf(i,j) = ((z1>adist && h>adist) ? pow((adist/(h-adist))*((h-z1)/z1), P) : 1.0)/s->rbe(i,j);
    continue;
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
    s->rbe(i,j) = 1.0;
    
    
    //cout<<"adist: "<<adist<<" ca: "<<ca<<" cbe: "<<s->cbe(i,j)<<endl;
    }
}

// mean over the bottom cell 0 < z < zt of the Rouse profile relative to c_a,
//     r(z) = ( a/(h-a) (h-z)/z )^P  for z > a,   r = 1 below a,
// Simpson rule in ln z (r falls like z^-P above a)
double bedconc_VR::rouse_mean(double a, double zt, double h, double P)
{
    if(zt<=a || h<=a)
    return 1.0;
    
    const int N = 32;
    const double x0 = log(a);
    const double dx = (log(zt) - x0)/double(N);
    double sum = 0.0;
    
    for(int q=0; q<=N; ++q)
    {
    const double z = MIN(exp(x0 + double(q)*dx), h);
    const double r = pow((a/(h-a))*((h-z)/z), P)*z;
    
    sum += ((q==0 || q==N) ? 1.0 : ((q%2==1) ? 4.0 : 2.0))*r;
    }
    
    const double R = (a + sum*dx/3.0)/zt;
    
    return MAX(MIN(R,1.0),1.0e-8);
}

