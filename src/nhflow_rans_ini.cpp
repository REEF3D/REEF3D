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

#include"nhflow_rans_io.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include"ghostcell.h"

void nhflow_rans_io::ini(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
	/*gcval_kin=20;
	gcval_eps=30;
	gcval_edv=24;
	
	if(p->B60>=1)
	uref=p->Ui;
	
	if(p->B90>0)
	uref=0.01;
	
    if(fabs(uref)<1.0e-6)
    uref=0.01;

    plain_wallfunc(p,a,pgc);*/
    
    /*if(p->B90==1)
    LOOP
    {
    KIN[IJK] = 0.0001;
    EPS[IJK] = 1000.0001;
    }*/
    
    
    if(p->B60>=1)
    {
    LOOP
    {
    beddist = p->ZSP[IJK] - d->bed(i,j);
    
    inflow_profile(p, beddist, d->WL(i,j), KIN[IJK], EPS[IJK], d->EV[IJK]);
    
    if(p->B11==0)
    EPS[IJK] = 1.0;
    }
    
    // inflow
    inflow(p,d,pgc);
    }
    
    pgc->start20V(p,KIN,20);
    pgc->start30V(p,EPS,30);
    pgc->start24V(p,d->EV,24);
    
    // EV0 (unlimited eddy viscosity) is used for k/eps/omega diffusion and PK0: initialise it consistently
    LOOP
    d->EV0[IJK] = d->EV[IJK];
    
    pgc->start24V(p,d->EV0,24);
    
    LOOP
    if(p->DF[IJK]<0)
    {
    KIN[IJK] = 0.0;
    EPS[IJK] = 0.0;
    }
}

void nhflow_rans_io::inflow(lexer* p, fdm_nhf *d, ghostcell* pgc)
{
    double evval,kinval,epsval;
    
    if(p->B60>=1)
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
    beddist = p->ZSP[IJK] - d->bed(i,j);
    
    inflow_profile(p, beddist, d->WL(i,j), kinval, epsval, evval);

    d->EV[Im1JK] = evval;
    d->EV[Im2JK] = evval;
    d->EV[Im3JK] = evval;
    
    d->EV0[Im1JK] = evval;   // start24V keeps inflow ghosts (B60 1), so EV0 needs the profile as well
    d->EV0[Im2JK] = evval;
    d->EV0[Im3JK] = evval;
    
    KIN[Im1JK] = kinval;
    KIN[Im2JK] = kinval;
    KIN[Im3JK] = kinval;
    
    EPS[Im1JK] = epsval;
    EPS[Im2JK] = epsval;
    EPS[Im3JK] = epsval;
    }
}

// equilibrium open-channel profile, same as CFD rans_io::inflow_turb:
//   u* = Ui/(2.5 ln(11 H/ks)),  k = u*^2/sqrt(cmu) (1-z/H),  eps = u*^3/(kappa z) (1-z/H),
//   omega = u*/(sqrt(cmu) kappa z),  nu_t = kappa u* z (1-z/H);  (1-z/H) >= 0.1, z >= half the bottom cell
void nhflow_rans_io::inflow_profile(lexer* p, double z, double H, double &kv, double &ev, double &nut)
{
    const double kappa = 0.4;
    double ks = (p->S10==0) ? p->B50 : p->S20*p->S21;
    
    if(ks<=0.0)
    ks=0.0001;
    
    H = MAX(H, 1.0e-6);
    z = MAX(z, 0.5*p->DZN[KP]*H);
    
    const double ustar = fabs(p->Ui)/(2.5*log(MAX(11.0*H/ks,2.0)));
    const double fz = MAX(1.0 - z/H, 0.1);
    
    kv  = ustar*ustar/sqrt(p->cmu)*fz;
    nut = kappa*ustar*z*fz;
    
    if(p->A560==1 || p->A560==21)
    ev = pow(ustar,3.0)/(kappa*z)*fz;
    
    else
    ev = ustar/(sqrt(p->cmu)*kappa*z);
}

