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

#include"nhflow_kepsilon_func.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"vrans.h"
#include"vrans_nhflow.h"

nhflow_kepsilon_func::nhflow_kepsilon_func(lexer* p, fdm_nhf *d, ghostcell *pgc) : nhflow_rans_io(p,d), nhflow_kepsilon_bc(p)
{
    sst_a1 = 0.31;
}

nhflow_kepsilon_func::~nhflow_kepsilon_func()
{
}

void  nhflow_kepsilon_func::clearfield(lexer *p, fdm_nhf *d, double *F)
{
	LOOP
	F[IJK]=0.0;
}

void nhflow_kepsilon_func::isource(lexer *p, fdm_nhf *d)
{
    if(p->T33==0)
	LOOP
	d->F[IJK]=0.0;
    
    // momentum is solved for UH: source is -(2/3) D dk/dx|z, with dk/dx|z = dk/dx|s + sigx dk/ds
    if(p->T33==1)
    LOOP
	d->F[IJK] = -(2.0/3.0)*d->WL(i,j)*((KIN[Ip1JK]-KIN[Im1JK])/(p->DXP[IP]+p->DXP[IM1])
              + 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(KIN[IJKp1]-KIN[IJKm1])/(p->DZP[KP]+p->DZP[KM1]));
}

void nhflow_kepsilon_func::jsource(lexer *p, fdm_nhf *d)
{
    if(p->T33==0)
	LOOP
	d->G[IJK]=0.0;
    
    if(p->T33==1)
    LOOP
	d->G[IJK] = -(2.0/3.0)*d->WL(i,j)*((KIN[IJp1K]-KIN[IJm1K])/(p->DYP[JP]+p->DYP[JM1])
              + 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(KIN[IJKp1]-KIN[IJKm1])/(p->DZP[KP]+p->DZP[KM1]))*p->y_dir;
}

void nhflow_kepsilon_func::ksource(lexer *p, fdm_nhf *d)
{
    if(p->T33==0)
	LOOP
	d->H[IJK]=0.0;
    
    // D * sigz * dk/ds = dk/ds
    if(p->T33==1)
    LOOP
	d->H[IJK] = -(2.0/3.0)*(KIN[IJKp1]-KIN[IJKm1])/(p->DZP[KP]+p->DZP[KM1]);
}

void nhflow_kepsilon_func::eddyvisc(lexer* p, fdm_nhf *d, ghostcell* pgc, vrans_nhflow* pvrans)
{
    // RANS (A560 1) and URANS (A560 21)
    // A564 0: nu_t = cmu k^2/eps
    // A564 1: + realizability,  nu_t <= T31 k/S   (same bound as k-omega A564 1)
    // A564 2: + SST/Bradshaw,   nu_t <= a1 k/(F2 S), with omega = eps/(cmu k)
    if(p->A560==1 || p->A560==21)
    LOOP
    {
        const double kval = MAX(KIN[IJK],0.0);
        const double Sval = strainterm(p,d);

        double ev = EPS[IJK]>1.0e-20 ? p->cmu*kval*kval/EPS[IJK] : 0.0;

        if(p->A564==1)
        ev = MIN(ev, p->T31*kval/(Sval>1.0e-20?Sval:1.0e-20));

        if(p->A564==2 && EPS[IJK]>1.0e-20)
        {
            const double wval = EPS[IJK]/(p->cmu*(kval>1.0e-20?kval:1.0e-20));

            double den = sst_a1*wval;
            const double sLim = sst_F2(p,d,kval,wval)*Sval;   // Bradshaw, blended
            const double rLim = (sst_a1/p->T31)*Sval;         // realizability, unblended
            if(sLim>den) den = sLim;
            if(rLim>den) den = rLim;

            ev = sst_a1*kval/(den>1.0e-20?den:1.0e-20);
        }

        // URANS filter
        if(p->A560==21)
        {
            if(p->j_dir==0)
            dxm = pow(p->DXN[IP]*p->DZN[KP]*d->WL(i,j), (1.0/2.0));

            if(p->j_dir==1)
            dxm = pow(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*d->WL(i,j), (1.0/3.0));

            f = MIN(1.0, dxm*p->A568*EPS[IJK]/pow((KIN[IJK]>(1.0e-20)?(KIN[IJK]):(1.0e20)),1.5));

            ev *= f;
        }

        d->EV0[IJK] = MAX(ev, 0.0001*d->VISC[IJK]);
    }

    // stabilization (Larsen & Fuhrman type): nu_t <= cmu k^2/eps * c1e Omega^2/(lambda2 c2e S^2)
    if(p->A565==0)
    LOOP
    d->EV[IJK] = d->EV0[IJK];

    if(p->A565==1)
    LOOP
    {
    double Sij2_val = Sij2(p,d);

	if(Sij2_val>1.0e-20)
    d->EV[IJK] = MIN(d->EV0[IJK], p->cmu*MAX(KIN[IJK]*KIN[IJK]
                     /((EPS[IJK])>(1.0e-20)?(EPS[IJK]):(1.0e20)),0.0)
                     *(ke_c_1e*Qij2(p,d))/(p->T42*ke_c_2e*Sij2_val));

    else
    d->EV[IJK] = d->EV0[IJK];
    }

    // VRANS: eddy viscosity treatment inside porous structures (B 295 2)
    pvrans->eddyv_func(p,d);

    LOOP
    if(p->DF[IJK]<0)
    {
    d->EV[IJK] = 0.0;
    d->EV0[IJK] = 0.0;
    }
}

void nhflow_kepsilon_func::kinsource(lexer *p, fdm_nhf *d, vrans_nhflow* pvrans)
{	
    int count=0;

    LOOP
    {
        if(WALLF[IJK]==0)
        {
        // dissipation linearised implicitly (eps/k * k): keeps k >= 0 for any dt
        d->M.p[count] += MAX(EPS[IJK],0.0)/(KIN[IJK]>1.0e-10?KIN[IJK]:1.0e-10);
        
        d->rhsvec.V[count] += PK[IJK];
        }
	++count;
    }
    
    count=0;
    if(p->A566==1)
    LOOP
    {
        d->rhsvec.V[count]  -= PK_b[IJK];
        
	++count;
    }
    
    // VRANS turbulence source (B 295 1)
    pvrans->ke_source(p,d,KIN,EPS);
}

void nhflow_kepsilon_func::epssource(lexer *p, fdm_nhf *d, vrans_nhflow* pvrans)
{
    count=0;
    
        LOOP
        {
		d->M.p[count] += ke_c_2e * MAX(EPS[IJK],0.0)/(KIN[IJK]>(1.0e-10)?(fabs(KIN[IJK])):(1.0e20));

        // c1 eps/k P(nu_t0) = c1 cmu k S^2 for nu_t0 = cmu k^2/eps; bounded by that value where the
        // 1e-4 nu floor makes nu_t0 > cmu k^2/eps, otherwise eps/k * nu_floor blows up for k -> 0
        const double ratio = MAX(EPS[IJK],0.0)/(KIN[IJK]>(1.0e-10)?(fabs(KIN[IJK])):(1.0e20));
        const double kpos  = MAX(KIN[IJK],0.0);

        if(ratio*d->EV0[IJK]<=p->cmu*kpos || d->EV0[IJK]<=1.0e-20)
        d->rhsvec.V[count] +=  ke_c_1e * ratio * PK0[IJK];
        
        else
        d->rhsvec.V[count] +=  ke_c_1e * p->cmu*kpos * PK0[IJK]/d->EV0[IJK];
        
        ++count;
        }
        
    // VRANS turbulence source (B 295 1)
    pvrans->eps_source(p,d,KIN,EPS);
}

void nhflow_kepsilon_func::epsfsf(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // free-surface value in the top sigma layer (A 567 > 0), turbulence length scale y' (Celik & Rodi 1984):
    //   A567 1: y' = T37 h (legacy k-eps meaning),  2: y' = T37,  3: y' = T37 h (h local water depth)
    //   (2 and 3 mean the same in k-eps and k-omega, as T 36 1 and 3 in CFD)
    // applied as a lower bound, so the free-surface value only ever lowers nu_t (as CFD T 36)
    k=p->knoz-1;
    
	if(p->A567>=1 && p->A567<=3)
	SLICELOOP4
	if(p->DF[IJK]>0)
	{
    double ly = (p->A567==1) ? p->T37*d->WL(i,j) : ((p->A567==2) ? p->T37 : p->T37*d->WL(i,j));
    
	EPS[IJK] = MAX(EPS[IJK], 2.5*pow(p->cmu,0.75)*pow(fabs(KIN[IJK]),1.5)/(ly>1.0e-20?ly:1.0e-20));
	}
}

double nhflow_kepsilon_func::sst_walldist(lexer *p, fdm_nhf *d)
{
    double y = p->ZSP[IJK] - d->bed(i,j);

    if(d->SOLID[IJK] > 0.0)
    y = MIN(y, d->SOLID[IJK]);

    return MAX(y, 1.0e-10);
}

double nhflow_kepsilon_func::sst_F2(lexer *p, fdm_nhf *d, double kval, double wval)
{
    double y = sst_walldist(p,d);

    kval = MAX(kval, 0.0);
    wval = MAX(wval, 1.0e-20);

    double arg2 = MAX( 2.0*sqrt(kval)/(p->cmu*wval*y),
                       500.0*d->VISC[IJK]/(y*y*wval) );

    arg2 = MIN(arg2, 25.0);

    return tanh(arg2*arg2);
}





