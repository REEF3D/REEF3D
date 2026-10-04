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

#include"rans_io.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"bc_noflux.h"

void rans_io::ini(lexer* p, fdm*a, ghostcell* pgc)
{
	gcval_kin=20;
	gcval_eps=30;
	gcval_edv=24;
	
	uref=0.0;
	
	if(p->B60>=1)
	uref=p->Ui;
	
	if(p->B90>0)
	uref=0.01;
	
    if(fabs(uref)<1.0e-6)
    uref=0.01;
    
    //cout<<"UREF: "<<uref<<endl;
    /*
    if(p->T10==2)
    {
    
    
    LOOP
    {
    kin(i,j,k) = (2.0/3.0)*uref*0.07;
    eps(i,j,k) = 0.16*pow(kin(i,j,k),1.5)/(0.07*p->F60);
        
        
    }
    }*/

    plain_wallfunc(p,a,pgc);
}


void rans_io::plain_wallfunc(lexer* p, fdm*a, ghostcell* pgc)
{
    double hmax=-1.0e20;
    double hmin=+1.0e20;

    // water depth
    LOOP
    if(a->phi(i,j,k)>0.0)
    {
        hmin=MIN(hmin,p->pos_z());
        hmax=MAX(hmax,p->pos_z());
    }

    hmax=pgc->globalmax(hmax);
    hmin=pgc->globalmin(hmin);

    depth=hmax-hmin;

	tau_calc(a,p,hmax);


	LOOP
	{
	a->eddyv(i,j,k)=sqrt(tau)*0.11*depth;

	kin(i,j,k)=0.5*kinbed;

	if(p->T10==1 || p->T10==11 || p->T10==21)
	eps(i,j,k)=(0.09*kin(i,j,k)*kin(i,j,k))/(a->eddyv(i,j,k)+1.0e-20);

	if(p->T10==2 || p->T10==12 || p->T10==22)
	eps(i,j,k)=(kin(i,j,k))/(a->eddyv(i,j,k)+1.0e-20);

	if(p->T10==3 || p->T10==13)
	eps(i,j,k)=(kin(i,j,k))/(a->eddyv(i,j,k));
	}

	GC4LOOP
	if((p->gcb4[n][4]==21 || p->gcb4[n][4]==5) && !bc_periodic_face(p,p->gcb4[n][0],p->gcb4[n][1],p->gcb4[n][2],p->gcb4[n][3]))
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];

        kin(i,j,k)=kinbed;

        if(p->T10==1 || p->T10==11 || p->T10==21)
        {
        eps(i,j,k)=(pow(0.09,0.75)*pow(kin(i,j,k),1.5))/(0.5*0.4*p->DXM);
        a->eddyv(i,j,k) = p->cmu*kin(i,j,k)*kin(i,j,k)/eps(i,j,k);
        }

        if(p->T10==2 || p->T10==12 || p->T10==22)
        {
        eps(i,j,k)=pow(kin(i,j,k),0.5)/(0.5*0.4*p->DXM*pow(0.09,0.25));
        a->eddyv(i,j,k) = kin(i,j,k)/eps(i,j,k);
        }

        if(p->T10==3 || p->T10==13)
        {
        eps(i,j,k)=pow(kin(i,j,k),0.5)/(0.5*0.4*p->DXM*pow(0.09,0.25));
        a->eddyv(i,j,k) = kin(i,j,k)/eps(i,j,k);
        }
	}

	pgc->start4(p,kin,20);
	pgc->start4(p,eps,30);
	pgc->start4(p,a->eddyv,24);

    // k-omega diffuses k and omega with the unlimited eddy viscosity eddyv0, set before the first step
    LOOP
    eddyv0(i,j,k)=a->eddyv(i,j,k);

    pgc->start4(p,eddyv0,24);

}

void rans_io::tau_calc(fdm* a, lexer* p, double maxwdist)
{
	ks=p->B50;	
	
	H=B=depth+p->DXM;
	M=26.0/pow((ks/3.0),(1.0/6.0));
	I=pow(uref/(M*pow(H,(2.0/3.0))),2.0);
	tau=(9.81*H*I);
	kinbed = tau/sqrt(0.09);
}

// inflow turbulence for discharge inflow (B 60 >= 1), written into the inflow ghost cells i-1..i-3 before
// each k and eps/omega solve; the bc routines take these values as Dirichlet data.
// Equilibrium open-channel profile, consistent with the log law and the wall functions:
//   u* = Ui/(2.5 ln(11 H/ks))          (as ioflow_f::inflow_log)
//   k = u*^2/sqrt(cmu) (1-z/H),  eps = u*^3/(kappa z) (1-z/H),  omega = u*/(sqrt(cmu) kappa z),
//   i.e. nu_t = kappa u* z (1-z/H); (1-z/H) is bounded below by 0.1 at the free surface.
// Air cells and all other inflow types (e.g. wave generation) get zero gradient.
void rans_io::inflow_turb(lexer* p, fdm* a, ghostcell* pgc)
{
    const double kappa = 0.4;
    double hmin=+1.0e20;
    double hmax=-1.0e20;
    double H=0.0, ks, ustar=0.0, z, fz, kval, eval;
    int n,q;
    
    const bool keps = (p->T10==1 || p->T10==11 || p->T10==21);
    
    if(p->B60>=1)
    {
        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
        
            if(a->phi(i,j,k)>0.0)
            {
            hmin=MIN(hmin,p->ZN[KP]);
            hmax=MAX(hmax,p->ZN[KP1]);
            }
        }
        hmax=pgc->globalmax(hmax);
        hmin=pgc->globalmin(hmin);
        
        H = hmax-hmin;
        
        ks = (p->S10==0) ? p->B50 : p->S20*p->S21;
        
        if(ks<=0.0)
        ks=0.0001;
        
        if(H>0.0)
        ustar = fabs(p->Ui)/(2.5*log(MAX(11.0*H/ks,2.0)));
    }
    
    for(n=0;n<p->gcin_count;++n)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
        kval = kin(i,j,k);
        eval = eps(i,j,k);
        
        if(p->B60>=1 && H>0.0 && ustar>0.0 && a->phi(i,j,k)>0.0)
        {
        z  = MAX(p->ZP[KP]-hmin, 0.5*p->DZN[KP]);
        fz = MAX(1.0 - z/H, 0.1);
        
        kval = ustar*ustar/sqrt(p->cmu)*fz;
        
        if(keps)
        eval = pow(ustar,3.0)/(kappa*z)*fz;
        
        if(!keps)
        eval = ustar/(sqrt(p->cmu)*kappa*z);
        }
        
        for(q=1;q<=3;++q)
        {
        kin(i-q,j,k) = kval;
        eps(i-q,j,k) = eval;
        }
    }
}

// local water depth at the column (i,j) for the free-surface damping T 36 3: a->WL from the ioflow
// waterlevel_update when it found the interface (ioflow_f, iowave), otherwise (ioflow_v and ioflow_gravity
// leave WL = 0; WL = 1e-4 means no interface in the local column) the interface and the lowest fluid face
// are searched in the local column. Returns -1 if this rank's column has no interface (e.g. z-decomposition).
double rans_io::fsf_depth(lexer* p, fdm *a)
{
    if(a->WL(i,j)>2.0e-4)
    return a->WL(i,j);
    
    const int k0 = k;
    double zint=-1.0e20, zbed=1.0e20;
    
    KLOOP
    PCHECK
    {
        if(a->topo(i,j,k)>0.0)
        zbed = MIN(zbed, p->ZN[KP]);
        
        if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
        zint = MAX(zint, p->ZP[KP] + a->phi(i,j,k)*p->DZP[KP]/(a->phi(i,j,k)-a->phi(i,j,k+1)));
    }
    k = k0;
    
    if(zint<-1.0e19 || zbed>1.0e19 || zint<=zbed)
    return -1.0;
    
    return zint-zbed;
}
