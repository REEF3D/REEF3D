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
#include<vector>

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


    // T 46 I L: the initial field takes the prescribed inflow turbulence (k, eps/omega, nu_t) in the whole fluid
    // (T 46 switches the turbulence initialisation I 13 on, control::control_parse)
    if(p->T46_I>0.0 && p->T46_L>0.0)
    {
    const bool keps = (p->T10==1 || p->T10==11 || p->T10==21);
    const double k46 = inflow_turb_k46(p);
    const double e46 = inflow_turb_e46(p,k46);
    const double nut46 = keps ? p->cmu*k46*k46/e46 : k46/e46;
    
    if(p->mpirank==0)
    cout<<"T 46: inflow and initial turbulence  I: "<<p->T46_I<<"  L: "<<p->T46_L<<"  k: "<<k46
        <<(keps?"  eps: ":"  omega: ")<<e46<<"  nu_t: "<<nut46<<endl;
    
    LOOP
    {
    kin(i,j,k) = k46;
    eps(i,j,k) = e46;
    a->eddyv(i,j,k) = nut46;
    }
    }
    
    if(!(p->T46_I>0.0 && p->T46_L>0.0))
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
// Default: equilibrium open-channel profile, per inflow column, consistent with the log law and the wall functions
//   u* = Ui/(2.5 ln(11 H/ks)),  k = u*^2/sqrt(cmu) (1-z/H),  eps = u*^3/(kappa z) (1-z/H),  omega = u*/(sqrt(cmu) kappa z)
//   i.e. nu_t = kappa u* z (1-z/H); (1-z/H) is bounded below by 0.1 at the top.
// H is the depth of the column, from the lowest wet face to
//   - the interpolated interface with a free surface (level set F 30 > 0 or VOF F 80 > 0),
//   - the top face of the highest fluid cell in single phase (F 30 0, F 80 0: phi = 1 everywhere, no interface;
//     before, no H was found and the inflow kept zero gradient, i.e. the initial field).
// Bed and top of a column are taken over all ranks (global column j), so z-decomposed columns get the full depth.
// Columns without an interface use the inflow-wide values. ks is the bed roughness of roughness::ks_val (B 55 or
// B 50; with sediment S 21 S 20 or S 20 by S 28).
// T 46 I L (I > 0): prescribed inflow turbulence instead of the profile, uniform over the wet inflow:
//   k = 1.5 (I Ui)^2,  eps = cmu^0.75 k^1.5/L,  omega = sqrt(k)/(cmu^0.25 L)   (nu_t = cmu^0.25 sqrt(k) L)
// Air cells and the wave-only inflow (B 60 0) keep zero gradient.
void rans_io::inflow_turb(lexer* p, fdm* a, ghostcell* pgc)
{
    const double kappa = 0.4;
    double zbed_g=+1.0e20, zint_g=-1.0e20;
    double ks=0.0001, z, fz, kval, eval, H, zb, ustar;
    double k46=0.0, e46=0.0;
    int n,q;
    
    const bool keps = (p->T10==1 || p->T10==11 || p->T10==21);
    const bool fsf  = (p->F30>0 || p->F80>0);
    const bool t46  = (p->T46_I>0.0 && p->T46_L>0.0);
    const int ny = p->gknoy;
    
    // col[jg] = -(bed), col[ny+jg] = top of the water (interface or top of the fluid) of the global column jg
    std::vector<double> col;
    
    if(p->B60>=1 && t46)
    {
        k46 = inflow_turb_k46(p);
        e46 = inflow_turb_e46(p,k46);
    }
    
    if(p->B60>=1 && !t46)
    {
        col.assign(2*ny, -1.0e20);
        
        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];
        k=p->gcin[n][2];
        
            const int jg = j + p->origin_j;
            
            if(jg<0 || jg>=ny)
            continue;
            
            if(a->phi(i,j,k)>0.0)
            col[jg] = MAX(col[jg], -p->ZN[KP]);
            
            if(fsf && a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
            col[ny+jg] = MAX(col[ny+jg], p->ZP[KP] + a->phi(i,j,k)*p->DZP[KP]/(a->phi(i,j,k)-a->phi(i,j,k+1)));
            
            if(!fsf)
            col[ny+jg] = MAX(col[ny+jg], p->ZN[KP1]);
        }
        
        pgc->globalmax(col.data(), 2*ny);
        
        for(int jg=0;jg<ny;++jg)
        {
            zbed_g = MIN(zbed_g, -col[jg]);
            zint_g = MAX(zint_g, col[ny+jg]);
        }
        
        ks = (p->B55>0.0) ? p->B55 : p->B50;
        
        if(p->S10>0)
        ks = (p->S28==1) ? p->S21*p->S20 : p->S20;
        
        if(ks<=0.0)
        ks=0.0001;
    }
    
    for(n=0;n<p->gcin_count;++n)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
        kval = kin(i,j,k);
        eval = eps(i,j,k);
        
        if(p->B60>=1 && t46 && a->phi(i,j,k)>0.0)
        {
            kval = k46;
            eval = e46;
        }
        
        if(p->B60>=1 && !t46 && a->phi(i,j,k)>0.0)
        {
            const int jg = MAX(0, MIN(ny-1, j + p->origin_j));
            
            zb = -col[jg];
            H  = col[ny+jg] - zb;
            
            if(col[ny+jg]<-1.0e19 || zb>1.0e19 || H<=0.0)
            {
            zb = zbed_g;
            H  = zint_g - zbed_g;
            }
            
            if(H>0.0 && zb<1.0e19)
            {
            ustar = fabs(p->Ui)/(2.5*log(MAX(11.0*H/ks,2.0)));
            
            z  = MAX(p->ZP[KP]-zb, 0.5*p->DZN[KP]);
            fz = MAX(1.0 - z/H, 0.1);
            
            kval = ustar*ustar/sqrt(p->cmu)*fz;
            
            if(keps)
            eval = pow(ustar,3.0)/(kappa*z)*fz;
            
            if(!keps)
            eval = ustar/(sqrt(p->cmu)*kappa*z);
            }
        }
        
        for(q=1;q<=3;++q)
        {
        kin(i-q,j,k) = kval;
        eps(i-q,j,k) = eval;
        }
    }
}

// T 46 I L: k = 1.5 (I U)^2 with the discharge inflow velocity U = |Ui| (0.01 m/s if there is none)
double rans_io::inflow_turb_k46(lexer* p)
{
    const double U = (fabs(p->Ui)>1.0e-6) ? fabs(p->Ui) : 0.01;
    
    return MAX(1.5*pow(p->T46_I*U,2.0), 1.0e-12);
}

// T 46 I L: eps = cmu^0.75 k^1.5/L (k-epsilon), omega = sqrt(k)/(cmu^0.25 L) (k-omega)
double rans_io::inflow_turb_e46(lexer* p, double kval)
{
    if(p->T10==1 || p->T10==11 || p->T10==21)
    return pow(p->cmu,0.75)*pow(kval,1.5)/p->T46_L;
    
    return sqrt(kval)/(pow(p->cmu,0.25)*p->T46_L);
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
