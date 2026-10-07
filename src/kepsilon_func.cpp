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

#include"kepsilon_func.h"
#include"ioflow.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"vrans.h"

kepsilon_func::kepsilon_func(lexer* p, fdm* a, ghostcell *pgc) : rans_io(p,a), kepsilon_bc(p)
{
}

kepsilon_func::~kepsilon_func()
{
}

void  kepsilon_func::clearfield(lexer *p, fdm*  a, field& b)
{
	LOOP
	b(i,j,k)=0.0;
}

void kepsilon_func::isource(lexer *p, fdm* a)
{
    if(p->T33==0)
	ULOOP
	a->F(i,j,k)=0.0;
    
    if(p->T33==1)
    ULOOP
	a->F(i,j,k) = -(2.0/3.0)*(kin(i+1,j,k)-kin(i,j,k))/p->DXP[IP];
}

void kepsilon_func::jsource(lexer *p, fdm* a)
{
    if(p->T33==0)
	VLOOP
	a->G(i,j,k)=0.0;
    
    if(p->T33==1)
    VLOOP
	a->G(i,j,k) = -(2.0/3.0)*(kin(i,j+1,k)-kin(i,j,k))/p->DYP[JP];
}

void kepsilon_func::ksource(lexer *p, fdm* a)
{
    if(p->T33==0)
	WLOOP
	a->H(i,j,k)=0.0;
    
    if(p->T33==1)
    WLOOP
	a->H(i,j,k) = -(2.0/3.0)*(kin(i,j,k+1)-kin(i,j,k))/p->DZP[KP];
}

void  kepsilon_func::eddyvisc(fdm* a, lexer* p, ghostcell* pgc, vrans* pvrans)
{
	double H;
	double factor,epsi;
	
	LOOP
    a->eddyv(i,j,k) = MAX(MIN(p->cmu*MAX(kin(i,j,k)*kin(i,j,k)
					  /((eps(i,j,k))>(1.0e-20)?(eps(i,j,k)):(1.0e20)),0.0),fabs(p->T31*kin(i,j,k))/(strainterm(p,a)+1.0e-20)),
					  0.0001*a->visc(i,j,k));

	
	if(p->T10==21)
	LOOP
	a->eddyv(i,j,k) = MIN(a->eddyv(i,j,k), p->DXM*p->cmu*pow((kin(i,j,k)>(1.0e-20)?(kin(i,j,k)):(1.0e20)),0.5));
    
    // active wave generation / absorption: no eddy viscosity in the air next to the boundary (as komega_func)
    if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
    {
		for(int q=0;q<5;++q)
		for(int n=0;n<p->gcin_count;++n)
		{
		i=p->gcin[n][0]+q;
		j=p->gcin[n][1];
		k=p->gcin[n][2];
        
        if(i>=p->knox || p->flag4[IJK]<0)   // stay inside the local subdomain and the fluid
        continue;

		if(a->phi(i,j,k)<0.0)
		a->eddyv(i,j,k)=MIN(a->eddyv(i,j,k),1.0e-4);
		}
    }
    
    // free surface eddyv minimum (T 39 1, as komega_func): Smagorinsky value in the interface band
    if(p->T39==1)
    {
    const double c_sgs=0.2;
    
        LOOP
        {
        epsi = p->T38*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

        if(p->j_dir==0)
        epsi = p->T38*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
        
        double dirac = 0.0;
        
        if(fabs(a->phi(i,j,k))<epsi)
        dirac = 0.5*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
        
        if(dirac>0.0)
        {
        const double sgs_val = pow(c_sgs,2.0)*(p->j_dir==1?pow(p->DXN[IP]*p->DYN[JP]*p->DZN[KP],2.0/3.0):p->DXN[IP]*p->DZN[KP])
                 *strainterm(p,a->u,a->v,a->w);
                 
        a->eddyv(i,j,k) = MAX(a->eddyv(i,j,k),dirac*sgs_val);
        }
        }
    }
    
    pvrans->eddyv_func(p,a);
    
	pgc->start4(p,a->eddyv,24);
}

void  kepsilon_func::kinsource(lexer *p, fdm* a, vrans* pvrans)
{
    count=0;

    LOOP
    {
	if(wallf(i,j,k)==0)
	{
    // dissipation implicit as (eps/k) k (keeps k positive, as cmu omega k in k-omega)
    a->M.p[count] += MAX(eps(i,j,k),0.0)/(fabs(kin(i,j,k))>(1.0e-10)?(fabs(kin(i,j,k))):(1.0e20));
	a->rhsvec.V[count]  += pk(p,a,a->eddyv);
	}
	
	++count;
    }

    // buoyancy (T 45 1), G_b = -pk_b: a sink (stable stratification) is taken implicitly as (-G_b/k) k, so it
    // cannot drive k negative (Patankar); a source goes to the right-hand side (as komega_func)
    count=0;
    if(p->T45==1)
    LOOP
    {
        const double gb = -pk_b(p,a,a->eddyv);
        
        if(gb<0.0)
        a->M.p[count] += -gb/MAX(kin(i,j,k),1.0e-10);
        
        if(gb>0.0)
        a->rhsvec.V[count] += gb;
        
	++count;
    }

    pvrans->ke_source(p,a,kin);
}

void  kepsilon_func::epssource(lexer *p, fdm* a, vrans* pvrans)
{
	double dirac;
    count=0;

	LOOP
	{
    a->M.p[count] += ke_c_2e * MAX((eps(i,j,k))/(fabs(kin(i,j,k))>(1.0e-10)?(fabs(kin(i,j,k))):(1.0e20)),0.0);

    // c1 eps/k P with eps >= 0: a negative eps would turn this into a sink scaled by 1/k and grow without bound;
    // bounded by c1 cmu k S^2 where the 1e-4 nu floor makes eddyv > cmu k^2/eps (as nhflow_kepsilon_func)
    const double ratio = MAX(eps(i,j,k),0.0)/(fabs(kin(i,j,k))>(1.0e-10)?(fabs(kin(i,j,k))):(1.0e20));
    const double kpos  = MAX(kin(i,j,k),0.0);
    const double pk_ev = pk(p,a,a->eddyv);

    if(ratio*a->eddyv(i,j,k)<=p->cmu*kpos || a->eddyv(i,j,k)<=1.0e-20)
	a->rhsvec.V[count] += ke_c_1e * ratio * pk_ev;
    
    else
	a->rhsvec.V[count] += ke_c_1e * p->cmu*kpos * pk_ev/a->eddyv(i,j,k);

    ++count;
	}
    
    pvrans->eps_source(p,a,kin,eps);
}

void  kepsilon_func::epsfsf(lexer *p, fdm* a,ghostcell *pgc, ioflow *pflow)
{
	double epsi;
	
    // free-surface damping (T 36 > 0): turbulence length scale y' at the interface (Celik & Rodi 1984)
    //   eps_s = cmu^0.75 k^1.5/(kappa y'),  T36 1: y' = T37,  2: 1/y' = 1/T37 + 1/walld,  3: y' = T37 h (h local water depth)
    // applied as a lower bound eps = max(eps, w eps_s) with the dimensionless weight
    // w = 0.5(1 + cos(pi phi/epsi)) in the band |phi| < epsi (w = 1 at the interface, 0 at the band edge),
    // so the damping only ever lowers nu_t and does not depend on the grid spacing
    if(p->T36==3)
    pflow->waterlevel_update(p,a,pgc);
    
	if(p->T36>0)
	LOOP
	{
            if(p->j_dir==0)
            epsi = p->T38*(1.0/2.0)*(p->DXN[IP]+p->DZN[KP]);
            
            if(p->j_dir==1)
            epsi = p->T38*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);
            
        double w = 0.0;
        
		if(fabs(a->phi(i,j,k))<epsi)
		w = 0.5*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
	
        if(w>0.0)
        {
        double ly = p->T37;
        
        if(p->T36==2)
        ly = 1.0/(1.0/p->T37 + 1.0/(a->walld(i,j,k)>1.0e-20?a->walld(i,j,k):1.0e20));
        
        if(p->T36==3)
        {
        const double h = fsf_depth(p,a);
        ly = (h>0.0) ? p->T37*h : -1.0;   // no interface in this rank's column: no damping
        }
        
        if(ly>0.0)
        {
        ly = MAX(ly, 0.5*p->DZN[KP]);
        
        const double eps_s = 2.5*pow(p->cmu,0.75)*pow(fabs(kin(i,j,k)),1.5)/(ly>1.0e-20?ly:1.0e-20);
        
        eps(i,j,k) = MAX(eps(i,j,k), w*eps_s);
        }
        }
	}
}
