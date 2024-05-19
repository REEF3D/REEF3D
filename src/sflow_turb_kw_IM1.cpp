/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"sflow_turb_kw_IM1.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_idiff.h"
#include"solver2D.h"
#include"sflow_iweno_hj.h"


#define HP (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)

sflow_turb_kw_IM1::sflow_turb_kw_IM1(lexer* p) : sflow_turb_io(p), kn(p), wn(p), Pk(p), S(p), ustar(p), cf(p),
                                                 wallf(p), Vw(p), Qw(p),
                                                 kw_alpha(5.0/9.0), kw_beta(3.0/40.0),kw_sigma_k(2.0),kw_sigma_w(2.0)
{
    gcval_kin=20;
	gcval_eps=30;
    
    pconvec = new sflow_iweno_hj(p);
    pdiff = new sflow_idiff(p);
}

sflow_turb_kw_IM1::~sflow_turb_kw_IM1()
{
}

void sflow_turb_kw_IM1::start(lexer *p, fdm2D *b, ghostcell *pgc, sflow_convection *pdisc, sflow_diffusion *pdiffmom, solver2D *psolv, ioflow *pflow)
{
    Pk_update(p,b,pgc);
    ustar_update(p,b,pgc);

//kin
    starttime=pgc->timer();
	clearrhs(p,b);
    pconvec->start(p,b,kin,4,b->P,b->Q);
    pdiff->diff_scalar(p,b,pgc,psolv,kin,kw_sigma_k,1.0);
	kin_source(p,b);
	timesource(p,b,kn);
    wall_law_kin(p,b);
    psolv->start(p,pgc,kin,b->M,b->xvec,b->rhsvec,4);
    pgc->gcsl_start4(p,kin,gcval_kin);
	p->kintime=pgc->timer()-starttime;
	p->kiniter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"kin_iter: "<<p->kiniter<<"  kin_time: "<<setprecision(3)<<p->kintime<<endl;

//omega
    starttime=pgc->timer();
	clearrhs(p,b);
    pconvec->start(p,b,eps,4,b->P,b->Q);
    pdiff->diff_scalar(p,b,pgc,psolv,eps,kw_sigma_w,1.0);
	omega_source(p,b);
	timesource(p,b,wn);
    wall_law_omega(p,b);
	psolv->start(p,pgc,eps,b->M,b->xvec,b->rhsvec,4);
    pgc->gcsl_start4(p,eps,gcval_eps);
	p->epstime=pgc->timer()-starttime;
	p->epsiter=p->solveriter;
	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"omega_iter: "<<p->epsiter<<"  omega_time: "<<setprecision(3)<<p->epstime<<endl;

	eddyvisc(p,b,pgc);
}

void sflow_turb_kw_IM1::ktimesave(lexer* p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    kn(i,j) = kin(i,j);
}

void sflow_turb_kw_IM1::etimesave(lexer* p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    wn(i,j) = eps(i,j);
}
    
void sflow_turb_kw_IM1::eddyvisc(lexer* p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    b->eddyv(i,j) = MAX(MIN(MAX(kin(i,j)
                        /((eps(i,j))>(1.0e-20)?(eps(i,j)):(1.0e20)),0.0),fabs(p->T31*kin(i,j))/S(i,j)),
                        0.0001*p->W2);

	pgc->gcsl_start4(p,b->eddyv,24);
}

void sflow_turb_kw_IM1::kin_source(lexer* p, fdm2D *b)
{
    count=0;
    SLICELOOP4
    {
    if(wallf(i,j)==0)
    b->M.p[count] +=  p->cmu * MAX(eps(i,j),0.0);
    
    if(wallf(i,j)==0)
	b->rhsvec.V[count]  += Pk(i,j)
                    
                        + (1.0/sqrt(fabs(cf(i,j))>1.0e-20?cf(i,j):1.0e20))*pow(ustar(i,j),3.0)/HP;

	++count;
    }
}

void sflow_turb_kw_IM1::omega_source(lexer* p, fdm2D *b)
{
    count=0;
    SLICELOOP4
    {
    b->M.p[count] += kw_beta * MAX(eps(i,j),0.0);

    b->rhsvec.V[count] +=   kw_alpha * (MAX(eps(i,j),0.0)/(kin(i,j)>(1.0e-10)?(fabs(kin(i,j))):(1.0e20)))*Pk(i,k)
    
                       + (3.456/(pow((fabs(cf(i,j))>1.0e-20?cf(i,j):1.0e20),0.75))*pow(p->cmu,1.0)) * (pow(ustar(i,j),4.0)/(HP*HP));
                       
                       //+ (ceg*ce2/pow((fabs(cf(i,j))>1.0e-20?cf(i,j):1.0e20),0.75))*pow(p->cmu,0.5)*pow(ustar(i,j),4.0)/(HP*HP);
    ++count;
    }

}

void sflow_turb_kw_IM1::Pk_update(lexer* p, fdm2D *b, ghostcell *pgc)
{
    double dudx,dvdy,dudy,dvdx;
    
    SLICELOOP4
    {
    dudx=dvdy=dudy=dvdx=0.0;
    
    
    dudx = (b->P(i,j) - b->P(i-1,j))/(p->DXM);
    dvdy = (b->Q(i,j) - b->Q(i,j-1))/(p->DXM);
    dudy = (0.5*(b->P(i,j+1)+b->P(i-1,j+1)) - 0.5*(b->P(i,j-1)+b->P(i-1,j-1)))/(2.0*p->DXM);
    dvdx = (0.5*(b->Q(i+1,j)+b->Q(i+1,j-1)) - 0.5*(b->Q(i-1,j)+b->Q(i-1,j-1)))/(2.0*p->DXM);

    Pk(i,j) = b->eddyv(i,j)*(2.0*pow(dudx,2.0) + 2.0*pow(dvdy,2.0) + pow(dudy+dvdx,2.0));
    
    S(i,j) = sqrt(pow(dudx,2.0) + pow(dvdy,2.0) + 0.5*pow(dudy+dvdx,2.0));
    
    Qw(i,j) = sqrt(MAX(0.5*(dudy - dvdx)*(dvdx - dudy),0.0));


    }
}

void sflow_turb_kw_IM1::ustar_update(lexer* p, fdm2D *b, ghostcell *pgc)
{
    double uvel,vvel,manning;
    
    SLICELOOP4
    {
    uvel = 0.5*(b->P(i,j) + b->P(i-1,j));
    vvel = 0.5*(b->Q(i,j) + b->Q(i,j-1));
    
    manning = pow(b->ks(i,j),1.0/6.0)/26.0;
    
    cf(i,j) = pow(manning,2.0)*9.81/pow(HP,1.0/3.0);
    
    ustar(i,j) = sqrt(cf(i,j)*(uvel*uvel + vvel*vvel));
    }
    
    
    int n;
	SLICELOOP4
	wallf(i,j)=0;
	
	GCSL4LOOP
	if(p->gcbsl4[n][4]==21)
	{
	i = p->gcbsl4[n][0];
	j = p->gcbsl4[n][1];
	
	wallf(i,j)=1;
	}
}

void sflow_turb_kw_IM1::timesource(lexer* p, fdm2D *b, slice &fn)
{
    count=0;
    SLICELOOP4
    {
        b->M.p[count] += 1.0/p->dt;

        b->rhsvec.V[count] += b->L(i,j) + fn(i,j)/p->dt;
		
	++count;
    }
}

void sflow_turb_kw_IM1::clearrhs(lexer* p, fdm2D *b)
{
    count=0;
    SLICELOOP4
    {
    b->rhsvec.V[count]=0.0;
	b->L(i,j)=0.0;
	++count;
    }
}

// ****************************
// WALL KIN
// ****************************
void sflow_turb_kw_IM1::wall_law_kin(lexer* p, fdm2D *b)
{
    double uvel,vvel;
    double dist=0.5*p->DXM;
    double u_abs,uplus,tau,kappa;
    kappa=0.4;
    
    n=0;
	SLICELOOP4
	{
        pip=1;
        uvel=0.5*(b->P(i,j)+b->P(i-1,j));
        pip=0;

        pip=2;
        vvel=0.5*(b->Q(i,j)+b->Q(i,j-1));
        pip=0;

        u_abs = sqrt(uvel*uvel + vvel*vvel);

		if(30.0*dist<b->ks(i,j))
		dist=b->ks(i,j)/30.0;

		uplus = (1.0/kappa)*log(30.0*(dist/b->ks(i,j)));

        tau=(u_abs*u_abs)/pow((uplus>0.0?uplus:(1.0e20)),2.0);
    
    
		if(p->flagslice4[Im1J]<0 || p->flagslice4[Ip1J]<0 || p->flagslice4[IJm1]<0 || p->flagslice4[IJp1]<0)
		{
		b->M.p[n] += (pow(p->cmu,0.75)*pow(fabs(kin(i,j)),0.5)*uplus)/dist;
        b->rhsvec.V[n] += (tau*u_abs)/dist;
		}
		
	++n;
	}
    
    
    n=0;
	SLICELOOP4
	{
		if(p->flagslice4[Im1J]<0)
		{
        b->rhsvec.V[n] -= b->M.s[n]*kin(i-1,j);
		b->M.s[n] = 0.0;
		}
        
        if(p->flagslice4[Ip1J]<0)
		{
        b->rhsvec.V[n] -= b->M.n[n]*kin(i+1,j);
		b->M.n[n] = 0.0;
		}
        
        if(p->flagslice4[IJm1]<0)
		{
        b->rhsvec.V[n] -= b->M.e[n]*kin(i,j-1);
		b->M.e[n] = 0.0;
		}
        
        if(p->flagslice4[IJp1]<0)
		{
        b->rhsvec.V[n] -= b->M.w[n]*kin(i,j+1);
		b->M.w[n] = 0.0;
		}
		
	++n;
	}

}

void sflow_turb_kw_IM1::wall_law_omega(lexer* p, fdm2D *b)
{

    double dist=0.5*p->DXM;
    
    SLICELOOP4
    if(p->flagslice4[Im1J]<0 || p->flagslice4[Ip1J]<0 || p->flagslice4[IJm1]<0 || p->flagslice4[IJp1]<0)
    eps(i,j) = pow((kin(i,j)>(0.0)?(kin(i,j)):(0.0)),0.5) / (0.4*dist*pow(p->cmu, 0.25));
    
    
    n=0;
	SLICELOOP4
	{
		if(p->flagslice4[Im1J]<0)
		{
        b->rhsvec.V[n] -= b->M.s[n]*eps(i-1,j);
		b->M.s[n] = 0.0;
		}
        
        if(p->flagslice4[Ip1J]<0)
		{
        b->rhsvec.V[n] -= b->M.n[n]*eps(i+1,j);
		b->M.n[n] = 0.0;
		}
        
        if(p->flagslice4[IJm1]<0)
		{
        b->rhsvec.V[n] -= b->M.e[n]*eps(i,j-1);
		b->M.e[n] = 0.0;
		}
        
        if(p->flagslice4[IJp1]<0)
		{
        b->rhsvec.V[n] -= b->M.w[n]*eps(i,j+1);
		b->M.w[n] = 0.0;
		}
		
	++n;
	}
}












