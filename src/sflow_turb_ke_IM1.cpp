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

#include"sflow_turb_ke_IM1.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_idiff.h"
#include"solver2D.h"
#include"sflow_ifou.h"

#define HP (fabs(b->hp(i,j))>(p->A244)?b->hp(i,j):1.0e20)

sflow_turb_ke_IM1::sflow_turb_ke_IM1(lexer* p) : sflow_turb_io(p), kn(p), en(p), Pk(p), S(p), ustar(p), cf(p), rk(p),
                                                 wallf(p),ce1(1.44),ce2(1.92),sigk(1.0),sige(1.3),ceg(p->A264)
{
    gcval_kin=20;
	gcval_eps=30;

    // implicit first-order upwind, as the k-omega model: k and epsilon stay positive at any Courant number. The
    // WENO-HJ stencil is not monotone and takes its outer points explicitly (unstable above a Courant number of
    // about one, A 210 11-13): epsilon falls to zero in single cells, k loses its sink there and grows without bound
    pconvec = new sflow_ifou(p);
    
    pdiff = new sflow_idiff(p);
}

sflow_turb_ke_IM1::~sflow_turb_ke_IM1()
{
}

void sflow_turb_ke_IM1::start(lexer *p, fdm2D *b, ghostcell *pgc, sflow_convection *pdisc, sflow_diffusion *pdiffmom, solver2D *psolv, ioflow *pflow)
{
    Pk_update(p,b,pgc);
    ustar_update(p,b,pgc);
    fresh_cells(p,b);

    // Two passes per step over the lagged coefficients eps/k of the sinks. Both equations of a pass use the same
    // eps/k (turbulence_rate): where the sources dominate (dt eps/k >> 1, the time steps of the semi-implicit
    // momentum) k = P_k k_old/eps_old and eps = P_eps k_old/(ce2 eps_old) reach the equilibrium in two passes
    // from any state, while the eps equation with the k of the same pass (Gauss-Seidel) diverges there, each
    // pass 1.6 times further from the equilibrium than the last.
    for(int pass=0; pass<2; ++pass)
    {
    turbulence_rate(p);

    //kin
        starttime=pgc->timer();
    	clearrhs(p,b);
        pconvec->start(p,b,kin,4,b->P,b->Q);
        pdiff->diff_scalar(p,b,pgc,psolv,kin,sigk,1.0);
    	kin_source(p,b);
    	timesource(p,b,kn);
        wall_bc(p,b,kin);
        psolv->start(p,pgc,kin,b->M,b->xvec,b->rhsvec,4);
        pgc->gcsl_start4(p,kin,gcval_kin);
    	p->kintime=pgc->timer()-starttime;
    	p->kiniter=p->solveriter;
    	if(p->mpirank==0 && pass==1 && (p->count%p->P12==0))
    	cout<<"kin_iter: "<<p->kiniter<<"  kin_time: "<<setprecision(3)<<p->kintime<<endl;

    //eps
        starttime=pgc->timer();
    	clearrhs(p,b);
        pconvec->start(p,b,eps,4,b->P,b->Q);
        pdiff->diff_scalar(p,b,pgc,psolv,eps,sige,1.0);
    	eps_source(p,b);
    	timesource(p,b,en);
        wall_bc(p,b,eps);
    	psolv->start(p,pgc,eps,b->M,b->xvec,b->rhsvec,4);
        pgc->gcsl_start4(p,eps,gcval_eps);
    	p->epstime=pgc->timer()-starttime;
    	p->epsiter=p->solveriter;
    	if(p->mpirank==0 && pass==1 && (p->count%p->P12==0))
    	cout<<"eps_iter: "<<p->epsiter<<"  eps_time: "<<setprecision(3)<<p->epstime<<endl;
    }

	eddyvisc(p,b,pgc);

    SLICELOOP4
    b->kin(i,j) = kin(i,j);
    
    pgc->gcsl_start4(p,b->kin,gcval_kin);
}

void sflow_turb_ke_IM1::ktimesave(lexer* p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    kn(i,j) = kin(i,j);
}

void sflow_turb_ke_IM1::etimesave(lexer* p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    en(i,j) = eps(i,j);
}

void sflow_turb_ke_IM1::eddyvisc(lexer* p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    b->eddyv(i,j) = p->cmu*MAX(MIN(MAX(MAX(kin(i,j),0.0)*MAX(kin(i,j),0.0)
                        /((eps(i,j))>(1.0e-20)?(eps(i,j)):(1.0e20)),0.0),fabs(p->A263*kin(i,j))/(S(i,j)>1.0e-20?S(i,j):1.0e-20)),
                        0.0001*p->W2);
                        
                        
    SLICELOOP4
    if(p->wet[IJ]==0)
    {
    kin(i,j) = 0.0;
    eps(i,j) = 0.0;
    b->eddyv(i,j) = 0.0;
    }

	pgc->gcsl_start4(p,b->eddyv,24);
}

// eps/k of the sinks and of the production of epsilon, for both equations of a pass
void sflow_turb_ke_IM1::turbulence_rate(lexer* p)
{
    SLICELOOP4
    rk(i,j) = MAX(eps(i,j),0.0)/(kin(i,j)>(1.0e-10)?(fabs(kin(i,j))):(1.0e20));
}

void sflow_turb_ke_IM1::kin_source(lexer* p, fdm2D *b)
{
    count=0;
    SLICELOOP4
    {
    // dissipation implicit as (eps/k)*k, keeps k positive (as in kepsilon_func::kinsource)
    b->M.p[count] += rk(i,j);
    
	b->rhsvec.V[count]  += Pk(i,j)

                        + (1.0/sqrt(fabs(cf(i,j))>1.0e-20?cf(i,j):1.0e20))*pow(ustar(i,j),3.0)/HP;
	++count;
    }
}

void sflow_turb_ke_IM1::eps_source(lexer* p, fdm2D *b)
{
    count=0;
    SLICELOOP4
    {
    b->M.p[count] += ce2 * rk(i,j);

    b->rhsvec.V[count] += ce1 * rk(i,j)*Pk(i,j)

                       + (ceg*ce2/pow((fabs(cf(i,j))>1.0e-20?cf(i,j):1.0e20),0.75))*pow(p->cmu,0.5)*pow(ustar(i,j),4.0)/(HP*HP);
    ++count;
    }
}

void sflow_turb_ke_IM1::Pk_update(lexer* p, fdm2D *b, ghostcell *pgc)
{
    double dudx,dvdy,dudy,dvdx;
    double uc,up,um,vc,vp,vm;

    SLICELOOP4
    {
    dudx=dvdy=dudy=dvdx=0.0;

    dudx = (b->P(i,j) - b->P(i-1,j))/(p->DXM);
    dvdy = (b->Q(i,j) - b->Q(i,j-1))/(p->DXM);
    
    
    // cross derivatives: walls (free-slip in SFLOW momentum), open boundaries and dry neighbours
    // are treated as zero-gradient (mirror) instead of no-slip; no y-derivatives in 1D
    uc = 0.5*(b->P(i,j)+b->P(i-1,j));
    up = (p->flagslice4[IJp1]<0 || p->wet[IJp1]==0) ? uc : 0.5*(b->P(i,j+1)+b->P(i-1,j+1));
    um = (p->flagslice4[IJm1]<0 || p->wet[IJm1]==0) ? uc : 0.5*(b->P(i,j-1)+b->P(i-1,j-1));
    dudy = (up - um)/(2.0*p->DXM)*p->y_dir;
    
    vc = 0.5*(b->Q(i,j)+b->Q(i,j-1));
    vp = (p->flagslice4[Ip1J]<0 || p->wet[Ip1J]==0) ? vc : 0.5*(b->Q(i+1,j)+b->Q(i+1,j-1));
    vm = (p->flagslice4[Im1J]<0 || p->wet[Im1J]==0) ? vc : 0.5*(b->Q(i-1,j)+b->Q(i-1,j-1));
    dvdx = (vp - vm)/(2.0*p->DXM)*p->y_dir;
    
    dvdy *= p->y_dir;

    Pk(i,j) = b->eddyv(i,j)*(2.0*pow(dudx,2.0) + 2.0*pow(dvdy,2.0) + pow(dudy+dvdx,2.0));

    S(i,j) = sqrt(2.0*pow(dudx,2.0) + 2.0*pow(dvdy,2.0) + pow(dudy+dvdx,2.0));
    }
}

void sflow_turb_ke_IM1::ustar_update(lexer* p, fdm2D *b, ghostcell *pgc)
{
    double uvel,vvel,manning;

    SLICELOOP4
    {
    uvel=0.0;
    vvel=0.0;
    
    //if(b->wet1(i,j)==0 && b->wet1(i-1,j))
    uvel = b->U(i,j);   // cell-centred velocity, as sflow_rough_manning (face averages are halved next to dry cells)
    
    //if(b->wet2(i,j)==0 && b->wet2(i,j-1))
    vvel = b->V(i,j);

    manning = pow(b->ks(i,j),1.0/6.0)/20.0;   // same as sflow_rough_manning

    cf(i,j) = pow(manning,2.0)*fabs(p->W22)/pow(HP,1.0/3.0);

    ustar(i,j) = sqrt(cf(i,j)*(uvel*uvel + vvel*vvel));
    }

    // wallf
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
    
    SLICELOOP4
    WETDRY
    if(p->wet[Im1J]<0 || p->wet[Ip1J]<0 || p->wet[IJm1]<0 || p->wet[IJp1]<0)
    wallf(i,j)=1;
}

// Cells without turbulence from the last step (fallen wet, or the start of the run) take the local equilibrium of
// the bed sources (Rastogi and Rodi): eps = P_kv, k = u*^2/(ceg sqrt(cmu) cf^(1/4)). From k = eps = 0 the sinks
// eps/k are zero in the first step and k = dt P_kv, far above the equilibrium for the time steps of the
// semi-implicit momentum (the cell at a wetting front).
void sflow_turb_ke_IM1::fresh_cells(lexer* p, fdm2D *b)
{
    SLICELOOP4
    if(p->wet[IJ]==1 && kn(i,j)<=0.0 && en(i,j)<=0.0 && ustar(i,j)>0.0 && cf(i,j)>1.0e-20)
    {
    kin(i,j) = kn(i,j) = pow(ustar(i,j),2.0)/(ceg*sqrt(p->cmu)*pow(cf(i,j),0.25));
    eps(i,j) = en(i,j) = pow(ustar(i,j),3.0)/(HP*sqrt(cf(i,j)));
    }
}

void sflow_turb_ke_IM1::clearrhs(lexer* p, fdm2D *b)
{
    n=0;
    SLICELOOP4
    {
    b->M.p[n]  = 0.0;

    b->M.n[n] = 0.0;
    b->M.s[n] = 0.0;

    b->M.w[n] = 0.0;
    b->M.e[n] = 0.0;
        
    b->rhsvec.V[n]=0.0;
    
	b->L(i,j)=0.0;
	++n;
    }
}

