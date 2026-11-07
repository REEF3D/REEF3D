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
Architect: Hans Bihs
--------------------------------------------------------------------*/

#include"fnpf_theta.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"solver.h"
#include"reini.h"
#include"fnpf_laplace_cds2.h"
#include"fnpf_fsfbc.h"
#include"fnpf_fsfbc_wd.h"
#include"fnpf_prof.h"
#include"reefmg.h"
#include"fnpf_body.h"
#include"fnpf_amr.h"
#include"fnpf_sz_weights.h"
#include"sflow_bicgstab.h"
#include"matrix2D.h"
#include"vec2D.h"
#include<Eigen/Sparse>
#include<Eigen/SparseLU>
#include<memory>

namespace
{
    // the sparse direct solver of the prototype: pattern analysed once, factorised every step
    struct theta_lu
    {
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> lu;
        bool analysed = false;
        int n = 0;
        
        // the rows of the last assembly: row of node (i,k) (-1 inactive), number of Laplace rows,
        // sz stencil, right-hand side and solution
        std::vector<int> rid;
        int nlap = 0;
        int np = 0;
        const double *wsz = nullptr;
        Eigen::SparseMatrix<double> A;
        Eigen::VectorXd b, x;
    };
    std::unique_ptr<theta_lu> tlu;
}

fnpf_theta::fnpf_theta(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_ini(p,c,pgc),fnpf_sigma(p,c,pgc),
                                                      ek(p),fk(p),Ne(p),Nf(p),Neo(p),Nfo(p),wn(p),erk1(p),frk1(p),
                                                      nEm(p),nFm(p),eko(p),fko(p),erk2(p),R0(p),Rs(p),relR(p),relT(p),relRe(p),relTe(p),fnew(p),wnew(p),en(p),
                                                      colD(p),colA(p),colB(p),fsold(p),nEo(p),nFo(p),pxe(p),pxf(p),kE1(p),kF1(p)
{
    gcval=250;
    if(p->j_dir==0)
    gcval=150;
   
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    gcval_eta = 55;
    gcval_fifsf = 60;
    
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }
    
    theta = p->A305;
    lsolv = (p->A306==2) ? 2 : 1;
    beta = theta*theta;
    betasum = 0.0;
    ctol = p->A307;
    cycsum0 = cycsumc = 0;
    npt = 0;
    wsz = nullptr;
    npass = MAX(p->A302,0);
    passtol = p->A304;
    passsum = stepsum = editsum = 0;
    padvsum = 0;
    csmaxsum = 0.0;
    nstep = 0;
    dt_prev = 0.0;
    solvetime = 0.0;
    
    // limits
    int bad=0;
    if(lsolv==2 && p->j_dir==1)     {bad=1; if(p->mpirank==0) cout<<"A 310 5: SparseLU (A 306 2) 2D only"<<endl;}
    if(lsolv==2 && p->M10>1)        {bad=1; if(p->mpirank==0) cout<<"A 310 5: SparseLU (A 306 2) serial only (M 10 1)"<<endl;}
    if(lsolv==1 && p->N10!=1)       {bad=1; if(p->mpirank==0) cout<<"A 310 5: needs REEFMG (N 10 1) or SparseLU (A 306 2)"<<endl;}
    if(p->A324>0)                   {bad=1; if(p->mpirank==0) cout<<"A 310 5: no Picard passes of the Laplace solve (A 324 0; the corrector passes iterate)"<<endl;}
    if(p->A343==3)                  {bad=1; if(p->mpirank==0) cout<<"A 310 5: wetting-drying A 343 0, 1 or 2 (A 343 3 not tested: unstable in Ting & Kirby)"<<endl;}
    if(p->X10>0 && lsolv==2)        {bad=1; if(p->mpirank==0) cout<<"A 310 5: resolved bodies (X 10) with REEFMG only (A 306 1)"<<endl;}
    if(p->X10>0 && p->A302<1)       {bad=1; if(p->mpirank==0) cout<<"A 310 5: resolved bodies (X 10) need corrector passes (A 302 >= 1)"<<endl;}
    // mesh refinement (fnpf_amr_theta.cpp): the patches without coastline or bodies
    if(p->G1>0 && lsolv==2)         {bad=1; if(p->mpirank==0) cout<<"A 310 5: mesh refinement (G 1) with REEFMG only (A 306 1)"<<endl;}
    if(p->G1>0 && p->X10>0)         {bad=1; if(p->mpirank==0) cout<<"A 310 5: mesh refinement (G 1) without resolved bodies (X 10 0)"<<endl;}
    if(p->G1>0 && p->A350>4)        {bad=1; if(p->mpirank==0) cout<<"A 310 5: mesh refinement (G 1) with breaking A 350 0-4 (as RK3)"<<endl;}
    if(p->G1>0 && p->G30==1)        {bad=1; if(p->mpirank==0) cout<<"A 310 5: mesh refinement (G 1) with the patches off the coastline (G 30 0)"<<endl;}
    if(p->G1>0 && p->G7==1)         {bad=1; if(p->mpirank==0) cout<<"A 310 5: mesh refinement (G 1) without subcycling (G 7 0)"<<endl;}
    if(p->G1>0 && p->A328==1)       {bad=1; if(p->mpirank==0) cout<<"A 310 5: mesh refinement (G 1) with explicit sigma cross terms (A 328 0)"<<endl;}
    if(p->A311==3 || p->A327!=1)    {bad=1; if(p->mpirank==0) cout<<"A 310 5: needs the Fornberg sz stencil (A 327 1) and A 311 1, 2, 4, 5 or 6"<<endl;}
    if(p->A328==1 && lsolv==2)      {bad=1; if(p->mpirank==0) cout<<"A 310 5: implicit sigma cross terms (A 328 1) with REEFMG only (A 306 1)"<<endl;}
    if(p->B98>2 || p->B99>2)        {bad=1; if(p->mpirank==0) cout<<"A 310 5: relaxation zones only (B 98 0-2, B 99 0-2)"<<endl;}
    if(theta<0.5)                   {bad=1; if(p->mpirank==0) cout<<"A 310 5: theta (A 305) >= 0.5"<<endl;}
    if(p->A308>0.0 && theta!=0.5)   {bad=1; if(p->mpirank==0) cout<<"A 310 5: the 4th-order phase (A 308) needs theta 0.5 (A 305)"<<endl;}
    if(p->A308>0.0 && p->N51<2)     {bad=1; if(p->mpirank==0) cout<<"A 310 5: the 4th-order phase (A 308) needs the time step N 51 2 or 3"<<endl;}
    
    if(bad==1)
    pgc->final(true);
    
    plap0 = new fnpf_laplace_cds2(p);
    plap = plap0;
    
    // A 343 1: the free-surface conditions with wetting-drying, as RK3
    if(p->A343==0)
    pf = new fnpf_fsfbc(p,c,pgc);
    
    if(p->A343>=1)
    pf = new fnpf_fsfbc_wd(p,c,pgc);
    
    pbody = fnpf_body::create(p,c,pgc);
    // the body keeps the plain solver for its psi solves; the decorated one (geometry before, body band
    // after the solve) for the Dirichlet solves of the step
    plap = pbody->laplace(plap0);
    
    c->vel3D_fill = [this,c](lexer *p, ghostcell *pgc)
    {
        const int state = c->vel_state;
        {fnpf_prof::scope ps_(fnpf_prof::velcalc); velcalc_sig(p,c,pgc,c->Fi);}
        c->vel_state = state;
    };
    
    if(p->mpirank==0)
    cout<<"FNPF time integration A 310 5: theta-Robin, theta = "<<theta<<", corrector passes A 302 "<<npass<<", "<<(lsolv==2 ? "SparseLU" : "REEFMG")<<endl;
}

fnpf_theta::~fnpf_theta()
{
}

void fnpf_theta::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini)
{
    const double g  = fabs(p->W22);
    const double dt = p->dt;
    
    // damping of eta and Fifsf with the viscosity vb of breaking and of the coastline (implicit over the step, as
    // in each RK3 stage): at the start of the step, i.e. as a source of the step, not after the solve; the
    // solve then gives Fi of the damped surface, without a further Laplace solve (until patch 0022 a consistency
    // solve followed the damping after the solve, in the surf zone in almost every step). w^n is that of the
    // surface before the damping, as the tendencies of the step.
    if(p->A350>0)
    {
    double vbm=0.0;
    SLICELOOP4
    vbm = MAX(vbm,fabs(c->vb(i,j)));
    vbm = pgc->globalmax(vbm);
    
    if(vbm>0.0)
    {
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,c->eta,gcval_eta,1.0);}
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,c->Fifsf,gcval_fifsf,1.0);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->eta,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);}
    }
    
    // mesh refinement: the patches level by level with the damped values of the coarser grid around them, then the
    // covered columns of level 0 from the patches
    bool damped_amr = false;
    if(pamr!=nullptr)
    {fnpf_prof::scope ps_(fnpf_prof::damping); damped_amr = pamr->theta_damping(p,pgc);}
    
    if(vbm>0.0 || damped_amr)
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);}
    }
    
    SLICELOOP4
    en(i,j) = c->eta(i,j);
    
    // wave generation targets at t^n+1, the time of the state they relax (the driver computed them at t^n,
    // as for the RK3 stages; a target at t^n shifts the generated wave by dt, lesson 2 of the hand-over)
    {
    const double st = p->simtime;
    p->simtime = st + dt;
    pflow->wavegen_precalc_fnpf(p,c,pgc);
    p->simtime = st;
    }
    
    // tendencies of the full free-surface conditions at t^n (Fx, Fz, Ex of the end of the last step)
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->kfsfbc(p,c,pgc);}
    SLICELOOP4
    ek(i,j) = c->K(i,j);
    
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->dfsfbc(p,c,pgc,c->eta);}
    SLICELOOP4
    fk(i,j) = c->K(i,j);
    
    // resolved bodies (X 10): the loads of the state at t^n and the predictor of the body (Heun, stage 0) to
    // t^n+1; the corrector stage with the loads of the first corrector state (corrector_state)
    psolv_step = psolv;
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->stage(p,c,pgc,psolv,pf,ek,fk,0);}
    
    // nonlinear parts: everything but eta_t = w and Fifsf_t = -g eta
    SLICELOOP4
    {
    wn(i,j) = c->Fz(i,j);
    Ne(i,j) = ek(i,j) - c->Fz(i,j);
    Nf(i,j) = fk(i,j) + g*c->eta(i,j);
    }

    
    // A 308 > 0: Newmark beta of the Fifsf update below 1/4 for a phase error of 4th order (beta 1/12) or
    // smaller than that of Crank-Nicolson, within the stability limit omega dt < 1/sqrt(1/4 - beta) of the
    // grid cut-off mode, omega_max dt = sqrt(3) dt/dt_crit (N 51 2: dt_crit = sqrt(3)/omega_max, with the
    // Doppler shift), with the safety factor c = A 308
    beta = theta*theta;
    if(p->A308>0.0 && p->dt_crit>0.0)
    {
    const double om = sqrt(3.0)*dt/p->dt_crit;
    beta = MIN(0.25, MAX(1.0/12.0, 0.25 - 1.0/(p->A308*p->A308*om*om)));
    }
    betasum += beta;
    
    // Adams-Bashforth 2 with variable step (Euler in the first step)
    double ab0=1.0, ab1=0.0;
    if(nstep>0 && dt_prev>0.0)
    {
    const double r = dt/dt_prev;
    ab0 = 1.0 + 0.5*r;
    ab1 = -0.5*r;
    }
    
    SLICELOOP4
    {
    // pass 0: N at t^n+1/2 extrapolated
    nEm(i,j) = ab0*Ne(i,j) + ab1*Neo(i,j);
    nFm(i,j) = ab0*Nf(i,j) + ab1*Nfo(i,j);
    
    R0(i,j) = c->Fifsf(i,j) - g*dt*c->eta(i,j) - g*(theta-beta)*dt*dt*wn(i,j);
    Rs(i,j) = R0(i,j) - g*theta*dt*dt*nEm(i,j) + dt*nFm(i,j);
    
    // predictor of eta^n+1 and Fifsf^n+1, Adams-Bashforth 2 of the full tendencies: the sigma grid of
    // the solve (an Euler predictor makes an O(dt^2) grid error per step) and the cross terms of pass 0
    erk1(i,j) = c->eta(i,j) + dt*(ab0*ek(i,j) + ab1*eko(i,j));
    frk1(i,j) = c->Fifsf(i,j) + dt*(ab0*fk(i,j) + ab1*fko(i,j));
    
    eko(i,j) = ek(i,j);
    fko(i,j) = fk(i,j);
    }
    
    // mesh refinement: the same on the patches
    if(pamr!=nullptr)
    pamr->theta_predict(p,pgc,theta,beta);
    
    {fnpf_prof::scope ps_(fnpf_prof::relax); relax_affine(p,c,pgc,pflow,0,relRe,relTe);}
    
    SLICELOOP4
    erk1(i,j) = relRe(i,j)*erk1(i,j) + relTe(i,j);
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,erk1,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->fifsf_relax(p,pgc,frk1);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,frk1,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,Rs,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->surface(p,c,pgc,erk1,frk1,gcval_eta,gcval_fifsf);}
    
    // mesh refinement: the predicted state of all grids consistent (covered columns, the cells around the patches)
    if(pamr!=nullptr)
    pamr->theta_surface(p,pgc,erk1,frk1,0);
    
    {fnpf_prof::scope ps_(fnpf_prof::inflow); pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,frk1,erk1);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,erk1,frk1);}
    {fnpf_prof::scope ps_(fnpf_prof::sigma); sigma_update(p,c,pgc,pf,erk1);}
    
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,frk1,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    
    // relaxation of Fifsf in the zones as f -> relR f + relT (affine in f: from f = 0 and f = 1), part of the
    // Robin rows; Fifsf is not relaxed again after the solve
    {fnpf_prof::scope ps_(fnpf_prof::relax); relax_affine(p,c,pgc,pflow,1,relR,relT);}
    
    // pass 0: one Laplace problem with the Robin condition: Fi, Fifsf^n+1 (fnew) and w^n+1 (wnew)
    if(lsolv==2)
    robin_solve(p,c,pgc,psolv,frk1,Rs,fnew);
    else
    robin_solve_mg(p,c,pgc,psolv,frk1,Rs,fnew,0);
    
    // corrector passes: N at t^n+1/2 trapezoidal, the sigma grid of the solve at eta^n+1 of the last
    // solution and the cross terms with the last solution; defect correction with the factorisation of pass 0
    for(int q=0; q<npass; ++q)
    {
    passq = q;
    corrector_state(p,c,pgc,pflow);
    const double chg = (lsolv==2) ? robin_correct(p,c,pgc,psolv,Rs,fnew) : robin_solve_mg(p,c,pgc,psolv,fnew,Rs,fnew,1);
    ++passsum;
    
    // A 304: stop when the change of the solution is below A 304 times its size
    if(passtol>0.0 && chg<=passtol)
    break;
    }
    ++stepsum;
    
    SLICELOOP4
    {
    c->eta(i,j) = en(i,j) + dt*(theta*wnew(i,j) + (1.0-theta)*wn(i,j) + nEm(i,j));
    c->Fifsf(i,j) = fnew(i,j);
    
    Neo(i,j) = Ne(i,j);
    Nfo(i,j) = Nf(i,j);
    }
    
    dt_prev = dt;
    ++nstep;
    
    // wave generation and absorption of eta (Fifsf: in the Robin rows)
    SLICELOOP4
    c->eta(i,j) = relRe(i,j)*c->eta(i,j) + relTe(i,j);
    
    // bodies: eta and Fifsf over the footprint (no physics there, the free-surface node lies in the body)
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->eta,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->surface(p,c,pgc,c->eta,c->Fifsf,gcval_eta,gcval_fifsf);}
    
    // mesh refinement: eta^n+1, Fifsf^n+1 of the patches, all grids consistent
    if(pamr!=nullptr)
    {
    pamr->theta_final(p,pgc,theta);
    pamr->theta_surface(p,pgc,c->eta,c->Fifsf,2);
    }
    
    // the surface of the solve: the edits below (breaking filters, wetting-drying, coastline) change it after the
    // solve, and Fi below would no longer be that of the edited surface (as with the relaxation zones before
    // patch 0006); then a Dirichlet solve on the final grid follows (as the last RK3 stage)
    SLICELOOP4
    {
    erk1(i,j) = c->eta(i,j);
    frk1(i,j) = c->Fifsf(i,j);
    }
    
    {fnpf_prof::scope ps_(fnpf_prof::wetdry); pf->wetdry(p,c,pgc,c->eta,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->eta,gcval_eta);}
    // A 343 >= 1: the coastline relaxation is part of the affine maps of the zones (relax_affine), i.e. of the
    // Robin rows (Fifsf) and of the eta update; here it would act twice
    if(p->A343==0)
    {
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_eta(p,c,pgc,c->eta);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_fi(p,c,pgc,c->Fifsf);}
    }
    // (Fifsf: relaxation in the Robin rows)
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);}
    
    {fnpf_prof::scope ps_(fnpf_prof::breaking); pf->breaking_zones(p,pflow);}
    {fnpf_prof::scope ps_(fnpf_prof::breaking); pf->breaking(p,c,pgc,c->eta,en,c->Fifsf,1.0);}
    if((p->A350==1 || p->A350==4) && p->A352>0)
    {
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->eta,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);}
    }
    
    // mesh refinement: the detection and the filters on the patches (an edit there: the consistency solve below)
    bool edit_amr = false;
    if(pamr!=nullptr)
    {fnpf_prof::scope ps_(fnpf_prof::breaking); edit_amr = pamr->theta_breaking(p,pgc);}
    
    {fnpf_prof::scope ps_(fnpf_prof::inflow); pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,c->Fifsf,c->eta);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::sigma); sigma_update(p,c,pgc,pf,c->eta);}
    
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    
    // consistency solve after edits of the surface
    {
    double emax=0.0, escale=0.0;
    SLICELOOP4
    {
    emax = MAX(emax, MAX(fabs(c->eta(i,j)-erk1(i,j)), fabs(c->Fifsf(i,j)-frk1(i,j))));
    escale = MAX(escale, MAX(fabs(erk1(i,j)), fabs(frk1(i,j))));
    }
    emax = pgc->globalmax(emax);
    escale = pgc->globalmax(escale);
    
    if(emax>1.0e-12*MAX(escale,1.0e-20) || edit_amr)
    {
    {fnpf_prof::scope ps_(fnpf_prof::laplace); plap->start(p,c,pgc,psolv,pf,c->Fi,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    ++editsum;
    }
    }
    
    // w^n+1 for the next step
    {fnpf_prof::scope ps_(fnpf_prof::fsfwvel); pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);}
    if(pamr!=nullptr)
    pamr->theta_end(p,pgc);
    
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    {fnpf_prof::scope ps_(fnpf_prof::velcalc); velcalc_sig(p,c,pgc,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->velocity(p,c,pgc);}
    
    if(c->U!=nullptr)
    {
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel(p,c,pgc,c->U);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel(p,c,pgc,c->V);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel(p,c,pgc,c->W);}
    }
    else
    {
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel_sl(p,c,pgc,c->Us);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel_sl(p,c,pgc,c->Vs);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_vel_sl(p,c,pgc,c->Ws);}
    }
    
    c->vel_state = 2;
    
    // mesh refinement: velocities of the patches, regrid
    if(pamr!=nullptr)
    pamr->step_end(p,c,pgc);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"theta-Robin solve: "<<setprecision(3)<<solvetime<<" s total, corrector passes per step "<<double(passsum)/double(MAX(stepsum,1L))<<", consistency solves "<<editsum<<", advection-preconditioned passes "<<padvsum<<" (max surface Courant "<<csmaxsum<<")";
    if(p->A308>0.0)
    cout<<", beta "<<beta<<" (mean "<<betasum/double(MAX(stepsum,1L))<<")";
    if(lsolv==1)
    cout<<", REEFMG cycles per step: pass 0 "<<double(cycsum0)/double(MAX(stepsum,1L))<<", corrector "<<double(cycsumc)/double(MAX(stepsum,1L));
    cout<<endl;
    }
}

// Corrector: eta^n+1 of the current solution (eta update with Fifsf^n+1 = fnew, w^n+1 = wnew), N^n+1 on it,
// N^n+1/2 = (N^n + N^n+1)/2, the right-hand side R of the Robin rows, and the sigma grid of eta^n+1 with the
// current solution as boundary values for the next assembly
void fnpf_theta::corrector_state(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow)
{
    const double g  = fabs(p->W22);
    const double dt = p->dt;
    
    SLICELOOP4
    {
    erk2(i,j) = en(i,j) + dt*(theta*wnew(i,j) + (1.0-theta)*wn(i,j) + nEm(i,j));
    c->Fz(i,j) = wnew(i,j);
    }
    
    SLICELOOP4
    erk2(i,j) = relRe(i,j)*erk2(i,j) + relTe(i,j);
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,erk2,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->surface(p,c,pgc,erk2,fnew,gcval_eta,gcval_fifsf);}
    
    // mesh refinement: the corrector state of the patches (eta^n+1, N^n+1/2, R), all grids consistent
    if(pamr!=nullptr)
    pamr->theta_corrector(p,pgc,erk2,fnew,theta,beta);
    
    // Fx, Ex, Exu of eta^n+1, Fifsf^n+1 (the upwinding of A315 uses the Fz just set)
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,erk2,fnew);}
    
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->kfsfbc(p,c,pgc);}
    // the surface Courant number of the wet columns, |Fx| dt/dx + |Fy| dt/dy, of this pass: above 0.5 the change of
    // N is preconditioned with the advection (below)
    double cs=0.0;
    SLICELOOP4
    if(p->wet[IJ]==1)
    cs = MAX(cs, fabs(c->Fx(i,j))*dt/p->DXN[IP] + (p->j_dir==1 ? fabs(c->Fy(i,j))*dt/p->DYN[JP] : 0.0));
    cs = pgc->globalmax(cs);
    const bool padv = cs>0.5;
    csmaxsum = MAX(csmaxsum,cs);
    if(padv) ++padvsum;
    SLICELOOP4
    {
    nEo(i,j) = nEm(i,j);
    nFo(i,j) = nFm(i,j);
    }
    
    SLICELOOP4
    {
    nEm(i,j) = 0.5*(Ne(i,j) + c->K(i,j) - wnew(i,j));
    kE1(i,j) = c->K(i,j);
    }
    
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->dfsfbc(p,c,pgc,erk2);}
    SLICELOOP4
    {
    nFm(i,j) = 0.5*(Nf(i,j) + c->K(i,j) + g*erk2(i,j));
    kF1(i,j) = c->K(i,j);
    }
    
    // the change of N of this pass through (I + dt/2 (Fx Dx + Fy Dy))^-1, Dx, Dy upwind: Newton for the advection
    // of the free-surface conditions (eta_t + Fx Ex, Fifsf_t + Fx^2/2: both advect with Fx), so that the passes
    // converge also where |Fx| dt/dx > 1. Without it they diverge at |Fx| dt/dx of about 2 (the swash of Ting &
    // Kirby at 2 x the RK3 step); with it Ting & Kirby runs at 3 x (validation 09).
    if(padv)
    {
    SLICELOOP4
    {
    pxe(i,j) = nEm(i,j) - nEo(i,j);
    pxf(i,j) = nFm(i,j) - nFo(i,j);
    }
    precond_adv(p,c,pgc,pxe);
    precond_adv(p,c,pgc,pxf);
    SLICELOOP4
    {
    nEm(i,j) = nEo(i,j) + pxe(i,j);
    nFm(i,j) = nFo(i,j) + pxf(i,j);
    }
    }
    
    SLICELOOP4
    Rs(i,j) = R0(i,j) - g*theta*dt*dt*nEm(i,j) + dt*nFm(i,j);
    
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,Rs,gcval_fifsf);}
    
    // the sigma grid of eta^n+1
    {fnpf_prof::scope ps_(fnpf_prof::sigma); sigma_update(p,c,pgc,pf,erk2);}
    
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,fnew,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    
    // resolved bodies: the corrector stage of the body (Heun) with the loads of the first corrector state at
    // t^n+1 (its tendencies kE1, kF1, Fi and the sigma grid of eta^n+1); the next passes solve with the body
    // at t^n+1 (the geometry before their assembly)
    if(passq==0)
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->stage(p,c,pgc,psolv_step,pf,kE1,kF1,1);}
}

// The rows of fnpf_laplace_cds2 for the interior nodes k = 0..knoz-1 (assembled only) plus one Robin row per column
// for the surface node:  Fs + alpha sum_m w_m Fi(knoz-m) = R,  alpha = g theta^2 dt^2 sigz, w_m the sz weights
// (fnpf_sz_weights, Fi(knoz) = Fs). The coupling mt of the top row to the surface node is recovered from two
// assemblies with the surface values f and f + 1 (rhs = ... - mt f; the explicit sigma cross terms only use
// differences of the surface values and do not change). M does not depend on f, so the second assembly gives
// M, and rhs is that of the first. The sigma cross terms are those of the current Fi (interior) and fs.
//
// Note: the calls to the other classes (SLICELOOP4 in 2D) leave the member index j = 1, so the loops over the
// rows set j = 0 first.
void fnpf_theta::robin_build(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, slice &fs, slice &R, slice &ftmp)
{
    const double g = fabs(p->W22);
    const int N = p->knoz;

    // --- assembly with Fifsf = fs: rhs
    plap0->assemble_only = true;
    {fnpf_prof::scope ps_(fnpf_prof::laplace); plap0->start(p,c,pgc,psolv,pf,c->Fi,fs);}
    rhsA = c->rhsvec.V;

    // --- assembly with Fifsf = fs + 1: M, and the coupling mt0 = rhs(fs) - rhs(fs+1) of the top rows
    SLICELOOP4
    ftmp(i,j) = fs(i,j) + 1.0;

    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,ftmp,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,ftmp,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    {fnpf_prof::scope ps_(fnpf_prof::laplace); plap0->start(p,c,pgc,psolv,pf,c->Fi,ftmp);}
    plap0->assemble_only = false;

    const int nrow = int(rhsA.size()) < int(c->rhsvec.V.size()) ? int(rhsA.size()) : int(c->rhsvec.V.size());
    mt0.assign(rhsA.size(),0.0);

    for(int r=0;r<nrow;++r)
    mt0[r] = rhsA[r] - c->rhsvec.V[r];

    // the surface values of fs back into Fi
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,fs,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    
    if(!tlu)
    tlu = std::make_unique<theta_lu>();
    
    theta_lu &L = *tlu;
    
    // --- row index of the nodes (as fnpf_laplace_cds2: ILOOP JLOOP KLOOP with flag4 > 0)
    L.rid.assign(p->knox*N, -1);
    int n=0;
    ILOOP
    JLOOP
    KLOOP
    if(p->flag4[IJK]>0)
    {
    L.rid[i*N+k] = n;
    ++n;
    }
    
    const int ncol = p->knox;
    const int ntot = n + ncol;
    L.nlap = n;
    
    // sz weights of the free-surface stencil
    int np = 3;
    if(p->A311==2 || p->A311==4 || p->A311==5) np = 5;
    if(p->A311==6) np = 7;
    while(np>3 && N<np-1) np -= 2;
    static const fnpf_sz_weights szw(p);
    const double *wsz = szw.weights(np);
    L.np = np;
    L.wsz = wsz;
    
    if(wsz==nullptr)
    {
    cout<<"A 310 5: no sz stencil with "<<np<<" points for "<<N<<" sigma layers"<<endl;
    pgc->final(true);
    }
    
    const double *Mp = c->M.p.data();
    const double *Mn = c->M.n.data();
    const double *Ms = c->M.s.data();
    const double *Mt = c->M.t.data();
    const double *Mb = c->M.b.data();
    const double *rhs = rhsA.data();
    const std::vector<int> &rid = L.rid;
    
    std::vector<Eigen::Triplet<double>> T;
    T.reserve(size_t(ntot)*7);
    L.b.resize(ntot);
    
    j=0;
    ILOOP
    KLOOP
    {
        const int r = rid[i*N+k];
        if(r<0)
        continue;
        
        T.emplace_back(r,r,Mp[r]);
        
        if(Ms[r]!=0.0 && i>0 && rid[(i-1)*N+k]>=0)  T.emplace_back(r,rid[(i-1)*N+k],Ms[r]);
        if(Mn[r]!=0.0 && i<ncol-1 && rid[(i+1)*N+k]>=0) T.emplace_back(r,rid[(i+1)*N+k],Mn[r]);
        if(Mb[r]!=0.0 && k>0 && rid[i*N+k-1]>=0) T.emplace_back(r,rid[i*N+k-1],Mb[r]);
        if(Mt[r]!=0.0 && k<N-1 && rid[i*N+k+1]>=0) T.emplace_back(r,rid[i*N+k+1],Mt[r]);
        
        double rr = rhs[r];
        
        if(k==N-1)
        {
        // surface node: unknown instead of Dirichlet
        T.emplace_back(r,n+i,mt0[r]);
        rr += mt0[r]*fs(i,j);
        }
        
        L.b(r) = rr;
    }
    
    // Robin rows
    ILOOP
    {
        const int s = n+i;
        const double alpha = relR(i,j)*g*beta*p->dt*p->dt*p->sigz[IJ];
        
        T.emplace_back(s,s,1.0 + alpha*wsz[0]);
        
        for(int m=1;m<np;++m)
        {
        const int rk = rid[i*N+N-m];
        if(rk>=0)
        T.emplace_back(s,rk,alpha*wsz[m]);
        }
        
        L.b(s) = relR(i,j)*R(i,j) + relT(i,j);
    }
    
    L.A.resize(ntot,ntot);
    L.A.setFromTriplets(T.begin(),T.end());
    L.A.makeCompressed();
}

// Pass 0: assembly on the predicted sigma grid, factorisation and solve
void fnpf_theta::robin_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, slice &fs, slice &R, slice &fout)
{
    const double t0 = pgc->timer();
    
    robin_build(p,c,pgc,psolv,fs,R,fout);
    
    theta_lu &L = *tlu;
    const int ntot = int(L.b.size());
    
    if(L.n!=ntot)
    {
    L.n = ntot;
    L.analysed = false;
    }
    
    if(!L.analysed)
    {
    L.lu.analyzePattern(L.A);
    L.analysed = true;
    }
    
    L.lu.factorize(L.A);
    
    if(L.lu.info()!=Eigen::Success)
    {
    cout<<"A 310 5: factorisation of the theta-Robin system failed"<<endl;
    pgc->final(true);
    }
    
    L.x = L.lu.solve(L.b);
    
    robin_unpack(p,c,pgc,fout);
    
    p->poissoniter = 1;
    solvetime += pgc->timer() - t0;
    p->poissontime = pgc->timer() - t0;
}

// Corrector pass: the system on the current sigma grid (eta^n+1 of the last solution) with the cross terms of
// the last solution and R; one defect-correction step x += A0^-1 (b - A x) with the factorisation A0 of pass 0.
// A - A0 is of the size of the change of the sigma grid, so the step converges fast.
double fnpf_theta::robin_correct(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, slice &R, slice &fout)
{
    const double t0 = pgc->timer();
    theta_lu &L = *tlu;
    const Eigen::VectorXd xold = L.x;
    
    // the factorisation stays that of pass 0; robin_build replaces A and b
    robin_build(p,c,pgc,psolv,fout,R,wnew);
    
    if(int(L.b.size())!=int(xold.size()))
    {
    cout<<"A 310 5: the rows of the corrector differ from those of the solve"<<endl;
    pgc->final(true);
    }
    
    const Eigen::VectorXd res = L.b - L.A*xold;
    const Eigen::VectorXd dx = L.lu.solve(res);
    L.x = xold + dx;
    
    const double change = dx.cwiseAbs().maxCoeff()/MAX(L.x.cwiseAbs().maxCoeff(),1.0e-20);
    
    robin_unpack(p,c,pgc,fout);
    
    p->poissoniter += 1;
    solvetime += pgc->timer() - t0;
    p->poissontime += pgc->timer() - t0;
    
    return change;
}

void fnpf_theta::robin_unpack(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &fout)
{
    const int N = p->knoz;
    const theta_lu &L = *tlu;
    const std::vector<int> &rid = L.rid;
    const int n = L.nlap;
    
    j=0;
    ILOOP
    {
        KLOOP
        {
        const int r = rid[i*N+k];
        if(r>=0)
        c->Fi[FIJK] = L.x(r);
        }
        
        fout(i,j) = L.x(n+i);
        
        double sz = L.wsz[0]*L.x(n+i);
        for(int m=1;m<L.np;++m)
        {
        const int rk = rid[i*N+N-m];
        if(rk>=0)
        sz += L.wsz[m]*L.x(rk);
        }
        
        wnew(i,j) = p->sigz[IJ]*sz;
    }
    
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,fout,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,wnew,gcval_fifsf);}
}

// The points and weights of the sz stencil (fnpf_sz_weights, as A 311): w = sigz sum_m w_m Fi(knoz-m), Fi(knoz) = Fs
void fnpf_theta::sz_stencil(lexer *p, ghostcell *pgc)
{
    if(wsz!=nullptr)
    return;
    
    const int N = p->knoz;
    int np = 3;
    if(p->A311==2 || p->A311==4 || p->A311==5) np = 5;
    if(p->A311==6) np = 7;
    while(np>3 && N<np-1) np -= 2;
    static const fnpf_sz_weights szw(p);
    wsz = szw.weights(np);
    npt = np;
    
    if(wsz==nullptr)
    {
    if(p->mpirank==0)
    cout<<"A 310 5: no sz stencil with "<<np<<" points for "<<N<<" sigma layers"<<endl;
    pgc->final(true);
    }
}

// A 306 1 (REEFMG). The Robin row of a column, with the surface node Fs = Fi(knoz) and Fi(knoz-m) its sz stencil,
//   D Fs + A sum_m>=1 w_m Fi(knoz-m) = B,   A = r g theta^2 dt^2 sigz,  D = 1 + A w_0,  B = r R + T,
// is eliminated into the top Laplace row (k = knoz-1), whose coupling to the surface node is mt:
//   Fs = (B - A sum_m>=1 w_m Fi(knoz-m))/D,
//   top row: diagonal - mt A w_1/D, k-1 coupling - mt A w_2/D, rhs + mt fs - mt B/D (the assembly had -mt fs),
//   and for stencils of more than three points the couplings - mt A w_m/D to k = knoz-m, m >= 3, in xcol:
//   REEFMG's BiCGStab multiplies with them, the V-cycle preconditions with the 7-point part.
// D >= 1 (A >= 0, w_0 > 0); r = 0 (Dirichlet zone): Fs = T. In the limit of large steps the top row becomes that
// of a Neumann condition, which the line solve treats exactly. The sigma cross terms of rhs are those of the
// current Fi (interior) and fs, as with SparseLU.
void fnpf_theta::robin_build_mg(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, slice &fs, slice &R)
{
    const double g = fabs(p->W22);
    const int N = p->knoz;
    
    sz_stencil(p,pgc);
    
    // resolved bodies: the geometry on the current sigma grid (hull faces, footprint) before the assembly
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->solve_pre(p,c,pgc);}
    slice *foot = pbody->footprint_mask();
    
    // one assembly with Fifsf = fs, the coupling mt of the top rows to the surface node stored in mts
    if(int(mts.size())!=int(c->M.p.size()))
    mts.assign(c->M.p.size(),0.0);
    
    plap0->assemble_only = true;
    plap0->mtsurf = mts.data();
    {fnpf_prof::scope ps_(fnpf_prof::laplace); plap0->start(p,c,pgc,psolv,pf,c->Fi,fs);}
    plap0->mtsurf = nullptr;
    plap0->assemble_only = false;
    
    double *Mp = c->M.p.data();
    double *Mb = c->M.b.data();
    double *rhs = c->rhsvec.V.data();
    
    c->M.xcol.clear();
    
    int n=0;
    ILOOP
    JLOOP
    {
        colD(i,j) = 0.0;
        
        KLOOP
        if(p->flag4[IJK]>0)
        {
            const int r = n;
            ++n;
            
            if(k!=N-1 || p->wet[IJ]==0 || p->flag7[FIJK]<=0)
            continue;
            
            // body footprint: the surface node is a Dirichlet node with the extended Fifsf (fs)
            if(foot!=nullptr && (*foot)(i,j)>0.5)
            continue;
            
            const double mt = mts[r];
            const double A = relR(i,j)*g*beta*p->dt*p->dt*p->sigz[IJ];
            const double D = 1.0 + A*wsz[0];
            const double B = relR(i,j)*R(i,j) + relT(i,j);
            
            colD(i,j) = D;
            colA(i,j) = A;
            colB(i,j) = B;
            
            rhs[r] += mt*fs(i,j) - mt*B/D;
            Mp[r] += -mt*A*wsz[1]/D;
            
            if(npt>2 && N>1)
            Mb[r] += -mt*A*wsz[2]/D;
            
            for(int m=3;m<npt;++m)
            if(N-m>=0)
            {
            const double v = -mt*A*wsz[m]/D;
            if(v!=0.0)
            c->M.xcol.push_back({i,j,N-1,N-m,v});
            }
        }
    }
}

// A 306 1: build and solve with REEFMG (warm start from the current Fi), Fs and w^n+1 from the Robin rows.
// Returns the change of Fi and Fs relative to their size.
double fnpf_theta::robin_solve_mg(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, slice &fs, slice &R, slice &fout, int pass)
{
    const double t0 = pgc->timer();
    const int N = p->knoz;
    
    robin_build_mg(p,c,pgc,psolv,fs,R);
    
    // Fi and Fs before the solve
    const int nf = p->imax*p->jmax*(p->kmax+2);
    if(int(fiold.size())!=nf)
    fiold.resize(nf);
    
    ILOOP
    JLOOP
    {
    KLOOP
    fiold[FIJK] = c->Fi[FIJK];
    fsold(i,j) = fs(i,j);
    }
    
    // A 307 > 0: inexact solves when corrector passes follow (and in them): the residual is reduced by A 307
    // relative to that of the warm start (REEFMG, tolerance < 0); the passes re-solve the whole problem,
    // so the error of a pass does not accumulate, and the last one is A 307 times its (small) initial residual
    const double tol0 = p->N44;
    if(ctol>0.0 && npass>0)
    p->N44 = -ctol;
    
    // corrector passes: REEFMG keeps the coarse levels of pass 0 (the rows change little), the fine level and the
    // Krylov operator are those of the pass
    if(pass>0)
    reefmg::keep_coarse = true;
    
    // mesh refinement: the rows of the patches and the composite solve of all grids (A 307 as above)
    double dpat=0.0, xpat=0.0;
    if(pamr!=nullptr && pamr->active())
    {fnpf_prof::scope ps_(fnpf_prof::laplace); pamr->theta_solve(p,c,pgc,c->Fi,pass>0 ? 1 : 0,beta,dpat,xpat);}
    else
    {fnpf_prof::scope ps_(fnpf_prof::laplace); psolv->startF(p,pgc,c->Fi,c->rhsvec,c->M,8);}
    reefmg::keep_coarse = false;
    
    // resolved bodies: the body band of Fi
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->solve_post(p,c,pgc,c->Fi);}
    slice *foot = pbody->footprint_mask();
    
    p->N44 = tol0;
    c->M.xcol.clear();
    
    if(pass==0)
    cycsum0 += p->solveriter;
    else
    cycsumc += p->solveriter;
    
    // Fs and w^n+1
    double dmax=0.0, xmax=0.0;
    
    ILOOP
    JLOOP
    {
        if(colD(i,j)>0.0)
        {
        double sm = 0.0;
        for(int m=1;m<npt;++m)
        {
        k = N-m;
        sm += wsz[m]*c->Fi[FIJK];
        }
        
        const double Fs = (colB(i,j) - colA(i,j)*sm)/colD(i,j);
        fout(i,j) = Fs;
        wnew(i,j) = p->sigz[IJ]*(wsz[0]*Fs + sm);
        
        dmax = MAX(dmax,fabs(Fs-fsold(i,j)));
        xmax = MAX(xmax,fabs(Fs));
        }
        else
        if(foot!=nullptr && (*foot)(i,j)>0.5)
        {
        // body footprint: the Dirichlet value of the solve (the extended Fifsf), no vertical velocity
        fout(i,j) = fs(i,j);
        wnew(i,j) = 0.0;
        }
        else
        {
        // dry column (A 343 >= 1, not in the solve), no vertical velocity: A 343 1 Fifsf of the map of the land
        // cells (0); A 343 2/3 Fifsf held (wetdry_dynamic sets the dry cells after the step, as for RK3)
        fout(i,j) = (p->A343>=2) ? c->Fifsf(i,j) : relR(i,j)*R(i,j) + relT(i,j);
        wnew(i,j) = 0.0;
        }
        
        KLOOP
        if(p->flag4[IJK]>0)
        {
        dmax = MAX(dmax,fabs(c->Fi[FIJK]-fiold[FIJK]));
        xmax = MAX(xmax,fabs(c->Fi[FIJK]));
        }
    }
    
    dmax = pgc->globalmax(MAX(dmax,dpat));
    xmax = pgc->globalmax(MAX(xmax,xpat));
    
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,fout,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,wnew,gcval_fifsf);}
    
    if(pass==0)
    p->poissoniter = p->solveriter;
    else
    p->poissoniter += p->solveriter;
    
    solvetime += pgc->timer() - t0;
    p->poissontime = (pass==0) ? pgc->timer() - t0 : p->poissontime + pgc->timer() - t0;
    
    return dmax/MAX(xmax,1.0e-20);
}

// The relaxation zones act as f -> r f + (1-r) f_target. iowave applies them to a field; here the affine map
// is read off by applying them to f = 0 and f = 1 (which = 0: eta, 1: Fifsf), and with A 303 > 0 the factors
// become r^(dt/A 303), so that over a fixed time the zones relax by the same amount whatever the step
// (as B 95 for SFLOW and NHFLOW). Not for zones with a state (B 98 6, B 99 6: refused).
void fnpf_theta::relax_affine(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, int which, slice &R, slice &T)
{
    SLICELOOP4
    {
    T(i,j) = 0.0;
    R(i,j) = 1.0;
    }
    
    // A 343 >= 1: the coastline relaxation (A 341/A 342/A 347, to 0 or to the running mean A 326) and Fifsf = 0 on
    // land, before the zones as in the RK3 stages; affine in f as the zones
    if(p->A343>=1)
    {
    if(which==0)
    {
    pf->coastline_eta(p,c,pgc,T);
    pf->coastline_eta(p,c,pgc,R);
    }
    else
    {
    pf->coastline_fi(p,c,pgc,T);
    pf->coastline_fi(p,c,pgc,R);
    }
    }
    
    if(which==0)
    {
    pflow->eta_relax(p,pgc,T);
    pflow->eta_relax(p,pgc,R);
    }
    else
    {
    pflow->fifsf_relax(p,pgc,T);
    pflow->fifsf_relax(p,pgc,R);
    }
    
    // A 303 > 0: factors r^(dt/A 303)
    // A 303 < 0 (default -0.1): the zones act per time as in an RK3 run at N 47 = -A 303, and per step at least
    // as in an RK3 step: an RK3 step relaxes the deviation from the target in its three stages by
    // F(r) = r (1/3 + 2/3 r (3/4 + r/4)); the theta step, s = dt/(-A 303 dt_crit) of those steps long, by
    // F(r)^max(s,1) (dt_crit: the step of the time-step criterion at N 47 1). Below the step of N 47 = -A 303 the
    // zones are those of an RK3 run with the same step, which gets stronger zones per time with smaller steps;
    // above it they do not weaken with the step.
    double fac = 1.0;
    int rk3eq = 0;
    if(p->A303>0.0)
    fac = p->dt/p->A303;
    if(p->A303<0.0 && p->dt_crit>0.0)
    {
    rk3eq = 1;
    fac = MAX(p->dt/(-p->A303*p->dt_crit), 1.0);
    }
    
    SLICELOOP4
    {
    double r = R(i,j) - T(i,j);
    double t = T(i,j);
    
    if((fac!=1.0 || rk3eq==1) && r>0.0 && r<1.0)
    {
    const double target = t/(1.0-r);
    const double F = (rk3eq==1) ? r*(1.0/3.0 + (2.0/3.0)*r*(0.75 + 0.25*r)) : r;
    r = pow(F,fac);
    t = (1.0-r)*target;
    }
    else
    if(rk3eq==1 && r<=0.0)
    r = 0.0;
    
    R(i,j) = r;
    T(i,j) = t;
    }
}

void fnpf_theta::inidisc_step1(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{
    SLICELOOP4
	c->eta(i,j) = 0.0;
    
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->eta,50);}
    
    etaloc_sig(p,c,pgc);
    sigma_ini(p,c,pgc,pf,c->eta);
}

void fnpf_theta::inidisc_step2(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{	
    for(int qn=0;qn<20;++qn)
    {
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_eta(p,c,pgc,c->eta);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_fi(p,c,pgc,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_fi_ini(p,c,pgc,c->Fifsf);}
    }
    
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->eta,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);}
    
    pf->fsfdisc_ini(p,c,pgc,c->eta,c->Fifsf);
    {fnpf_prof::scope ps_(fnpf_prof::wetdry); pf->wetdry(p,c,pgc,c->eta,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::sigma); sigma_update(p,c,pgc,pf,c->eta);}
    pbody->initialize(p,c,pgc);
    
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
    
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    
    if(p->I30>=1)
    {
    {fnpf_prof::scope ps_(fnpf_prof::fsfwvel); pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::velcalc); velcalc_sig(p,c,pgc,c->Fi);}
    }
}

// mesh refinement (G 1): the patches take the theta step (fnpf_amr_theta.cpp); the Dirichlet solves (consistency
// solve, initial state) on all grids through the decorated Laplace solver, as fnpf_RK3
void fnpf_theta::attach_amr(fnpf_amr *a)
{
    pamr = a;
    a->attach_level0(pf,this,this,this);
    a->attach_body(pbody);
    pbody->amr_attach(a);
    plap = a->laplace(plap,plap0);
}

void fnpf_theta::ini_wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
    {fnpf_prof::scope ps_(fnpf_prof::wetdry); pf->wetdry(p,c,pgc,c->eta,c->Fifsf);}
}

// x <- (I + dt/2 (Fx Dx + Fy Dy))^-1 x, first-order upwind differences with the sign of Fx, Fy; no coupling across
// dry or boundary cells (the correction is zero there)
void fnpf_theta::precond_adv(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &x)
{
    if(!psolv2d)
    {
    psolv2d = std::make_unique<sflow_bicgstab>(p,pgc);
    int nr=0;
    SLICELOOP4
    ++nr;
    Nadv = std::make_unique<matrix2D>(p,nr);
    xadv = std::make_unique<vec2D>(p,nr);
    radv = std::make_unique<vec2D>(p,nr);
    }
    matrix2D &NA = *Nadv;
    
    const double h = 0.5*p->dt;
    int n=0;
    SLICELOOP4
    {
    const double ux = c->Fx(i,j), uy = (p->j_dir==1) ? c->Fy(i,j) : 0.0;
    const double ax = h*fabs(ux)/p->DXN[IP];
    const double ay = (p->j_dir==1) ? h*fabs(uy)/p->DYN[JP] : 0.0;
    const int wet = p->wet[IJ];
    
    NA.p[n] = 1.0 + (wet ? ax + ay : 0.0);
    NA.s[n] = (wet && ux>0.0 && p->flagslice4[Im1J]>0 && p->wet[Im1J]) ? -ax : 0.0;
    NA.n[n] = (wet && ux<0.0 && p->flagslice4[Ip1J]>0 && p->wet[Ip1J]) ? -ax : 0.0;
    NA.e[n] = (wet && uy>0.0 && p->flagslice4[IJm1]>0 && p->wet[IJm1]) ? -ay : 0.0;
    NA.w[n] = (wet && uy<0.0 && p->flagslice4[IJp1]>0 && p->wet[IJp1]) ? -ay : 0.0;
    radv->V[n] = x(i,j);
    ++n;
    }
    
    pgc->gcsl_start4(p,x,gcval_fifsf);
    const double tol0 = p->N44;
    p->N44 = 1.0e-10;
    psolv2d->start(p,pgc,x,NA,*xadv,*radv,4);
    p->N44 = tol0;
    pgc->gcsl_start4(p,x,gcval_fifsf);
}

