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
#include"fnpf_body.h"
#include"fnpf_sz_weights.h"
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
    };
    std::unique_ptr<theta_lu> tlu;
}

fnpf_theta::fnpf_theta(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_ini(p,c,pgc),fnpf_sigma(p,c,pgc),
                                                      ek(p),fk(p),Ne(p),Nf(p),Neo(p),Nfo(p),wn(p),erk1(p),frk1(p),
                                                      Rs(p),fnew(p),wnew(p),en(p)
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
    
    theta = p->A300;
    nstep = 0;
    dt_prev = 0.0;
    solvetime = 0.0;
    
    // prototype limits
    int bad=0;
    if(p->j_dir==1)                 {bad=1; if(p->mpirank==0) cout<<"A 310 5: 2D only (prototype)"<<endl;}
    if(p->M10>1)                    {bad=1; if(p->mpirank==0) cout<<"A 310 5: serial only (prototype, M 10 1)"<<endl;}
    if(p->A343!=0)                  {bad=1; if(p->mpirank==0) cout<<"A 310 5: no wetting-drying (A 343 0)"<<endl;}
    if(p->X10>0)                    {bad=1; if(p->mpirank==0) cout<<"A 310 5: no resolved bodies (X 10)"<<endl;}
    if(p->G1>0)                     {bad=1; if(p->mpirank==0) cout<<"A 310 5: no mesh refinement (G 1)"<<endl;}
    if(p->A311==3 || p->A327!=1)    {bad=1; if(p->mpirank==0) cout<<"A 310 5: needs the Fornberg sz stencil (A 327 1) and A 311 1, 2, 4, 5 or 6"<<endl;}
    if(p->A328==1)                  {bad=1; if(p->mpirank==0) cout<<"A 310 5: sigma cross terms explicit only (A 328 0)"<<endl;}
    if(theta<0.5)                   {bad=1; if(p->mpirank==0) cout<<"A 310 5: theta (A 300) >= 0.5"<<endl;}
    
    if(bad==1)
    pgc->final(true);
    
    plap0 = new fnpf_laplace_cds2(p);
    plap = plap0;
    
    pf = new fnpf_fsfbc(p,c,pgc);
    
    pbody = fnpf_body::create(p,c,pgc);
    
    c->vel3D_fill = [this,c](lexer *p, ghostcell *pgc)
    {
        const int state = c->vel_state;
        velcalc_sig(p,c,pgc,c->Fi);
        c->vel_state = state;
    };
    
    if(p->mpirank==0)
    cout<<"FNPF time integration A 310 5: theta-Robin (prototype), theta = "<<theta<<", nonlinear terms Adams-Bashforth 2"<<endl;
}

fnpf_theta::~fnpf_theta()
{
}

void fnpf_theta::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, convection *pconvec, ioflow *pflow, reini *preini)
{
    const double g  = fabs(p->W22);
    const double dt = p->dt;
    
    SLICELOOP4
    en(i,j) = c->eta(i,j);
    
    // tendencies of the full free-surface conditions at t^n (Fx, Fz, Ex of the end of the last step)
    pf->kfsfbc(p,c,pgc);
    SLICELOOP4
    ek(i,j) = c->K(i,j);
    
    pf->dfsfbc(p,c,pgc,c->eta);
    SLICELOOP4
    fk(i,j) = c->K(i,j);
    
    // nonlinear parts: everything but eta_t = w and Fifsf_t = -g eta
    SLICELOOP4
    {
    wn(i,j) = c->Fz(i,j);
    Ne(i,j) = ek(i,j) - c->Fz(i,j);
    Nf(i,j) = fk(i,j) + g*c->eta(i,j);
    }
    
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
    const double nE = ab0*Ne(i,j) + ab1*Neo(i,j);
    const double nF = ab0*Nf(i,j) + ab1*Nfo(i,j);
    
    Rs(i,j) = c->Fifsf(i,j) - g*dt*c->eta(i,j) - g*theta*(1.0-theta)*dt*dt*wn(i,j) - g*theta*dt*dt*nE + dt*nF;
    
    // predictor (sigma grid of the solve, lagged cross terms and boundary values)
    erk1(i,j) = c->eta(i,j) + dt*ek(i,j);
    frk1(i,j) = c->Fifsf(i,j) + dt*fk(i,j);
    
    // keep the effective nonlinear terms for the eta update
    ek(i,j) = nE;
    fk(i,j) = nF;
    }
    
    pflow->eta_relax(p,pgc,erk1);
    pgc->gcsl_start4(p,erk1,gcval_eta);
    pflow->fifsf_relax(p,pgc,frk1);
    pgc->gcsl_start4(p,frk1,gcval_fifsf);
    pgc->gcsl_start4(p,Rs,gcval_fifsf);
    
    pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,frk1,erk1);
    pf->fsfdisc(p,c,pgc,erk1,frk1);
    sigma_update(p,c,pgc,pf,erk1);
    
    fsfbc_sig(p,c,pgc,frk1,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    
    // one Laplace problem with the Robin condition: Fi, Fifsf^n+1 (fnew) and w^n+1 (wnew)
    robin_solve(p,c,pgc,psolv,frk1,Rs,fnew);
    
    SLICELOOP4
    {
    c->eta(i,j) = en(i,j) + dt*(theta*wnew(i,j) + (1.0-theta)*wn(i,j) + ek(i,j));
    c->Fifsf(i,j) = fnew(i,j);
    
    Neo(i,j) = Ne(i,j);
    Nfo(i,j) = Nf(i,j);
    }
    
    dt_prev = dt;
    ++nstep;
    
    // damping, wave generation and absorption, as the last RK3 stage
    pf->damping(p,c,pgc,c->eta,gcval_eta,1.0);
    pf->damping(p,c,pgc,c->Fifsf,gcval_fifsf,1.0);
    
    pflow->eta_relax(p,pgc,c->eta);
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pf->coastline_eta(p,c,pgc,c->eta);
    pf->coastline_fi(p,c,pgc,c->Fifsf);
    pflow->fifsf_relax(p,pgc,c->Fifsf);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    
    pf->breaking_zones(p,pflow);
    pf->breaking(p,c,pgc,c->eta,en,c->Fifsf,1.0);
    if((p->A350==1 || p->A350==4) && p->A352>0)
    {
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    }
    
    pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,c->Fifsf,c->eta);
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    
    // w^n+1 for the next step (the relaxation zones changed Fifsf after the solve)
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
    
    bedbc_sig(p,c,pgc,c->Fi,pf);
    velcalc_sig(p,c,pgc,c->Fi);
    
    if(c->U!=nullptr)
    {
    pf->coastline_vel(p,c,pgc,c->U);
    pf->coastline_vel(p,c,pgc,c->V);
    pf->coastline_vel(p,c,pgc,c->W);
    }
    else
    {
    pf->coastline_vel_sl(p,c,pgc,c->Us);
    pf->coastline_vel_sl(p,c,pgc,c->Vs);
    pf->coastline_vel_sl(p,c,pgc,c->Ws);
    }
    
    c->vel_state = 2;
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"theta-Robin solve: "<<setprecision(3)<<solvetime<<" s total"<<endl;
}

// The rows of fnpf_laplace_cds2 for the interior nodes k = 0..knoz-1 (assembled only) plus one Robin row per column
// for the surface node:  Fs + alpha sum_m w_m Fi(knoz-m) = R,  alpha = g theta^2 dt^2 sigz, w_m the sz weights
// (fnpf_sz_weights, Fi(knoz) = Fs). The coupling mt of the top row to the surface node is recovered from two
// assemblies with the surface values f and f + 1 (rhs = ... - mt f; the explicit sigma cross terms only use
// differences of the surface values and do not change). M does not depend on f, so the second assembly gives
// M, and rhs is that of the first.
void fnpf_theta::robin_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, slice &fs, slice &R, slice &fout)
{
    const double t0 = pgc->timer();
    const double g = fabs(p->W22);
    const int N = p->knoz;

    // --- assembly with Fifsf = fs: rhs of the solve
    plap0->assemble_only = true;
    plap0->start(p,c,pgc,psolv,pf,c->Fi,fs);
    rhsA = c->rhsvec.V;

    // --- assembly with Fifsf = fs + 1: M, and the coupling mt0 = rhs(fs) - rhs(fs+1) of the top rows
    SLICELOOP4
    fout(i,j) = fs(i,j) + 1.0;

    pgc->gcsl_start4(p,fout,gcval_fifsf);
    fsfbc_sig(p,c,pgc,fout,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    plap0->start(p,c,pgc,psolv,pf,c->Fi,fout);
    plap0->assemble_only = false;

    const int nrow = int(rhsA.size()) < int(c->rhsvec.V.size()) ? int(rhsA.size()) : int(c->rhsvec.V.size());
    mt0.assign(rhsA.size(),0.0);

    for(int r=0;r<nrow;++r)
    mt0[r] = rhsA[r] - c->rhsvec.V[r];

    // the surface values of fs back into Fi
    fsfbc_sig(p,c,pgc,fs,c->Fi);
    pgc->start7V(p,c->Fi,c->bc,gcval);
    
    // --- row index of the nodes (as fnpf_laplace_cds2: ILOOP JLOOP KLOOP with flag4 > 0)
    std::vector<int> rid(p->knox*N, -1);
    int n=0;
    ILOOP
    JLOOP
    KLOOP
    if(p->flag4[IJK]>0)
    {
    rid[i*N+k] = n;
    ++n;
    }
    
    const int ncol = p->knox;
    const int ntot = n + ncol;
    
    // sz weights of the free-surface stencil
    int np = 3;
    if(p->A311==2 || p->A311==4 || p->A311==5) np = 5;
    if(p->A311==6) np = 7;
    while(np>3 && N<np-1) np -= 2;
    fnpf_sz_weights szw(p);
    const double *wsz = szw.weights(np);
    
    const double *Mp = c->M.p.data();
    const double *Mn = c->M.n.data();
    const double *Ms = c->M.s.data();
    const double *Mt = c->M.t.data();
    const double *Mb = c->M.b.data();
    const double *rhs = rhsA.data();

    if(wsz==nullptr)
    {
    cout<<"A 310 5: no sz stencil with "<<np<<" points for "<<N<<" sigma layers"<<endl;
    pgc->final(true);
    }
    
    std::vector<Eigen::Triplet<double>> T;
    T.reserve(size_t(ntot)*7);
    Eigen::VectorXd b(ntot), x(ntot);
    
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
        
        b(r) = rr;
    }
    
    // Robin rows
    ILOOP
    {
        const int s = n+i;
        const double alpha = g*theta*theta*p->dt*p->dt*p->sigz[IJ];
        
        T.emplace_back(s,s,1.0 + alpha*wsz[0]);
        
        for(int m=1;m<np;++m)
        {
        const int rk = rid[i*N+N-m];
        if(rk>=0)
        T.emplace_back(s,rk,alpha*wsz[m]);
        }
        
        b(s) = R(i,j);
    }
    
    Eigen::SparseMatrix<double> A(ntot,ntot);
    A.setFromTriplets(T.begin(),T.end());
    A.makeCompressed();
    
    if(!tlu || tlu->n!=ntot)
    {
    tlu = std::make_unique<theta_lu>();
    tlu->n = ntot;
    }
    
    if(!tlu->analysed)
    {
    tlu->lu.analyzePattern(A);
    tlu->analysed = true;
    }
    
    tlu->lu.factorize(A);
    
    if(tlu->lu.info()!=Eigen::Success)
    {
    cout<<"A 310 5: factorisation of the theta-Robin system failed"<<endl;
    pgc->final(true);
    }
    
    x = tlu->lu.solve(b);
    
    // --- back to Fi, Fifsf^n+1 and w^n+1
    ILOOP
    {
        KLOOP
        {
        const int r = rid[i*N+k];
        if(r>=0)
        c->Fi[FIJK] = x(r);
        }
        
        fout(i,j) = x(n+i);
        
        double sz = wsz[0]*x(n+i);
        for(int m=1;m<np;++m)
        {
        const int rk = rid[i*N+N-m];
        if(rk>=0)
        sz += wsz[m]*x(rk);
        }
        
        wnew(i,j) = p->sigz[IJ]*sz;
    }
    
    pgc->gcsl_start4(p,fout,gcval_fifsf);
    pgc->gcsl_start4(p,wnew,gcval_fifsf);
    
    p->poissoniter = 1;
    solvetime += pgc->timer() - t0;
    p->poissontime = pgc->timer() - t0;
}

void fnpf_theta::inidisc_step1(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{
    SLICELOOP4
	c->eta(i,j) = 0.0;
    
    pgc->gcsl_start4(p,c->eta,50);
    
    etaloc_sig(p,c,pgc);
    sigma_ini(p,c,pgc,pf,c->eta);
}

void fnpf_theta::inidisc_step2(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{	
    for(int qn=0;qn<20;++qn)
    {
    pf->coastline_eta(p,c,pgc,c->eta);
    pf->coastline_fi(p,c,pgc,c->Fifsf);
    pf->coastline_fi_ini(p,c,pgc,c->Fifsf);
    }
    
    pgc->gcsl_start4(p,c->eta,gcval_eta);
    pgc->gcsl_start4(p,c->Fifsf,gcval_fifsf);
    
    pf->fsfdisc_ini(p,c,pgc,c->eta,c->Fifsf);
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);   
    pf->fsfdisc(p,c,pgc,c->eta,c->Fifsf);
    sigma_update(p,c,pgc,pf,c->eta);
    pbody->initialize(p,c,pgc);
    
    fsfbc_sig(p,c,pgc,c->Fifsf,c->Fi);
    bedbc_sig(p,c,pgc,c->Fi,pf);
    
    pgc->start7V(p,c->Fi,c->bc,gcval);
    
    if(p->I30>=1)
    {
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
    velcalc_sig(p,c,pgc,c->Fi);
    }
}

void fnpf_theta::ini_wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{	
    pf->wetdry(p,c,pgc,c->eta,c->Fifsf);
}
