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

#include"momentum_rk.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"pressure.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"vrans.h"
#include"solver.h"
#include"reini.h"
#include"picard.h"
#include"picard_f.h"
#include"picard_lsm.h"
#include"picard_void.h"
#include"fluid_update_fsf.h"
#include"fluid_update_fsf_heat.h"
#include"fluid_update_fsf_heat_Bouss.h"
#include"fluid_update_fsf_comp.h"
#include"fluid_update_fsf_concentration.h"
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_df.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"density_rheo.h"

// loop over the cells of velocity component c (0 u, 1 v, 2 w); cpor is the porosity factor CPOR1/2/3
#define MOMRK_LOOP(c, ...) \
    if((c)==0)      {ULOOP {const double cpor=CPOR1; (void)cpor; __VA_ARGS__}} \
    else if((c)==1) {VLOOP {const double cpor=CPOR2; (void)cpor; __VA_ARGS__}} \
    else            {WLOOP {const double cpor=CPOR3; (void)cpor; __VA_ARGS__}}

momentum_rk::momentum_rk(lexer *p, fdm *a, ghostcell *pgc, convection *pconvection, convection *ppfsfdisc, diffusion *pdiffusion,
                         pressure* ppressure, poisson* ppoisson, turbulence *pturbulence, solver *psolver, solver *ppoissonsolver,
                         ioflow *pioflow, heat *&pheat, concentration *&pconc, reini *ppreini, fsi *ppfsi)
                        :momentum_forcing(p),bcmom(p),urk1(p),urk2(p),fx(p),vrk1(p),vrk2(p),fy(p),wrk1(p),wrk2(p),fz(p),
                         ls(nullptr),frk1(nullptr),frk2(nullptr),ur(nullptr),vr(nullptr),wr(nullptr),pd(nullptr),
                         Dtmp(nullptr),Du(nullptr),Dv(nullptr),Dw(nullptr),
                         pupdate(nullptr),ppicard(nullptr)
{
    for(int q=0; q<3; ++q)
    {
    Mx[q]=rox[q]=nullptr;
    My[q]=roy[q]=nullptr;
    Mz[q]=roz[q]=nullptr;
    }

	gcval_u=10;
	gcval_v=11;
	gcval_w=12;

    gcval_phi=51;

    if(p->F50==2)
    gcval_phi=52;

    if(p->F50==3)
    gcval_phi=53;

    if(p->F50==4)
    gcval_phi=54;

	pconvec=pconvection;
    pfsfdisc=ppfsfdisc;
	pdiff=pdiffusion;
	ppress=ppressure;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
    preini=ppreini;
    pfsi=ppfsi;
    p6dof=nullptr;

    // time integration
    scheme = SSP;
    stages = 3;
    levelset = false;

    if(p->N40==2 || p->N40==22 || p->N40==12)
    stages = 2;

    if(p->N40==4 || p->N40==24 || p->N40==44)
    scheme = LOWSTORAGE;

    if(p->N40==2 || p->N40==22 || p->N40==3 || p->N40==23 || p->N40==4 || p->N40==24 || p->N40==33)
    levelset = true;

    conservative = (p->N40==33);

    // SSP-RK2 (Heun) and SSP-RK3 (Shu-Osher)
    if(stages==2)
    {
    ssp_a[0]=0.0;       ssp_b[0]=1.0;
    ssp_a[1]=0.5;       ssp_b[1]=0.5;
    ssp_a[2]=0.0;       ssp_b[2]=0.0;
    }

    if(stages==3)
    {
    ssp_a[0]=0.0;       ssp_b[0]=1.0;
    ssp_a[1]=0.75;      ssp_b[1]=0.25;
    ssp_a[2]=1.0/3.0;   ssp_b[2]=2.0/3.0;
    }

    // low-storage RK3 (Spalart, Moser, Rogers 1991)
    ls_alpha[0]=4.0/15.0;   ls_alpha[1]=1.0/15.0;    ls_alpha[2]=1.0/6.0;
    ls_gamma[0]=8.0/15.0;   ls_gamma[1]=5.0/12.0;    ls_gamma[2]=3.0/4.0;
    ls_zeta[0] =0.0;        ls_zeta[1] =-17.0/60.0;  ls_zeta[2] =-5.0/12.0;

    // level set inside the stages
    if(levelset)
    {
    ls = new field4(p);
    frk1 = new field4(p);
    frk2 = new field4(p);

        if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90==0)
        pupdate = new fluid_update_fsf(p,a,pgc);

        if(p->F30>0 && p->H10==0 && p->W30==1 && p->F300==0 && p->W90==0)
        pupdate = new fluid_update_fsf_comp(p,a,pgc);

        if(p->F30>0 && p->H10>0 && p->W90==0 && p->F300==0 && p->H3==1)
        pupdate = new fluid_update_fsf_heat(p,a,pgc,pheat);

        if(p->F30>0 && p->H10>0 && p->W90==0 && p->F300==0 && p->H3==2)
        pupdate = new fluid_update_fsf_heat_Bouss(p,a,pgc,pheat);

        if(p->F30>0 && p->C10>0 && p->W90==0 && p->F300==0)
        pupdate = new fluid_update_fsf_concentration(p,a,pgc,pconc);

        if(p->F30>0 && p->H10==0 && p->W30==0 && p->F300==0 && p->W90>0)
        pupdate = new fluid_update_rheology(p);

        // single phase (F 30 0) or multiphase: no density/viscosity update from the level set
        if(p->F300>0 || p->F30==0 || pupdate==nullptr)
        pupdate = new fluid_update_void();

        if(p->F46==2)
        ppicard = new picard_f(p);

        if(p->F46==3)
        ppicard = new picard_lsm(p);

        if(p->F46!=2 && p->F46!=3)
        ppicard = new picard_void(p);
    }

    // second-order implicit diffusion
    imex = (p->D20==2 && p->D23==2 && pdiff->apply_available());
    imex_valid = false;

    if(imex)
    {
    Dtmp = new field1(p);
    Du = new field1(p);
    Dv = new field2(p);
    Dw = new field3(p);
    }

    // conservative form
    gcval_ro=1;
    ro_threshold = p->F91*p->W1 + p->W3;

    if(conservative)
    {
    ur = new field1(p);
    vr = new field2(p);
    wr = new field3(p);

        for(int q=0; q<3; ++q)
        {
        Mx[q] = new field1(p);
        My[q] = new field2(p);
        Mz[q] = new field3(p);
        rox[q] = new field1(p);
        roy[q] = new field2(p);
        roz[q] = new field3(p);
        }

        // face density
        if((p->F80==0) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==0)
        pd = new density_f(p);

        if((p->F80==0) && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0 && p->X10==1)
        pd = new density_df(p);

        if(p->F80==0 && p->H10==0 && p->W30==1  && p->F300==0 && p->W90==0)
        pd = new density_comp(p);

        if(p->F80==0 && p->H10>0 && p->F300==0 && p->W90==0)
        pd = new density_heat(p,pheat);

        if(p->F80==0 && p->C10>0 && p->F300==0 && p->W90==0)
        pd = new density_conc(p,pconc);

        if(p->F80>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90==0)
        pd = new density_vof(p);

        if((p->F30>0 && p->H10==0 && p->W30==0  && p->F300==0 && p->W90>0) || p->F300>=1)
        pd = new density_rheo(p);
    }
}

momentum_rk::~momentum_rk()
{
    delete ls;
    delete frk1;
    delete frk2;
    delete ur;
    delete vr;
    delete wr;
    delete Dtmp;
    delete Du;
    delete Dv;
    delete Dw;

    for(int q=0; q<3; ++q)
    {
    delete Mx[q];
    delete My[q];
    delete Mz[q];
    delete rox[q];
    delete roy[q];
    delete roz[q];
    }
}

void momentum_rk::start(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *pp6dof)
{
    p6dof = pp6dof;

    imex_step_setup(p);

    if(scheme==SSP)
    step_ssp(p,a,pgc,pvrans,p6dof);

    if(scheme==LOWSTORAGE)
    step_lowstorage(p,a,pgc,pvrans,p6dof);
}

// ---------------------------------------------------------------------------------------------
// SSP Runge-Kutta
// ---------------------------------------------------------------------------------------------

void momentum_rk::step_ssp(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{
    amr_step_begin(p,a,pgc);

    for(int s=0; s<stages; ++s)
    {
        amr_ls_transport(p,a,pgc,s);

        if(levelset)
        amr_ls_finish(p,a,pgc,s,true,-1);

        amr_momentum(p,a,pgc,pvrans,p6dof,s,true);

        projection(p,a,pgc,velout(a,0,s),velout(a,1,s),velout(a,2,s),ssp_b[s]);

        amr_stage_end(p,a,pgc,s);
    }
}

// ---------------------------------------------------------------------------------------------
// the parts of an SSP step (step_ssp; the mesh refinement cfd_amr calls them grid by grid, with the
// fills of the patches and the composite projection in between)
// ---------------------------------------------------------------------------------------------

void momentum_rk::amr_step_begin(lexer *p, fdm *a, ghostcell *pgc)
{
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);
	pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);

    if(stages==3)
	pflow->rkinflow(p,a,pgc,urk2,vrk2,wrk2);

    // reference volume of the level set at the start of the step (F 46)
    if(levelset)
    ppicard->volcalc(p,a,pgc,a->phi);
}

// face density of the stage (conservative form) and the level-set transport into the stage output
void momentum_rk::amr_ls_transport(lexer *p, fdm *a, ghostcell *pgc, int s)
{
    // face density of the stage, before the level set is advanced
    if(conservative)
    {
    pgc->start4(p,a->ro,gcval_ro);
    face_density(p,a,pgc,RO(0,s),RO(1,s),RO(2,s));
    }

    if(levelset)
    levelset_transport(p,a,pgc,s);
}

void momentum_rk::amr_ls_finish(lexer *p, fdm *a, ghostcell *pgc, int s, bool picard, int iters)
{
    levelset_finish(p,a,pgc,s,picard,iters);
}

// iterations of the reinitialisation of the stage output (iters > 0)
void momentum_rk::amr_ls_reini(lexer *p, fdm *a, ghostcell *pgc, int s, int iters)
{
    p->reini_iter = iters;
    preini->start(a,p,amr_phi_out(s),pgc,pflow);
}

int momentum_rk::amr_reini_iters(lexer *p, int s) const
{
    return (s==stages-1) ? p->F44 : MAX(p->F44-1,1);
}

// convection, sources, diffusion of the three components and the direct forcing of stage s
void momentum_rk::amr_momentum(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof, int s, bool forcing)
{
    const bool final = (s==stages-1);

    if(conservative)
    convection_conservative(p,a,pgc,s);

    for(int c=0; c<3; ++c)
    component_ssp(p,a,pgc,pvrans,c,s);

    field &uout = velout(a,0,s);
    field &vout = velout(a,1,s);
    field &wout = velout(a,2,s);

    if(conservative)
    {
    pgc->start1(p,uout,gcval_u);
    pgc->start2(p,vout,gcval_v);
    pgc->start3(p,wout,gcval_w);
    }

    if(forcing)
    momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                           uout, vout, wout, fx, fy, fz, s, ssp_b[s], final);
}

// after the projection: the level set of the stage, density and viscosity, the explicit diffusion
// of the second-order implicit scheme
void momentum_rk::amr_stage_end(lexer *p, fdm *a, ghostcell *pgc, int s)
{
    const bool final = (s==stages-1);

    if(conservative)
    clear_FGH(p,a);

    if(levelset)
    {
        field4 &fout = final ? *ls : (s==0 ? *frk1 : *frk2);

        LOOP
        a->phi(i,j,k) = fout(i,j,k);

        pgc->start4(p,a->phi,gcval_phi);
        pupdate->start(p,a,pgc,velout(a,0,s),velout(a,1,s),velout(a,2,s));
    }

    imex_accumulate(p,a,pgc,s,velout(a,0,s),velout(a,1,s),velout(a,2,s));
}

field& momentum_rk::amr_vel(fdm *a, int c, int s)
{
    return vel(a,c,s);
}

field& momentum_rk::amr_velout(fdm *a, int c, int s)
{
    return velout(a,c,s);
}

field4& momentum_rk::amr_phi_in(int s)
{
    return (s==0) ? *ls : (s==1 ? *frk1 : *frk2);
}

field4& momentum_rk::amr_phi_out(int s)
{
    return (s==stages-1) ? *ls : (s==0 ? *frk1 : *frk2);
}

void momentum_rk::levelset_transport(lexer *p, fdm *a, ghostcell *pgc, int s)
{
    const bool final = (s==stages-1);
    field4 &fout = final ? *ls : (s==0 ? *frk1 : *frk2);

    if(s==0)
    {
        LOOP
        {
        a->L(i,j,k)=0.0;
        (*ls)(i,j,k)=a->phi(i,j,k);
        }

        pfsfdisc->start(p,a,*ls,4,a->u,a->v,a->w);

        LOOP
        fout(i,j,k) = (*ls)(i,j,k)
                    + p->dt*a->L(i,j,k);
    }

    if(s>0)
    {
        field4 &fin = (s==1) ? *frk1 : *frk2;

        LOOP
        a->L(i,j,k)=0.0;

        pfsfdisc->start(p,a,fin,4,vel(a,0,s),vel(a,1,s),vel(a,2,s));

        const double as=ssp_a[s];
        const double bs=ssp_b[s];

        LOOP
        fout(i,j,k) = as*(*ls)(i,j,k)
                    + bs*fin(i,j,k)
                    + bs*p->dt*a->L(i,j,k);
    }
}

void momentum_rk::levelset_finish(lexer *p, fdm *a, ghostcell *pgc, int s, bool picard, int iters)
{
    const bool final = (s==stages-1);
    field4 &fout = final ? *ls : (s==0 ? *frk1 : *frk2);

    pflow->phi_relax(p,pgc,fout);

    pgc->start4(p,fout,gcval_phi);

    // F44 iterations in the final stage, one less in the others (default 3 / 2)
    p->reini_iter = final ? p->F44 : MAX(p->F44-1,1);
    if(iters>0)
    p->reini_iter = iters;
    preini->start(a,p,fout,pgc,pflow);
    
    // volume correction once per time step (it also adds Qi*dt-Qo*dt to the target volume)
    if(final && picard)
    ppicard->correct_ls(p,a,pgc,fout);
}

void momentum_rk::component_ssp(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, int c, int s)
{
	starttime=pgc->timer();

    field &un  = comp(a,a->u,a->v,a->w,c);     // u(n)
    field &uin = vel(a,c,s);                    // u(s)
    field &out = velout(a,c,s);                 // u(s+1)
    field &U = vel(a,0,s);
    field &V = vel(a,1,s);
    field &W = vel(a,2,s);
    field &F = FGH(a,c);

    const double as=ssp_a[s];
    const double bs=ssp_b[s];

	sources(p,a,pgc,pvrans,c,uin);
	rhs(p,a,c);

	// combine the explicit terms first, then solve the implicit diffusion with the stage weight
    if(!conservative)
    {
        convection_start(p,a,c,uin,U,V,W);

        if(s==0)
        {
        MOMRK_LOOP(c, out(i,j,k) = un(i,j,k) + p->dt*cpor*F(i,j,k); )
        }

        if(s>0)
        {
        MOMRK_LOOP(c, out(i,j,k) = as*un(i,j,k) + bs*uin(i,j,k) + bs*p->dt*cpor*F(i,j,k); )
        }
    }

    // conservative form: the convection is already in the reconstructed velocity ur
    if(conservative)
    {
        field &urc = UR(c);

        if(s==0)
        {
        MOMRK_LOOP(c, out(i,j,k) = urc(i,j,k) + p->dt*cpor*F(i,j,k); )
        }

        if(s>0)
        {
        MOMRK_LOOP(c, out(i,j,k) = urc(i,j,k) + bs*p->dt*cpor*F(i,j,k); )
        }
    }

    // VRANS resistance, point-implicit with the stage weight (B 268 1)
    pvrans->implicit_drag(p,a,c,bs,out,U,V,W);

    MOMRK_LOOP(c, F(i,j,k) = 0.0; )

    imex_pre(p,a,c,s,out);
    diffusion_start(p,a,pgc,c,out,U,V,W,imex ? imex_g[s] : bs);

	// explicit diffusion (D 20 1) adds to F; zero for the implicit schemes
    if(s==0)
    {
    MOMRK_LOOP(c, out(i,j,k) += p->dt*cpor*F(i,j,k); )
    }

    if(s>0)
    {
    MOMRK_LOOP(c, out(i,j,k) += bs*p->dt*cpor*F(i,j,k); )
    }

    double &t = (c==0) ? p->utime : ((c==1) ? p->vtime : p->wtime);

    if(s==0)
    t = pgc->timer()-starttime;

    if(s>0)
    t += pgc->timer()-starttime;
}

// conservative form: advance u, rho u and the face density rho with the convection of the stage,
// then reconstruct the velocity ur = rho u / rho (limited in the interface region)
void momentum_rk::convection_conservative(lexer *p, fdm *a, ghostcell *pgc, int s)
{
    field &U = vel(a,0,s);
    field &V = vel(a,1,s);
    field &W = vel(a,2,s);

    const double as=ssp_a[s];
    const double bs=ssp_b[s];
    const int so = (s==stages-1) ? 0 : s+1;     // index of the stage output in M, RO

    if(s==0)
    {
    pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
    }

    // u: convected velocity, into the stage output
    for(int c=0; c<3; ++c)
    convection_start(p,a,c,vel(a,c,s),U,V,W);

    for(int c=0; c<3; ++c)
    {
    field &un  = comp(a,a->u,a->v,a->w,c);
    field &uin = vel(a,c,s);
    field &out = velout(a,c,s);
    field &F = FGH(a,c);

        if(s==0)
        {
        MOMRK_LOOP(c, out(i,j,k) = un(i,j,k) + p->dt*cpor*F(i,j,k); )
        }

        if(s>0)
        {
        MOMRK_LOOP(c, out(i,j,k) = as*un(i,j,k) + bs*uin(i,j,k) + bs*p->dt*cpor*F(i,j,k); )
        }
    }

    clear_FGH(p,a);

    // rho u
    for(int c=0; c<3; ++c)
    {
    field &Ms = M(c,s);
    field &ros = RO(c,s);
    field &uin = vel(a,c,s);

    MOMRK_LOOP(c, Ms(i,j,k) = ros(i,j,k)*uin(i,j,k); )
    }

    pgc->start1(p,M(0,s),gcval_u);
	pgc->start2(p,M(1,s),gcval_v);
	pgc->start3(p,M(2,s),gcval_w);

    for(int c=0; c<3; ++c)
    convection_start(p,a,c,M(c,s),U,V,W);

    for(int c=0; c<3; ++c)
    {
    field &Mn = M(c,0);
    field &Ms = M(c,s);
    field &Mo = M(c,so);
    field &F = FGH(a,c);

        if(s==0)
        {
        MOMRK_LOOP(c, Mo(i,j,k) = Mn(i,j,k) + p->dt*cpor*F(i,j,k); )
        }

        if(s>0)
        {
        MOMRK_LOOP(c, Mo(i,j,k) = as*Mn(i,j,k) + bs*Ms(i,j,k) + bs*p->dt*cpor*F(i,j,k); )
        }
    }

    clear_FGH(p,a);

    // rho
    pgc->start1(p,RO(0,s),gcval_u);
    pgc->start2(p,RO(1,s),gcval_v);
    pgc->start3(p,RO(2,s),gcval_w);

    for(int c=0; c<3; ++c)
    convection_start(p,a,c,RO(c,s),U,V,W);

    for(int c=0; c<3; ++c)
    {
    field &ron = RO(c,0);
    field &ros = RO(c,s);
    field &roo = RO(c,so);
    field &F = FGH(a,c);

        if(s==0)
        {
        MOMRK_LOOP(c, roo(i,j,k) = ron(i,j,k) + p->dt*cpor*F(i,j,k); )
        }

        if(s>0)
        {
        MOMRK_LOOP(c, roo(i,j,k) = as*ron(i,j,k) + bs*ros(i,j,k) + bs*p->dt*cpor*F(i,j,k); )
        }
    }

    clear_FGH(p,a);

    // reconstruct u
    for(int c=0; c<3; ++c)
    {
    field &urc = UR(c);
    field &out = velout(a,c,s);
    field &Mo = M(c,so);
    field &roo = RO(c,so);
    field &ros = RO(c,s);

    MOMRK_LOOP(c, urc(i,j,k) = vel_limiter(p,a,out,Mo,roo,ros); )
    }

    pgc->start1(p,*ur,gcval_u);
	pgc->start2(p,*vr,gcval_v);
	pgc->start3(p,*wr,gcval_w);
}

// ---------------------------------------------------------------------------------------------
// low-storage Runge-Kutta
// ---------------------------------------------------------------------------------------------

void momentum_rk::step_lowstorage(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, sixdof *p6dof)
{
    pflow->discharge(p,a,pgc);
    pflow->inflow(p,a,pgc,a->u,a->v,a->w);

    // reference volume of the level set at the start of the step (F 46)
    if(levelset)
    ppicard->volcalc(p,a,pgc,a->phi);

    for(int s=0; s<3; ++s)
    {
        const bool final = (s==2);
        const double al2 = 2.0*ls_alpha[s];

        pflow->rkinflow(p,a,pgc,urk1,vrk1,wrk1);

        if(levelset)
        levelset_lowstorage(p,a,pgc,s);

        for(int c=0; c<3; ++c)
        component_lowstorage(p,a,pgc,pvrans,c,s);

        momentum_forcing_start(a, p, pgc, p6dof, pfsi,
                               urk1, vrk1, wrk1, fx, fy, fz, s, al2, final);

        ULOOP
        {
        a->u(i,j,k) = urk1(i,j,k);

        if(p->count<10)
        a->maxF = MAX(fabs(al2*CPOR1*fx(i,j,k)), a->maxF);

        p->sfmax = MAX(fabs(al2*CPOR1*fx(i,j,k)), p->sfmax);
        }

        VLOOP
        {
        a->v(i,j,k) = vrk1(i,j,k);

        if(p->count<10)
        a->maxG = MAX(fabs(al2*CPOR2*fy(i,j,k)), a->maxG);

        p->sfmax = MAX(fabs(al2*CPOR2*fy(i,j,k)), p->sfmax);
        }

        WLOOP
        {
        a->w(i,j,k) = wrk1(i,j,k);

        if(p->count<10)
        a->maxH = MAX(fabs(al2*CPOR3*fz(i,j,k)), a->maxH);

        p->sfmax = MAX(fabs(al2*CPOR3*fz(i,j,k)), p->sfmax);
        }

        projection(p,a,pgc,a->u,a->v,a->w,al2);

        imex_accumulate(p,a,pgc,s,a->u,a->v,a->w);
    }
}

void momentum_rk::levelset_lowstorage(lexer *p, fdm *a, ghostcell *pgc, int s)
{
    field4 &Cf = *frk1;     // convection of the previous stage

    // the convection adds to L, and the reinitialisation leaves its last right-hand side in L
    LOOP
    a->L(i,j,k)=0.0;

    pfsfdisc->start(p,a,a->phi,4,a->u,a->v,a->w);

    LOOP
    a->phi(i,j,k) += ls_gamma[s]*p->dt*a->L(i,j,k) + ls_zeta[s]*p->dt*Cf(i,j,k);

    LOOP
    Cf(i,j,k)=a->L(i,j,k);

    pflow->phi_relax(p,pgc,a->phi);

    pgc->start4(p,a->phi,gcval_phi);

    // F44 iterations in the final stage, one less in the others (default 3 / 2)
    p->reini_iter = (s==2) ? p->F44 : MAX(p->F44-1,1);
    preini->start(a,p,a->phi,pgc,pflow);
    
    // volume correction once per time step (it also adds Qi*dt-Qo*dt to the target volume)
    if(s==2)
    ppicard->correct_ls(p,a,pgc,a->phi);

    pupdate->start(p,a,pgc,a->u,a->v,a->w);
}

void momentum_rk::component_lowstorage(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, int c, int s)
{
    starttime=pgc->timer();

    field &un = comp(a,a->u,a->v,a->w,c);      // u(s)
    field &uk = comp(a,urk1,vrk1,wrk1,c);       // u(s+1)
    field &C  = comp(a,urk2,vrk2,wrk2,c);       // convection of the previous stage
    field &F  = FGH(a,c);

    const double al2 = 2.0*ls_alpha[s];
    const double gs = ls_gamma[s];
    const double zs = ls_zeta[s];

    sources(p,a,pgc,pvrans,c,un);
    rhs(p,a,c);

    // combine the explicit terms first (sources, pressure, convection), then solve the implicit diffusion
    MOMRK_LOOP(c, uk(i,j,k) = un(i,j,k) + al2*p->dt*cpor*F(i,j,k); )

    MOMRK_LOOP(c, F(i,j,k) = 0.0; )

    convection_start(p,a,c,un,a->u,a->v,a->w);

    MOMRK_LOOP(c, uk(i,j,k) += gs*p->dt*cpor*F(i,j,k) + zs*p->dt*cpor*C(i,j,k); )

    // VRANS resistance, point-implicit with the source weight 2 alpha_s (B 268 1)
    pvrans->implicit_drag(p,a,c,al2,uk,a->u,a->v,a->w);

    MOMRK_LOOP(c, C(i,j,k) = F(i,j,k); )

    MOMRK_LOOP(c, F(i,j,k) = 0.0; )

    imex_pre(p,a,c,s,uk);
    diffusion_start(p,a,pgc,c,uk,a->u,a->v,a->w,imex ? imex_g[s] : al2);

    // explicit diffusion (D 20 1) adds to F; zero for the implicit schemes
    MOMRK_LOOP(c, uk(i,j,k) += al2*p->dt*cpor*F(i,j,k); )

    double &t = (c==0) ? p->utime : ((c==1) ? p->vtime : p->wtime);

    if(s==0)
    t = pgc->timer()-starttime;

    if(s>0)
    t += pgc->timer()-starttime;
}

// ---------------------------------------------------------------------------------------------
// common parts of a stage
// ---------------------------------------------------------------------------------------------

// ust: the velocity component c of the stage, used by the explicit wall shear (was u(n) in every stage)
void momentum_rk::sources(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, int c, field &ust)
{
    if(c==0)
    {
	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,ust,gcval_u);
	ppress->upgrad(p,a,a->eta,a->eta_n);
    }

    if(c==1)
    {
	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,ust,gcval_v);
	ppress->vpgrad(p,a,a->eta,a->eta_n);
    }

    if(c==2)
    {
	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,ust,gcval_w);
	ppress->wpgrad(p,a,a->eta,a->eta_n);
    }
}

// gravity, external forces and the source terms collected in rhsvec, added to F, G, H
void momentum_rk::rhs(lexer *p, fdm *a, int c)
{
	n=0;

    if(c==0)
	ULOOP
	{
    a->maxF=MAX(fabs(a->rhsvec.V[n] + a->gi),a->maxF);
	a->F(i,j,k) += (a->rhsvec.V[n] + a->gi + p->W29_x + a->Fext(i,j,k))*PORVAL1;

	a->rhsvec.V[n]=0.0;
    a->Fext(i,j,k)=0.0;
	++n;
	}

    if(c==1)
	VLOOP
	{
    a->maxG=MAX(fabs(a->rhsvec.V[n] + a->gj),a->maxG);
	a->G(i,j,k) += (a->rhsvec.V[n] + a->gj + p->W29_y + a->Gext(i,j,k))*PORVAL2;

	a->rhsvec.V[n] = 0.0;
    a->Gext(i,j,k) = 0.0;
	++n;
	}

    if(c==2)
	WLOOP
	{
    a->maxH=MAX(fabs(a->rhsvec.V[n] + a->gk),a->maxH);
	a->H(i,j,k) += (a->rhsvec.V[n] + a->gk + p->W29_z + a->Hext(i,j,k))*PORVAL3;

	a->rhsvec.V[n] = 0.0;
    a->Hext(i,j,k) = 0.0;
	++n;
	}
}

void momentum_rk::convection_start(lexer *p, fdm *a, int c, field &f, field &U, field &V, field &W)
{
    pconvec->start(p,a,f,c+1,U,V,W);
}

void momentum_rk::diffusion_start(lexer *p, fdm *a, ghostcell *pgc, int c, field &f, field &U, field &V, field &W, double weight)
{
    if(c==0)
    pdiff->diff_u(p,a,pgc,psolv,f,f,U,V,W,weight);

    if(c==1)
    pdiff->diff_v(p,a,pgc,psolv,f,f,U,V,W,weight);

    if(c==2)
    pdiff->diff_w(p,a,pgc,psolv,f,f,U,V,W,weight);
}

void momentum_rk::projection(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w, double weight)
{
    pflow->pressure_io(p,a,pgc);
	ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,u,v,w,weight);

    amr_project_after(p,a,pgc,u,v,w);
}

// after the pressure correction: relaxation zones and the ghost cells of the velocities
void momentum_rk::amr_project_after(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
	pflow->u_relax(p,a,pgc,u);
	pflow->v_relax(p,a,pgc,v);
	pflow->w_relax(p,a,pgc,w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,u,gcval_u);
	pgc->start2(p,v,gcval_v);
	pgc->start3(p,w,gcval_w);
}

// ---------------------------------------------------------------------------------------------
// second-order implicit diffusion
// ---------------------------------------------------------------------------------------------

void momentum_rk::imex_step_setup(lexer *p)
{
    if(!imex)
    return;

    for(int s=0; s<3; ++s)
    {
    imex_use[s]=false;
    imex_w[s]=0.0;
    }

    if(scheme==SSP && stages==3)
    {
    // (u1, u2, u3): [1], [1/4, 1/4], [-1/2, 1, 1/2]; Shu-Osher inherits b3*(row 2) = [1/6, 1/6]
    imex_g[0]=1.0;      imex_w[0]=-2.0/3.0;
    imex_g[1]=0.25;     imex_w[1]=5.0/6.0;
    imex_g[2]=0.5;      imex_use[2]=true;
    }

    if(scheme==SSP && stages==2)
    {
        // (u0, u1, u2): [0, 1], [1/2, -1/2, 1]; Shu-Osher inherits b2*(row 1) = [0, 1/2]
        if(imex_valid)
        {
        imex_g[0]=1.0;      imex_w[0]=-1.0;
        imex_g[1]=1.0;      imex_use[1]=true;
        }

        if(!imex_valid)
        {
        imex_g[0]=1.0;
        imex_g[1]=0.5;
        }

    imex_w[1]=0.5;      // D(u(n+1)) for the next step
    }

    if(scheme==LOWSTORAGE)
    {
    // stage s: beta_s D(u_s) + (2 alpha_s - beta_s) D(u_s+1)
    const double beta[3] = {4.0/15.0, 0.0, 29.0/150.0};

        for(int s=0; s<3; ++s)
        {
        double b = beta[s];

        if(s==0 && !imex_valid)
        b = 0.0;

        imex_g[s] = 2.0*ls_alpha[s] - b;
        imex_use[s] = (b!=0.0);
        imex_w[s] = beta[(s+1)%3];
        }
    }

    imex_valid = (scheme==LOWSTORAGE || stages==2);
}

// before the implicit solve of stage s: add the accumulated explicit diffusion of earlier stages
void momentum_rk::imex_pre(lexer *p, fdm *a, int c, int s, field &f)
{
    if(!imex || !imex_use[s])
    return;

    field &D = DC(c);

    MOMRK_LOOP(c, f(i,j,k) += p->dt*D(i,j,k); )

    // used: clear all cells (also those that are not fluid in this stage)
    ILOOP
    JLOOP
    KLOOP
    D(i,j,k) = 0.0;
}

// end of stage s: accumulate the weighted diffusion of the projected stage velocities for later stages
void momentum_rk::imex_accumulate(lexer *p, fdm *a, ghostcell *pgc, int s, field &U, field &V, field &W)
{
    if(!imex || imex_w[s]==0.0)
    return;

    const double ws = imex_w[s];

    for(int c=0; c<3; ++c)
    {
        if(c==1 && p->j_dir==0)
        continue;

        if(c==0)
        pdiff->apply_u(p,a,pgc,*Dtmp,U,V,W);

        if(c==1)
        pdiff->apply_v(p,a,pgc,*Dtmp,U,V,W);

        if(c==2)
        pdiff->apply_w(p,a,pgc,*Dtmp,U,V,W);

    field &D = DC(c);
    field &T = *Dtmp;

    MOMRK_LOOP(c, D(i,j,k) += ws*T(i,j,k); )
    }
}

field& momentum_rk::DC(int c)
{
    if(c==0)
    return *Du;

    if(c==1)
    return *Dv;

    return *Dw;
}

// ---------------------------------------------------------------------------------------------
// stage fields
// ---------------------------------------------------------------------------------------------

field& momentum_rk::comp(fdm *a, field &u, field &v, field &w, int c)
{
    if(c==0)
    return u;

    if(c==1)
    return v;

    return w;
}

field& momentum_rk::vel(fdm *a, int c, int s)
{
    if(s==0)
    return comp(a,a->u,a->v,a->w,c);

    if(s==1)
    return comp(a,urk1,vrk1,wrk1,c);

    return comp(a,urk2,vrk2,wrk2,c);
}

field& momentum_rk::velout(fdm *a, int c, int s)
{
    if(s==stages-1)
    return comp(a,a->u,a->v,a->w,c);

    if(s==0)
    return comp(a,urk1,vrk1,wrk1,c);

    return comp(a,urk2,vrk2,wrk2,c);
}

field& momentum_rk::FGH(fdm *a, int c)
{
    if(c==0)
    return a->F;

    if(c==1)
    return a->G;

    return a->H;
}

void momentum_rk::clear_FGH(lexer *p, fdm *a)
{
    ULOOP
    a->F(i,j,k) = 0.0;

    VLOOP
    a->G(i,j,k) = 0.0;

    WLOOP
    a->H(i,j,k) = 0.0;
}

// ---------------------------------------------------------------------------------------------
// conservative form
// ---------------------------------------------------------------------------------------------

field& momentum_rk::M(int c, int idx)
{
    if(c==0)
    return *Mx[idx];

    if(c==1)
    return *My[idx];

    return *Mz[idx];
}

field& momentum_rk::RO(int c, int idx)
{
    if(c==0)
    return *rox[idx];

    if(c==1)
    return *roy[idx];

    return *roz[idx];
}

field& momentum_rk::UR(int c)
{
    if(c==0)
    return *ur;

    if(c==1)
    return *vr;

    return *wr;
}

void momentum_rk::face_density(lexer *p, fdm *a, ghostcell *pgc, field &rx, field &ry, field &rz)
{
    ULOOP
    rx(i,j,k) = pd->roface(p,a,1,0,0);

    VLOOP
    ry(i,j,k) = pd->roface(p,a,0,1,0);

    WLOOP
    rz(i,j,k) = pd->roface(p,a,0,0,1);

    pgc->start1(p,rx,50);
	pgc->start2(p,ry,50);
	pgc->start3(p,rz,50);
}

double momentum_rk::vel_limiter(lexer *p, fdm *a, field &vel, field &Mf, field &ro, field &ro_n)
{
    if(ro(i,j,k)>=ro_threshold)
    val = Mf(i,j,k)/ro(i,j,k);

    else
    if(ro(i,j,k)>p->W3 && ro(i,j,k)<ro_threshold && ro(i,j,k)<ro_n(i,j,k))
    val = (Mf(i,j,k)/ro_filter(p,a,ro))*(ro_filter(p,a,ro)/ro_threshold) + vel(i,j,k)*(ro_threshold-ro_filter(p,a,ro))/ro_threshold;

    else
    val = vel(i,j,k);

    return val;
}

double momentum_rk::ro_filter(lexer *p, fdm *a, field &ro)
{
    if(ro(i,j,k)<p->W3)
    return p->W3;

    return ro(i,j,k);
}
