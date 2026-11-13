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

#include"fnpf_rk.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"solver.h"
#include"fnpf_laplace_cds2.h"
#include"fnpf_fsfbc.h"
#include"fnpf_fsfbc_wd.h"
#include"fnpf_body.h"
#include"fnpf_amr.h"
#include"fnpf_prof.h"
#include"reefmg.h"

fnpf_rk::fnpf_rk(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_ini(p,c,pgc),fnpf_sigma(p,c,pgc),en(p)
{
    gcval=250;
    if(p->j_dir==0)
    gcval=150;

    // 3D
    gcval_eta = 55;
    gcval_fifsf = 60;

    // 2D
    if(p->j_dir==0)
    {
    gcval_eta = 155;
    gcval_fifsf = 160;
    }

    plap0 = new fnpf_laplace_cds2(p);
    plap = plap0;

    if(p->A343==0)
    pf = new fnpf_fsfbc(p,c,pgc);

    if(p->A343>=1)
    pf = new fnpf_fsfbc_wd(p,c,pgc);

    // resolved bodies (X 10): no-op unless present
    pbody = fnpf_body::create(p,c,pgc);
    plap = pbody->laplace(plap);

    // the 3D velocities for an output when they are not stored (fnpf_memory.h): as the last update
    c->vel3D_fill = [this,c](lexer *p, ghostcell *pgc)
    {
        const int state = c->vel_state;

        velcalc_sig(p,c,pgc,c->Fi);

        if(state==2)
        {
        pf->coastline_vel(p,c,pgc,c->U);
        pf->coastline_vel(p,c,pgc,c->V);
        pf->coastline_vel(p,c,pgc,c->W);
        }

        c->vel_state = state;
    };
}

fnpf_rk::~fnpf_rk()
{
}

void fnpf_rk::stage_tendency(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, slice &eta, slice &ek, slice &fk, int s)
{
    // fsf eta
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->kfsfbc(p,c,pgc);}
    SLICELOOP4
    ek(i,j) = c->K(i,j);

    // fsf Fi (dfsfbc only uses the eta passed to it)
    {fnpf_prof::scope ps_(fnpf_prof::fsfbc); pf->dfsfbc(p,c,pgc,eta);}
    SLICELOOP4
    fk(i,j) = c->K(i,j);

    // A 350 5: absorbing surface pressure of the deep-water breakers (the coefficients of t^n)
    if(p->A350==5)
    pf->absorb(p,c,ek,fk);

    // mesh refinement: tendencies of the patches (the body loads need them)
    if(pamr!=nullptr)
    pamr->stage_tendency(p,c,pgc,s);

    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->stage(p,c,pgc,psolv,pf,ek,fk,s);}
}

void fnpf_rk::stage_surface(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, slice &e, slice &f, double alpha, double t, double w,
                            double balpha, int s, bool zones)
{
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,e,gcval_eta,alpha);}
    {fnpf_prof::scope ps_(fnpf_prof::damping); pf->damping(p,c,pgc,f,gcval_fifsf,alpha);}

    // wavegen and coastline (B 95: the targets at the time of the stage output)
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->wavegen_stage_fnpf(p,c,pgc,t,w);}
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->eta_relax(p,pgc,e);}
    {fnpf_prof::scope ps_(fnpf_prof::wetdry); pf->wetdry(p,c,pgc,e,f);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,e,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_eta(p,c,pgc,e);}
    {fnpf_prof::scope ps_(fnpf_prof::coastline); pf->coastline_fi(p,c,pgc,f);}
    {fnpf_prof::scope ps_(fnpf_prof::relax); pflow->fifsf_relax(p,pgc,f);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,f,gcval_fifsf);}
    {fnpf_prof::scope ps_(fnpf_prof::body); pbody->surface(p,c,pgc,e,f,gcval_eta,gcval_fifsf);}

    // mesh refinement: stage values of the patches, restriction, cells around the patches
    if(pamr!=nullptr)
    pamr->stage_surface(p,c,pgc,e,f,s);

    // breaking
    if(zones)
    {fnpf_prof::scope ps_(fnpf_prof::breaking); pf->breaking_zones(p,pflow);}
    {fnpf_prof::scope ps_(fnpf_prof::breaking); pf->breaking(p,c,pgc,e,en,f,balpha);}
    // A 352: the filter changed eta and Fifsf in place after their halo exchange
    if((p->A350==1 || p->A350==4 || p->A350==5) && p->A352>0)
    {
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,e,gcval_eta);}
    {fnpf_prof::scope ps_(fnpf_prof::gc2d); pgc->gcsl_start4(p,f,gcval_fifsf);}
    }

    // fsfdisc and sigma update
    {fnpf_prof::scope ps_(fnpf_prof::inflow); pflow->inflow_fnpf(p,c,pgc,c->Fi,c->Uin,f,e);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfdisc); pf->fsfdisc(p,c,pgc,e,f);}
    {fnpf_prof::scope ps_(fnpf_prof::sigma); sigma_update(p,c,pgc,pf,e);}

    // boundary conditions of Fi
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); fsfbc_sig(p,c,pgc,f,c->Fi);}
    {fnpf_prof::scope ps_(fnpf_prof::bc3d); bedbc_sig(p,c,pgc,c->Fi,pf);}
}

void fnpf_rk::stage_solve(lexer *p, fdm_fnpf *c, ghostcell *pgc, solver *psolv, slice &e, slice &f, int svq, double svc, bool keep)
{
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    if(svq>0) plap0->start_value(p,svq,svc);
    if(keep) reefmg::keep_coarse = true;
    {fnpf_prof::scope ps_(fnpf_prof::laplace); plap->start(p,c,pgc,psolv,pf,c->Fi,f);}
    {fnpf_prof::scope ps_(fnpf_prof::gc3d); pgc->start7V(p,c->Fi,c->bc,gcval);}
    {fnpf_prof::scope ps_(fnpf_prof::fsfwvel); pf->fsfwvel(p,c,pgc,e,f);}
}

void fnpf_rk::step_velocity(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
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
}

void fnpf_rk::inidisc_step1(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
{
    SLICELOOP4
	c->eta(i,j) = 0.0;

    pgc->gcsl_start4(p,c->eta,50);

    etaloc_sig(p,c,pgc);
    sigma_ini(p,c,pgc,pf,c->eta);
}

void fnpf_rk::inidisc_step2(lexer *p, fdm_fnpf *c, ghostcell *pgc, ioflow *pflow, solver *psolv)
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

    // velocity
    if(p->I30>=1)
    {
    pf->fsfwvel(p,c,pgc,c->eta,c->Fifsf);
    velcalc_sig(p,c,pgc,c->Fi);
    pbody->velocity(p,c,pgc);
    }
}
