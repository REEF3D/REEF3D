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

#include "fnpf_fsfbc.h"
#include "lexer.h"
#include "fdm_fnpf.h"
#include "ghostcell.h"
#include "field4.h"
#include "convection.h"
#include "convection.h"
#include "ioflow.h"
#include "solver.h"
#include "reini.h"

#include "sflow_bicgstab.h"

#include "wind_f.h"
#include "wind_v.h"

using namespace std;

fnpf_fsfbc::fnpf_fsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_breaking(p,c,pgc), ef(p), df(p), pconvec(std::in_place_type<fnpf_voiddisc>, p),
                                                                                                            pddx(std::in_place_type<fnpf_ddx_cds2>, p)
{
    if(p->A311==1)
        pconvec.emplace<fnpf_cds2>(p);
    else if(p->A311==2)
        pconvec.emplace<fnpf_cds4>(p);
    else if(p->A311==3)
        pconvec.emplace<fnpf_weno3>(p);
    else if(p->A311==4 || p->A311==5)
    {
        pconvec.emplace<fnpf_weno5>(p);
        dqF.emplace(p);
        dqE.emplace(p);
    }
    else if(p->A311==6)
        pconvec.emplace<fnpf_cds6>(p);

    if(p->A312==3)
        pddx.emplace<fnpf_ddx_cds4>(p);

    if(p->A350>0)
    psolv = new sflow_bicgstab(p,pgc);

    // wind forcing
    if(p->A370>0)
    pwind = new wind_f(p);
    else //p->A370==0
    pwind = new wind_v(p);

    c->wd_criterion=p->A344;

    // 3D
    if(p->j_dir)
    {
        gcval_eta = 55;
        gcval_fifsf = 60;
    }
    // 2D
    else
    {
        gcval_eta = 155;
        gcval_fifsf = 160;
    }
}

fnpf_fsfbc::~fnpf_fsfbc()
{
}

void fnpf_fsfbc::fsfdisc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    c->WL(i,j) = MAX(0.0, eta(i,j) + p->wd - c->bed(i,j));

    pgc->gcsl_start4(p,c->WL,50);

    // ef/df: removed. Both slices were copied, halo-exchanged 4x and passed to
    // filter(), which has no (i,j) loop and so only touched one stale cell.
    // Neither ef nor df is read afterwards (df is only read in fsfdisc_ini,
    // which runs before the first fsfdisc call), so this was dead work
    // costing 4 slice halo exchanges per RK stage.

    std::visit([&](auto &conv, auto &ddx)
    {
        // WENO5: face divided differences once per field and direction, inlined stencils
        if constexpr(std::is_same_v<std::decay_t<decltype(conv)>, fnpf_weno5>)
        {
            conv.dsdiffx(Fifsf,*dqF);
            conv.dsdiffx(eta,*dqE);

            SLICELOOP4
            {
                const double uvel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);

                c->Fx(i,j) = conv.dswenox_dq(*dqF,uvel);
                c->Ex(i,j) = conv.dswenox_dq(*dqE,uvel);

                c->Exx(i,j) = ddx.sxx(p,eta);
            }

            // 3D
            if(p->j_dir)
            {
                conv.dsdiffy(Fifsf,*dqF);
                conv.dsdiffy(eta,*dqE);

                SLICELOOP4
                {
                    const double vvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);

                    c->Fy(i,j) = conv.dswenoy_dq(*dqF,vvel);
                    c->Ey(i,j) = conv.dswenoy_dq(*dqE,vvel);

                    c->Eyy(i,j) = ddx.syy(p,eta);
                }
            }
        }
        // 3D
        else if(p->j_dir)
        {
            SLICELOOP4
            {
                ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
                jvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);

                c->Fx(i,j) = conv.sx(p,Fifsf,ivel);
                c->Fy(i,j) = conv.sy(p,Fifsf,jvel);

                c->Ex(i,j) = conv.sx(p,eta,ivel);
                c->Ey(i,j) = conv.sy(p,eta,jvel);

                c->Exx(i,j) = ddx.sxx(p,eta);
                c->Eyy(i,j) = ddx.syy(p,eta);
            }
        }
        // 2D
        else
        {
            SLICELOOP4
            {
                ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);

                c->Fx(i,j) = conv.sx(p,Fifsf,ivel);
                c->Ex(i,j) = conv.sx(p,eta,ivel);

                c->Exx(i,j) = ddx.sxx(p,eta);
            }
        }
    }, pconvec, pddx);

    pgc->gcsl_start4(p,c->Ex,1);
    pgc->gcsl_start4(p,c->Ey,1);
}

void fnpf_fsfbc::fsfdisc_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    std::visit([&](auto &conv, auto &ddx)
    {
        // 3D
        if(p->j_dir)
        {
            SLICELOOP4
            {
                c->Bx(i,j) = conv.sx(p,c->depth,1.0);
                c->By(i,j) = conv.sy(p,c->depth,1.0);

                c->Bxx(i,j) = ddx.sxx(p,df);
                c->Byy(i,j) = ddx.syy(p,df);
            }
        }
        // 2D
        else
        {
            SLICELOOP4
            {
                c->Bx(i,j) = conv.sx(p,c->depth,1.0);
                c->Bxx(i,j) = ddx.sxx(p,df);
            }
        }
    }, pconvec, pddx);

    pgc->gcsl_start4(p,c->Bx,1);
    pgc->gcsl_start4(p,c->By,1);

    SLICELOOP4
    p->wet[IJ]=1;

    pgc->gcsl_start4Vint(p,p->wet,50);
}

void fnpf_fsfbc::fsfwvel(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // fi
    std::visit([&](auto &conv)
    {
        FFILOOP4
        {
            c->Fz(i,j) = p->sigz[IJ]*conv.sz(p,c->Fi);

            if(p->wet[IJ]==0)
            c->Fz(i,j) = 0.0;
        }
    }, pconvec);
}

void fnpf_fsfbc::kfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->A314==1)
    {
        SLICELOOP4
        c->K(i,j) =  c->Fz(i,j);
    }
    else if(p->A314==2)
    {
        SLICELOOP4
        c->K(i,j) =  - c->Fx(i,j)*c->Ex(i,j) - c->Fy(i,j)*c->Ey(i,j)
                    + c->Fz(i,j)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0));
    }
}

void fnpf_fsfbc::dfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta)
{
    if(p->A314==1)
    {
        SLICELOOP4
        c->K(i,j) =   - fabs(p->W22)*eta(i,j);
    }
    else if(p->A314==2)
    {
        SLICELOOP4
        c->K(i,j) =  - 0.5*c->Fx(i,j)*c->Fx(i,j) - 0.5*c->Fy(i,j)*c->Fy(i,j)
                    + 0.5*pow(c->Fz(i,j),2.0)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0)) - fabs(p->W22)*eta(i,j);
    }

    // Wind
    pwind->wind_forcing_fnpf(p,c,pgc,c->K,eta);
}

void fnpf_fsfbc::wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    p->wet[IJ]=1;

    pgc->gcsl_start4Vint(p,p->wet,50);
}
