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

#include "fnpf_fsfbc_wd.h"
#include "lexer.h"
#include "fdm_fnpf.h"
#include "ghostcell.h"
#include "field4.h"
#include "convection.h"
#include "convection.h"
#include "ioflow.h"
#include "solver.h"
#include "reini.h"
#include "fnpf_coastline.h"
#include "sflow_bicgstab.h"
#include "wind_f.h"
#include "wind_v.h"

fnpf_fsfbc_wd::fnpf_fsfbc_wd(lexer *p, fdm_fnpf *c, ghostcell *pgc) : fnpf_breaking(p,c,pgc),wetcoast(p),
                                                                      ef(p),df(p),
                                                                      pconvec(std::in_place_type<fnpf_voiddisc>, p),
                                                                      pdx(std::in_place_type<fnpf_hires>, p),
                                                                      pddx(std::in_place_type<fnpf_ddx_cds2>)
{
    if(p->A311==1)
        pconvec.emplace<fnpf_cds2_wd>(p,c);
    else if(p->A311==2)
        pconvec.emplace<fnpf_cds4_wd>(p);
    else if(p->A311==3)
        pconvec.emplace<fnpf_weno3>(p);
    else if(p->A311==4  || p->A311==5)
    {
        pconvec.emplace<fnpf_weno5_wd>(p);
        pconeta.emplace(p);
        dqF.emplace(p);
        dqE.emplace(p);
    }
    else if(p->A311==6)
        pconvec.emplace<fnpf_cds6_wd>(p);

    if(p->A312==3)
    {
        pdx.emplace<fnpf_cds4>(p);
        pddx.emplace<fnpf_ddx_cds4>();
    }

    pcoast = new fnpf_coastline(p);

    if(p->A350>0)
        psolv = new sflow_bicgstab(p,pgc);

    // wind forcing
    if(p->A370>0)
        pwind = new wind_f(p);
    else // p->A370==0
        pwind = new wind_v(p);

    coastline_count = 0;

    c->wd_criterion=p->A344;

    dist3=0.0;

    if(p->A341>0.0 && p->j_dir==0)
        dist3=p->A341*(p->DXD);
    else if(p->A341>0.0 && p->j_dir==1)
        dist3=p->A341*0.5*(p->DXD+p->DYD);

    if(p->A342>0.0)
        dist3=p->A342;

    dist4=1.0*dist3;

    dist5=3.0*dist4;

    expinverse = 1.0/(exp(1.0)-1.0);

    count_n=0;

    p->Iarray(temp,p->imax*p->jmax);

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

fnpf_fsfbc_wd::~fnpf_fsfbc_wd()
{
}

void fnpf_fsfbc_wd::fsfdisc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    c->WL(i,j) = MAX(c->wd_criterion, eta(i,j) + p->wd - c->bed(i,j));

    pgc->gcsl_start4(p,c->WL,50);

    // ef/df: removed. Both slices were copied, halo-exchanged 4x and passed to
    // filter(), which has no (i,j) loop and so only touched one stale cell.
    // Neither ef nor df is read afterwards (df is only read in fsfdisc_ini,
    // which runs before the first fsfdisc call), so this was dead work
    // costing 4 slice halo exchanges per RK stage.

    std::visit([&](auto &conv, auto &ddx, auto &dx)
    {
        // WENO5: Fifsf through the wet-dry aware fnpf_weno5_wd, eta through fnpf_weno5,
        // both on face divided differences computed once per direction
        if constexpr(std::is_same_v<std::decay_t<decltype(conv)>, fnpf_weno5_wd>)
        {
            conv.dsdiffx(Fifsf,*dqF);
            pconeta->dsdiffx(eta,*dqE);

            SLICELOOP4
            WETDRY
            {
                const double uvel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);

                c->Fx(i,j) = conv.dswenox_dq(*dqF,uvel);
                c->Ex(i,j) = pconeta->dswenox_dq(*dqE,uvel);

                c->Exx(i,j) = ddx.sxx(p,eta);

                c->Bx(i,j) = dx.sx(p,c->depth,1.0);
            }

            // 3D
            if(p->j_dir)
            {
                conv.dsdiffy(Fifsf,*dqF);
                pconeta->dsdiffy(eta,*dqE);

                SLICELOOP4
                WETDRY
                {
                    const double vvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);

                    c->Fy(i,j) = conv.dswenoy_dq(*dqF,vvel);
                    c->Ey(i,j) = pconeta->dswenoy_dq(*dqE,vvel);

                    c->Eyy(i,j) = ddx.syy(p,eta);

                    c->By(i,j) = dx.sy(p,c->depth,1.0);
                }
            }
        }
        // 3D
        else if(p->j_dir)
        {
            SLICELOOP4
            WETDRY
            {
                ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
                jvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);

                c->Fx(i,j) = conv.sx(p,Fifsf,ivel);
                c->Fy(i,j) = conv.sy(p,Fifsf,jvel);

                c->Ex(i,j) = conv.sx(p,eta,ivel);
                c->Ey(i,j) = conv.sy(p,eta,jvel);

                c->Exx(i,j) = ddx.sxx(p,eta);
                c->Eyy(i,j) = ddx.syy(p,eta);

                c->Bx(i,j) = dx.sx(p,c->depth,1.0);
                c->By(i,j) = dx.sy(p,c->depth,1.0);
            }
        }
        // 2D
        else
        {
            SLICELOOP4
            WETDRY
            {
                ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);

                c->Fx(i,j) = conv.sx(p,Fifsf,ivel);
                c->Ex(i,j) = conv.sx(p,eta,ivel);

                c->Exx(i,j) = ddx.sxx(p,eta);

                c->Bx(i,j) = dx.sx(p,c->depth,1.0);
            }
        }
    }, pconvec, pddx, pdx);
}

void fnpf_fsfbc_wd::fsfdisc_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    {
        c->Ex(i,j) = 0.0;
        c->Ey(i,j) = 0.0;

        c->Fx(i,j) = 0.0;
        c->Fy(i,j) = 0.0;

        c->K(i,j) = 0.0;
    }

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

void fnpf_fsfbc_wd::fsfwvel(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    // fi
    std::visit([&](auto &conv)
    {
        FFILOOP4
        {
            if(p->wet[IJ]==1)
                c->Fz(i,j) = p->sigz[IJ]*conv.sz(p,c->Fi);
            else if(p->wet[IJ]==0)
                c->Fz(i,j) = 0.0;
        }
    }, pconvec);

    if(p->count>0)
        coastline_Fz(p,c,pgc,c->Fz);
}

void fnpf_fsfbc_wd::kfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    if(p->A314==1)
    {
        SLICELOOP4
        {
            if(p->wet[IJ]==1)
                c->K(i,j) =  c->Fz(i,j);
            else if(p->wet[IJ]==0)
                c->K(i,j) = 0.0;
        }
    }
    else if(p->A314==2)
    {
        SLICELOOP4
        {
            if(p->wet[IJ]==1)
                c->K(i,j) = - c->Fx(i,j)*c->Ex(i,j) - c->Fy(i,j)*c->Ey(i,j)
                            + c->Fz(i,j)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0));
            else if(p->wet[IJ]==0)
                c->K(i,j) = 0.0;
        }
    }
}

void fnpf_fsfbc_wd::dfsfbc(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta)
{
    if(p->A314==1)
    {
        SLICELOOP4
        {
            if(p->wet[IJ]==1)
                c->K(i,j) = - fabs(p->W22)*eta(i,j);

            else if(p->wet[IJ]==0)
                c->K(i,j) = 0.0;
        }
    }
    else if(p->A314==2)
    {
        SLICELOOP4
        {
            if(p->wet[IJ]==1)
                c->K(i,j) = - 0.5*c->Fx(i,j)*c->Fx(i,j) - 0.5*c->Fy(i,j)*c->Fy(i,j)
                            + 0.5*pow(c->Fz(i,j),2.0)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0))
                            - fabs(p->W22)*eta(i,j);
            else if(p->wet[IJ]==0)
                c->K(i,j) = 0.0;
        }
    }

    pwind->wind_forcing_fnpf(p,c,pgc,c->K,eta);
}

