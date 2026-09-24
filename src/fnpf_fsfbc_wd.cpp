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
                                                                      ef(p),df(p),wetage(p),wdfront(p),wd_dvol(p),wd_nwet(p),
                                                                      pconvec(std::in_place_type<fnpf_voiddisc>, p),
                                                                      pdx(std::in_place_type<fnpf_hires>),
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
    wd_flagcount = -1;

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

    // A315==0 legacy: Ex,Ey upwinded by sign of Fx,Fy (zero gradient at zero speed) and used everywhere.
    // A315>=1:
    //   - Exu,Eyu (kinematic FSBC only): upwinded by dH/deta_x = Fx - 2*Fz*Ex, the characteristic speed
    //     of eta_t + Fx*Ex - Fz*(1 + Ex^2 + Ey^2) = 0; previous-stage Fz, Ex used for the sign.
    //   - zero speed gives the symmetric gradient instead of 0.
    //   - Ex,Ey (sigma metrics, breaking, wind, dynamic FSBC): upwinded by sign of Fx,Fy as before (A315==1)
    //     or symmetric (A315==2).
    // A316==1: WENO5 Fx,Fy fall back to a first-order gradient from wet neighbours where the
    //          stencil is not fully wet (legacy A316==0: zero).
    // A317==1: bed curvature Bxx,Byy updated together with Bx,By (legacy A317==0: Bxx=Byy=0).
    const int zerosym = (p->A315==0) ? 0 : 1;
    const int fallback = (p->A316==0) ? 0 : 1;

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

                c->Fx(i,j) = conv.dswenox_dq_wd(*dqF,uvel,zerosym,fallback);

                if(p->A315==0)
                {
                    c->Ex(i,j) = pconeta->dswenox_dq(*dqE,uvel);
                    c->Exu(i,j) = c->Ex(i,j);
                }
                else
                {
                    const double evel = uvel - 2.0*c->Fz(i,j)*c->Ex(i,j);

                    c->Exu(i,j) = pconeta->dswenox_dq_upsym(*dqE,evel);
                    c->Ex(i,j) = (p->A315==2) ? pconeta->dswenox_dq_sym(*dqE) : pconeta->dswenox_dq_upsym(*dqE,uvel);
                }

                c->Exx(i,j) = ddx.sxx(p,eta);

                c->Bx(i,j) = dx.sx(p,c->depth,1.0);

                if(p->A317==1)
                    c->Bxx(i,j) = ddx.sxx(p,c->depth);
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

                    c->Fy(i,j) = conv.dswenoy_dq_wd(*dqF,vvel,zerosym,fallback);

                    if(p->A315==0)
                    {
                        c->Ey(i,j) = pconeta->dswenoy_dq(*dqE,vvel);
                        c->Eyu(i,j) = c->Ey(i,j);
                    }
                    else
                    {
                        const double evel = vvel - 2.0*c->Fz(i,j)*c->Ey(i,j);

                        c->Eyu(i,j) = pconeta->dswenoy_dq_upsym(*dqE,evel);
                        c->Ey(i,j) = (p->A315==2) ? pconeta->dswenoy_dq_sym(*dqE) : pconeta->dswenoy_dq_upsym(*dqE,vvel);
                    }

                    c->Eyy(i,j) = ddx.syy(p,eta);

                    c->By(i,j) = dx.sy(p,c->depth,1.0);

                    if(p->A317==1)
                        c->Byy(i,j) = ddx.syy(p,c->depth);
                }
            }
        }
        else
        {
            // symmetric gradient: identical to sx/sy for the central schemes, mean of both biases for WENO3
            auto sxs = [&](slice &f) {return 0.5*(conv.sx(p,f,1.0) + conv.sx(p,f,-1.0));};
            auto sys = [&](slice &f) {return 0.5*(conv.sy(p,f,1.0) + conv.sy(p,f,-1.0));};
            auto sxu = [&](slice &f, double vel) {return vel!=0.0 ? conv.sx(p,f,vel) : sxs(f);};
            auto syu = [&](slice &f, double vel) {return vel!=0.0 ? conv.sy(p,f,vel) : sys(f);};

            SLICELOOP4
            WETDRY
            {
                ivel = (Fifsf(i+1,j) - Fifsf(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);

                if(p->A315==0)
                {
                    c->Fx(i,j) = conv.sx(p,Fifsf,ivel);
                    c->Ex(i,j) = conv.sx(p,eta,ivel);
                    c->Exu(i,j) = c->Ex(i,j);
                }
                else
                {
                    const double evel = ivel - 2.0*c->Fz(i,j)*c->Ex(i,j);

                    c->Fx(i,j) = sxu(Fifsf,ivel);
                    c->Exu(i,j) = sxu(eta,evel);
                    c->Ex(i,j) = (p->A315==2) ? sxs(eta) : sxu(eta,ivel);
                }

                c->Exx(i,j) = ddx.sxx(p,eta);

                c->Bx(i,j) = dx.sx(p,c->depth,1.0);

                if(p->A317==1)
                    c->Bxx(i,j) = ddx.sxx(p,c->depth);

                // 3D
                if(p->j_dir)
                {
                    jvel = (Fifsf(i,j+1) - Fifsf(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);

                    if(p->A315==0)
                    {
                        c->Fy(i,j) = conv.sy(p,Fifsf,jvel);
                        c->Ey(i,j) = conv.sy(p,eta,jvel);
                        c->Eyu(i,j) = c->Ey(i,j);
                    }
                    else
                    {
                        const double evel = jvel - 2.0*c->Fz(i,j)*c->Ey(i,j);

                        c->Fy(i,j) = syu(Fifsf,jvel);
                        c->Eyu(i,j) = syu(eta,evel);
                        c->Ey(i,j) = (p->A315==2) ? sys(eta) : syu(eta,jvel);
                    }

                    c->Eyy(i,j) = ddx.syy(p,eta);

                    c->By(i,j) = dx.sy(p,c->depth,1.0);

                    if(p->A317==1)
                        c->Byy(i,j) = ddx.syy(p,c->depth);
                }
            }
        }
    }, pconvec, pddx, pdx);

    // A343 2/3, A336 1: at the wet-dry front (dry cell inside the +-3 stencil)
    // the eta gradients are taken from wet cells only and Exx=Eyy=0. The eta WENO has no
    // wet-dry check and would otherwise read the film level of dry cells
    // (bed + A344, i.e. the land slope for A343 2), which enters Fz*(1+Ex^2+Ey^2),
    // the sigma metrics and the steepness breaking criterion. Same rule as the
    // A316 Fifsf fallback: upwind wet face, else the other wet face, else 0.
    if(p->A343>=2 && p->A336==1)
    {
        auto onesided = [](bool bw, bool fw, double gb, double gf, double vel)
        {
            if(vel>0.0 && bw) return gb;
            if(vel<0.0 && fw) return gf;
            if(bw && fw)      return 0.5*(gb+gf);
            if(bw)            return gb;
            if(fw)            return gf;
            return 0.0;
        };

        SLICELOOP4
        if(p->wet[IJ]==1 && wdfront(i,j)==1)
        {
            // x
            const bool bw = (p->wet[Im1J]==1);
            const bool fw = (p->wet[Ip1J]==1);
            const double gb = bw ? (eta(i,j)-eta(i-1,j))/p->DXP[IM1] : 0.0;
            const double gf = fw ? (eta(i+1,j)-eta(i,j))/p->DXP[IP] : 0.0;

            c->Ex(i,j) = onesided(bw,fw,gb,gf, (p->A315==2) ? 0.0 : c->Fx(i,j));
            c->Exu(i,j) = (p->A315==0) ? c->Ex(i,j) : onesided(bw,fw,gb,gf, c->Fx(i,j) - 2.0*c->Fz(i,j)*c->Ex(i,j));
            // no surface curvature in the sigma metrics at the front: the film
            // depth W jumps between columns there and Exx/W (sigxx) turned 2dx
            // noise of a few-mm film into a Laplace blow-up
            c->Exx(i,j) = 0.0;

            // y
            if(p->j_dir==1)
            {
                const bool bs = (p->wet[IJm1]==1);
                const bool fs = (p->wet[IJp1]==1);
                const double hb = bs ? (eta(i,j)-eta(i,j-1))/p->DYP[JM1] : 0.0;
                const double hf = fs ? (eta(i,j+1)-eta(i,j))/p->DYP[JP] : 0.0;

                c->Ey(i,j) = onesided(bs,fs,hb,hf, (p->A315==2) ? 0.0 : c->Fy(i,j));
                c->Eyu(i,j) = (p->A315==0) ? c->Ey(i,j) : onesided(bs,fs,hb,hf, c->Fy(i,j) - 2.0*c->Fz(i,j)*c->Ey(i,j));
                c->Eyy(i,j) = 0.0;
            }
        }
    }
}

void fnpf_fsfbc_wd::fsfdisc_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf)
{
    SLICELOOP4
    {
        c->Ex(i,j) = 0.0;
        c->Ey(i,j) = 0.0;

        c->Exu(i,j) = 0.0;
        c->Eyu(i,j) = 0.0;

        c->Fx(i,j) = 0.0;
        c->Fy(i,j) = 0.0;

        c->K(i,j) = 0.0;
    }

    // A317==0 legacy: Bx left-biased (sx with speed 1.0, overwritten by fsfdisc), Bxx from df,
    // which is never filled, so Bxx=0.
    // A317==1: Bx from the same symmetric operator as in fsfdisc, Bxx from the depth.
    slice &bdd = (p->A317==0) ? static_cast<slice&>(df) : static_cast<slice&>(c->depth);

    std::visit([&](auto &conv, auto &ddx, auto &dx)
    {
        auto bsx = [&](slice &f) {return p->A317==0 ? conv.sx(p,f,1.0) : dx.sx(p,f,1.0);};
        auto bsy = [&](slice &f) {return p->A317==0 ? conv.sy(p,f,1.0) : dx.sy(p,f,1.0);};

        // 3D
        if(p->j_dir)
        {
            SLICELOOP4
            {
                c->Bx(i,j) = bsx(c->depth);
                c->By(i,j) = bsy(c->depth);

                c->Bxx(i,j) = ddx.sxx(p,bdd);
                c->Byy(i,j) = ddx.syy(p,bdd);
            }
        }
        // 2D
        else
        {
            SLICELOOP4
            {
                c->Bx(i,j) = bsx(c->depth);
                c->Bxx(i,j) = ddx.sxx(p,bdd);
            }
        }
    }, pconvec, pddx, pdx);

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
                // Exu,Eyu: upwinded eta gradient (equal to Ex,Ey for A315==0)
                c->K(i,j) = - c->Fx(i,j)*c->Exu(i,j) - c->Fy(i,j)*c->Eyu(i,j)
                            + c->Fz(i,j)*(1.0 + pow(c->Exu(i,j),2.0) + pow(c->Eyu(i,j),2.0));
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

