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

#include"sediment_exner.h"
#include"sediment_curv.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include<math.h>

/*--------------------------------------------------------------------
Lags of the depth-averaged bend morphology (SFLOW), both steady relaxations along the depth-averaged flow

        f + L s.grad(f) = f_e ,     s = (u,v)/|U|

solved as the non-equilibrium bed load (sediment_exner_noneq.cpp): first-order upwind along s, Gauss-Seidel
with alternating sweep directions, to convergence in every sediment step (S 49 sweeps at most).

S 95 2, spiral-flow intensity with adaptation (Kalkwijk & Booij 1986; the form of Delft3D and HEC-RAS 2D):
     I + L_a s.grad(I) = I_e,   I_e = h |U| K_c,   L_a = S97 (m - 1) h / (2 kappa^2),   m = kappa C/sqrt(g)
     (K_c the signed curvature of the streamlines, C = sqrt(g)|U|/u* from the bed shear, at most 100; L_a = 0 where
     |U| < 1e-4 m/s, I_e = 0 there); the spiral flow of a
     bend builds up over L_a after the entrance and decays over L_a after the exit.  The bed-load direction
     turns by tan(delta) = alpha_I I/|U| (alpha_I of S 95 1), which is S 95 1 where I = I_e.
S 64 f_s > 0, the suspended part of the transport (single fraction, SFLOW): the transport formula gives the
     capacity q_e; the part (1 - f_s) is bed load (equilibrium, or S 33), the part f_s is suspended load and
     relaxes over the adaptation length of the suspension (Galappatti & Vreugdenhil 1985, depth-averaged)
     q_s + L_s s.grad(q_s) = f_s q_e,   L_s = h |U| / (S65 w_s)   (S 65 the ratio c_bed/c_mean, default 1).
     The spiral flow (S 95) and the transverse slope (S 86) turn the total load (bed load plus the lagged suspended
     load; turning the bed load only put the overshoot of the LFM bend of Talmon & de Graaff 1989 far downstream,
     sediment_curvilinear/t7).  The Exner fluxes carry the sum, so the bed volume stays conserved.
--------------------------------------------------------------------*/

int sediment_exner::relax_along(lexer *p, ghostcell *pgc, sediment_fdm *s, slice4 &src, slice4 &out, slice4 &len, int &ini)
{
    const double ydir = p->y_dir;
    const int itermax = MAX(p->S49,1);
    const double tol = lag_tol;
    double uvel,vvel,umag,sgx,sgy,Lsc,cx,cy,num,den,fnew;
    double fmax=0.0, dmax=0.0;
    int iter=0;

    SEDSLICELOOP
    {
        if(p->geo_curv==1)
        sediment_curv::cellvel(p,s,i,j,uvel,vvel);
        else
        {
        uvel = 0.5*(s->P(i,j) + s->P(i-1,j));
        vvel = 0.5*(s->Q(i,j) + s->Q(i,j-1))*ydir;
        }

        umag = sqrt(uvel*uvel + vvel*vvel);
        sgx = umag>1.0e-10?uvel/umag:0.0;
        sgy = umag>1.0e-10?vvel/umag:0.0;

        Lsc = MAX(len(i,j),0.0);

        cx = cy = 0.0;

        // the upwind neighbour: a sediment cell of the grid (or of the next rank, or a cell around the grid); not
        // the ghost cells beyond walls and open boundaries (the cell then takes its own value there)
        auto up = [&](int a, int b) { const int q = (a-p->imin)*p->jmax + (b-p->jmin); return p->DFBED[q]>0 && p->flagslice4[q]>0; };

        if(p->geo_curv==1)
        {
        const double sxi = sediment_curv::contra(p,1,i,j,sgx,sgy);
        const double seta = sediment_curv::contra(p,2,i,j,sgx,sgy);
        if(sxi>0.0 && up(i-1,j))
        cx = Lsc*sxi;
        if(sxi<0.0 && up(i+1,j))
        cx = Lsc*sxi;
        if(seta>0.0 && up(i,j-1))
        cy = Lsc*seta;
        if(seta<0.0 && up(i,j+1))
        cy = Lsc*seta;
        }
        else
        {
        if(sgx>0.0 && up(i-1,j) && p->DXP[IM1]>1.0e-20)
        cx = Lsc*sgx/p->DXP[IM1];
        if(sgx<0.0 && up(i+1,j) && p->DXP[IP]>1.0e-20)
        cx = Lsc*sgx/p->DXP[IP];
        if(sgy>0.0 && up(i,j-1) && p->DYP[JM1]>1.0e-20)
        cy = Lsc*sgy/p->DYP[JM1];
        if(sgy<0.0 && up(i,j+1) && p->DYP[JP]>1.0e-20)
        cy = Lsc*sgy/p->DYP[JP];
        }

        cxn(i,j) = cx;
        cyn(i,j) = cy;

        fmax = MAX(fmax,fabs(src(i,j)));
    }

    fmax = pgc->globalmax(fmax);

    // initial guess: the equilibrium, afterwards the previous sediment step
    if(ini==0)
    {
    SEDSLICELOOP
    out(i,j) = src(i,j);
    ini = 1;
    }

    lag_hold(p,pgc,out,0);

    for(iter=0; iter<itermax; ++iter)
    {
        const int sweep = iter%4;
        const bool xrev = (sweep==1 || sweep==2);
        const bool yrev = (sweep==2 || sweep==3);

        dmax = 0.0;

        for(int ii=0; ii<p->knox; ++ii)
        for(int jj=0; jj<p->knoy; ++jj)
        {
            i = xrev?(p->knox-1-ii):ii;
            j = yrev?(p->knoy-1-jj):jj;

            if(p->flagslice4[IJ]>0 && p->DFBED[IJ]>0 && lagfix(i,j)<=0.0)
            {
            cx = cxn(i,j);
            cy = cyn(i,j);
            num = src(i,j);
            den = 1.0;

            if(cx>0.0)
            {
            num += cx*out(i-1,j);
            den += cx;
            }
            if(cx<0.0)
            {
            num -= cx*out(i+1,j);
            den -= cx;
            }
            if(cy>0.0)
            {
            num += cy*out(i,j-1);
            den += cy;
            }
            if(cy<0.0)
            {
            num -= cy*out(i,j+1);
            den -= cy;
            }

            fnew = num/den;
            dmax = MAX(dmax,fabs(fnew-out(i,j)));
            out(i,j) = fnew;
            }
        }

        lag_hold(p,pgc,out,1);

        dmax = pgc->globalmax(dmax);
        // the change of the first sweep: how far the held cells (the cells around the grid) moved the solution
        if(iter==0)
        {
        lag_d0 = dmax;
        lag_f0 = fmax;
        }

        if(dmax<=tol*fmax)
        {
        ++iter;
        break;
        }
    }

    return iter;
}

// the ghost cells of out; the held cells (lagfix, mesh refinement and multi-block grids) keep their values: the boundary
// lists of a grid can write the cells around it that stand next to a wall
void sediment_exner::lag_hold(lexer *p, ghostcell *pgc, slice4 &out, int again)
{
    if(again==0)
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    if(lagfix(i,j)>0.0)
    lagval(i,j) = out(i,j);

    pgc->gcsl_start4(p,out,1);

    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    if(lagfix(i,j)>0.0)
    out(i,j) = lagval(i,j);
}

// S 64: the suspended part of qb relaxed over L_s, the bed part (1 - f_s) qb; the relaxed suspended load qsn for the
// cells around the grids (amr_field 6)
void sediment_exner::susp_split(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    if(p->S64<=0.0)
    return;

    const double fs = MIN(p->S64,1.0);
    const double ws = MAX(s->ws,1.0e-10);
    const double ydir = p->y_dir;
    double uvel,vvel;

    SLICEBASELOOP
    {
    qse(i,j) = 0.0;
    lss(i,j) = 0.0;
    qbb(i,j) = 0.0;
    }

    SEDSLICELOOP
    {
        if(p->geo_curv==1)
        sediment_curv::cellvel(p,s,i,j,uvel,vvel);
        else
        {
        uvel = 0.5*(s->P(i,j) + s->P(i-1,j));
        vvel = 0.5*(s->Q(i,j) + s->Q(i,j-1))*ydir;
        }

        qse(i,j) = fs*s->qbe(i,j);
        qbb(i,j) = (1.0-fs)*s->qb(i,j);
        lss(i,j) = MAX(s->waterlevel(i,j),0.0)*sqrt(uvel*uvel + vvel*vvel)/(MAX(p->S65,1.0e-10)*ws);
    }

    susp_apply(p,pgc,s);
}

void sediment_exner::susp_apply(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const int it = relax_along(p,pgc,s,qse,qsn,lss,susp_ini);
    lag_d[0] = lag_d0;
    lag_f[0] = lag_f0;

    double lmax = 0.0;
    SEDSLICELOOP
    {
        const double qb = qbb(i,j);
        const double qs = MAX(qsn(i,j),0.0);
        s->qb(i,j) = qb + qs;
        if(p->S62!=2)
        s->qbs(i,j) = qs;     // the output ST_qbs (S 62 2 would add qbs to qb a second time)
        lmax = MAX(lmax,lss(i,j));
    }
    pgc->gcsl_start4(p,s->qb,1);

    lmax = pgc->globalmax(lmax);
    if(p->mpirank==0 && lag_print)
    cout<<"suspended load lag  L_s max: "<<setprecision(4)<<lmax<<"  sweeps: "<<it<<endl;
}

// S 95 2: the spiral-flow intensity I relaxed over L_a towards I_e = h|U|K_c (bend_prep); tan(delta) = alpha_I I/|U|
void sediment_exner::spiral_lag(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const double kappa = 0.4;
    const double g = fabs(p->W22);

    SEDSLICELOOP
    {
        const double u = bu(i,j), v = bv(i,j);
        const double U = sqrt(u*u + v*v);
        const double ustar = s->shearvel_eff(i,j);
        // C between sqrt(g)/kappa (m = 1, no lag) and 100: still water (u* -> 0) would give an unbounded length
        const double C = (ustar>1.0e-10) ? MIN(sqrt(g)*U/ustar,100.0) : 100.0;
        const double m = kappa*C/sqrt(g);
        const double h = MAX(s->waterlevel(i,j),0.0);
        lsa(i,j) = (U<1.0e-4) ? 0.0 : MAX(0.0, p->S97*(m - 1.0)*h/(2.0*kappa*kappa));
    }

    spiral_apply(p,pgc,s);
}

void sediment_exner::spiral_apply(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const double kappa = 0.4;
    const double g = fabs(p->W22);
    const int it = relax_along(p,pgc,s,Ie,Isp,lsa,spiral_ini);
    lag_d[1] = lag_d0;
    lag_f[1] = lag_f0;

    double lmax = 0.0;
    SEDSLICELOOP
    {
        const double u = bu(i,j), v = bv(i,j);
        const double U = sqrt(u*u + v*v);
        btd(i,j) = 0.0;
        lmax = MAX(lmax,lsa(i,j));

        if(U<1.0e-4)
        continue;

        const double ustar = s->shearvel_eff(i,j);
        const double C = (ustar>1.0e-10) ? sqrt(g)*U/ustar : 1.0e10;
        const double alpha = MAX(0.0, p->S96*2.0/(kappa*kappa)*(1.0 - sqrt(g)/(kappa*C)));

        btd(i,j) = MAX(-1.0, MIN(1.0, alpha*Isp(i,j)/U));
    }
    pgc->gcsl_start4(p,btd,1);

    lmax = pgc->globalmax(lmax);
    if(p->mpirank==0 && lag_print)
    cout<<"spiral flow lag  L_a max: "<<setprecision(4)<<lmax<<"  sweeps: "<<it<<endl;
}

// mesh refinement and multi-block grids (sflow_amr::sediment_step_sflow): after the cells around the grids took the
// lagged fields of their sources, the relaxations again from them, so that the lag runs on through the interfaces
// (passes until the cells around the grids no longer move the solution); d[k], f[k]: the change of the first sweep and
// the largest equilibrium of the suspended load (k = 0) and of the spiral flow (k = 1)
void sediment_exner::lag_refresh(lexer *p, ghostcell *pgc, sediment_fdm *s, double *d, double *f)
{
    lag_print = 0;
    lag_d[0] = lag_d[1] = lag_f[0] = lag_f[1] = 0.0;

    if(p->S64>0.0)
    {
    susp_apply(p,pgc,s);
    qb_clear(p,s);

    if(p->S86>0)
    bend_K(p,pgc,s);
    }

    if(p->S95==2)
    spiral_apply(p,pgc,s);

    lag_print = 1;

    for(int k=0; k<2; ++k)
    {
    d[k] = MAX(d[k],lag_d[k]);
    f[k] = MAX(f[k],lag_f[k]);
    }
}
