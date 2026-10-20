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

// Bend morphology: two corrections of the bed-load direction (both off by default).
//
// S 86 1, transverse bed slope effect: gravity turns the bed load down the bed slope,
//     q = qb ( s - grad(z)/f(theta) ),  f(theta) = A_sh theta^B_sh (d/h)^C_sh   (S 87, S 88, S 89;
//     Talmon, van Mierlo & Struiksma 1995: 9, 0.5, 0.3). The part -qb/f dz/dn is a face flux of its own
//     (central, the mean of qb/f of the two cells, dz/dn between the cell centres): conservative,
//     diffusive, the same on both sides of a face (the faces between grids of the mesh refinement and
//     between blocks see one value).
// S 95 1, spiral-flow correction (depth-averaged SFLOW): the helical flow of a bend turns the near-bed
//     flow towards the centre of curvature of the streamlines by delta,
//     tan(delta) = alpha_I h / R_s,  alpha_I = E_spir 2/kappa^2 (1 - sqrt(g)/(kappa C))   (S 96 E_spir;
//     Engelund 1974, Struiksma et al. 1985), 1/R_s the signed curvature of the streamlines of the
//     depth-averaged velocity (positive turning left), C = sqrt(g)|U|/u* from the bed shear. The
//     transport direction of every face is turned by delta (the mean of its two cells).
// NHFLOW resolves the secondary flow itself (the bed velocity of the lowest layer): S 95 is refused there.

// the cell fields of the bend terms in the sediment cells, the ghost cells exchanged; on the mesh refinement
// and the multi-block grids the cells around a grid take them from their sources (amr_field) as qb, so that
// the two grids of a face compute the same flux
void sediment_exner::bend_prep(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    if(p->S86>0)
    {
        // the bed of the slope flux and the cells whose bed is held (the inflow guard S 78): the cells around
        // a grid of the mesh refinement take them from their source (their flow bed can lag behind after a
        // regrid, their own guard does not apply)
        SLICEBASELOOP
        bK(i,j) = bz(i,j) = bg(i,j) = 0.0;

        SLICELOOP4
        {
        bz(i,j) = s->bedzh(i,j);
        bg(i,j) = (s->guard(i,j)>0.0) ? 1.0 : 0.0;
        }

        bend_K(p,pgc,s);

        pgc->gcsl_start4(p,bz,1);
        pgc->gcsl_start4(p,bg,1);
    }

    if(p->S95>0)
    {
        const double kappa = 0.4;
        const double g = fabs(p->W22);
        double u,v;

        SLICEBASELOOP
        {
            if(p->geo_curv==1)
            sediment_curv::cellvel(p,s,i,j,u,v);
            else
            {
            u = 0.5*(s->P(i,j)+s->P(i-1,j));
            v = 0.5*(s->Q(i,j)+s->Q(i,j-1));
            }
        bu(i,j) = u;
        bv(i,j) = v;
        btd(i,j) = 0.0;
        }

        pgc->gcsl_start4(p,bu,1);
        pgc->gcsl_start4(p,bv,1);

        // gradients: curvilinear the chain rule of sediment_curv::bedgrad, Cartesian central (one-sided
        // next to cells without a bed)
        auto grad = [&](slice &f, double &fx, double &fy)
        {
            if(p->geo_curv==1)
            {
            sediment_curv::bedgrad(p,f,i,j,0,fx,fy);
            return;
            }

            fx = (f(i+1,j)-f(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
            if(p->DFBED[Im1J]<0)
            fx = (f(i+1,j)-f(i,j))/p->DXP[IP];
            if(p->DFBED[Ip1J]<0)
            fx = (f(i,j)-f(i-1,j))/p->DXP[IM1];
            if(p->DFBED[Im1J]<0 && p->DFBED[Ip1J]<0)
            fx = 0.0;

            fy = 0.0;
            if(p->j_dir==1)
            {
            fy = (f(i,j+1)-f(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
            if(p->DFBED[IJm1]<0)
            fy = (f(i,j+1)-f(i,j))/p->DYP[JP];
            if(p->DFBED[IJp1]<0)
            fy = (f(i,j)-f(i,j-1))/p->DYP[JM1];
            if(p->DFBED[IJm1]<0 && p->DFBED[IJp1]<0)
            fy = 0.0;
            }
        };

        double ux,uy,vx,vy;

        SEDSLICELOOP
        {
            u = bu(i,j);
            v = bv(i,j);
            const double U = sqrt(u*u + v*v);

            if(U<1.0e-4)
            continue;

            grad(bu,ux,uy);
            grad(bv,vx,vy);

            // signed curvature of the streamlines: u x (u.grad u)/|u|^3
            const double curv = (u*u*vx + u*v*vy - u*v*ux - v*v*uy)/(U*U*U);

            const double ustar = s->shearvel_eff(i,j);
            const double C = (ustar>1.0e-10) ? sqrt(g)*U/ustar : 1.0e10;
            const double alpha = MAX(0.0, p->S96*2.0/(kappa*kappa)*(1.0 - sqrt(g)/(kappa*C)));
            const double h = MAX(s->waterlevel(i,j),0.0);

            // S 95 2: the equilibrium intensity, relaxed below (sediment_exner_lag.cpp)
            if(p->S95==2)
            Ie(i,j) = h*U*curv;
            else
            btd(i,j) = MAX(-1.0, MIN(1.0, alpha*h*curv));
        }

        if(p->S95==2)
        {
        SEDSLICELOOP
        if(sqrt(bu(i,j)*bu(i,j) + bv(i,j)*bv(i,j))<1.0e-4)
        Ie(i,j) = 0.0;
        spiral_lag(p,pgc,s);
        }

        pgc->gcsl_start4(p,btd,1);

    }
}

// S 86: qb/f(theta) of the slope flux and the largest rate of its explicit diffusion (the time step); again after
// the lags of the suspended load moved qb (lag_refresh)
void sediment_exner::bend_K(lexer *p, ghostcell *pgc, sediment_fdm *s)
{
    const double d = (s->pmix!=nullptr) ? s->dk : p->S20;

    SEDSLICELOOP
    {
        const double theta = MAX(s->shields_eff(i,j),1.0e-6);
        const double h = MAX(s->waterlevel(i,j),1.0e-4);
        const double f = p->S87*pow(theta,p->S88)*pow(d/h,p->S89);
        bK(i,j) = (f>1.0e-10) ? s->qb(i,j)/f : 0.0;
    }

    bend_rate = 0.0;
    SEDSLICELOOP
    {
        double w1, w2;
        if(p->geo_curv==1)
        {
        w1 = sediment_curv::width(p,1,i,j);
        w2 = sediment_curv::width(p,2,i,j);
        }
        else
        {
        w1 = p->DXN[IP];
        w2 = p->DYN[JP];
        }
        const double r = 1.0/(w1*w1) + ((p->j_dir==1 && p->gknoy>1) ? 1.0/(w2*w2) : 0.0);
        bend_rate = MAX(bend_rate, bK(i,j)*r);
    }

    pgc->gcsl_start4(p,bK,1);
}

// the transport direction un/|u| of a face turned to the left by delta (tan(delta) = td), ut the component
// along the face (to the left of the normal)
double sediment_exner::bend_turn(double un, double ut, double mag, double td)
{
    if(fabs(mag)<=1.0e-10)
    return 0.0;

    if(td==0.0)
    return un/fabs(mag);

    return (un - td*ut)/(sqrt(1.0 + td*td)*fabs(mag));
}

// a cell with an erodible bed inside the domain (the cells around a patch or block: from the flags of their
// source, outside the domain closed)
bool sediment_exner::bend_open(lexer *p, int a, int b)
{
    const int c = (a-p->imin)*p->jmax + (b-p->jmin);
    return a+p->origin_i>=0 && a+p->origin_i<p->gknox && b+p->origin_j>=0 && b+p->origin_j<p->gknoy
        && p->flagslice4[c]>0 && p->DFBED[c]>0;
}

// tan(delta) of the face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (ii,jj): the mean of its two cells; 0 at
// closed faces (walls, the domain edge), where the turn would open the face, and where the two cells cancel
// to round-off (a symmetry line: the upwind side of S 31 2 / S 32 1 must not follow the sign of round-off)
double sediment_exner::bend_td(lexer *p, int dir, int ii, int jj)
{
    if(p_S95==0)
    return 0.0;

    const int ia = ii + (dir==1?1:0);
    const int ja = jj + (dir==2?1:0);

    if(!bend_open(p,ii,jj) || !bend_open(p,ia,ja))
    return 0.0;

    const double a = btd(ii,jj), b = btd(ia,ja);
    const double t = 0.5*(a + b);

    return (fabs(t) > 1.0e-12*(fabs(a) + fabs(b))) ? t : 0.0;
}

// the slope flux per unit lengththrough the face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (ii,jj), positive
// towards +i / +j: -(qb/f)_face dz/dn; closed to cells without an erodible bed and at the domain edges
double sediment_exner::bend_flux(lexer *p, sediment_fdm *s, int dir, int ii, int jj)
{
    if(p->S86==0)
    return 0.0;

    const int ia = ii + (dir==1?1:0);
    const int ja = jj + (dir==2?1:0);

    if(dir==2 && !(p->j_dir==1 && p->gknoy>1))
    return 0.0;

    // closed also to the cells whose bed is held (the inflow guard S 78: no bed change there)
    if(!bend_open(p,ii,jj) || !bend_open(p,ia,ja) || bg(ii,jj)<0.5 || bg(ia,ja)<0.5)
    return 0.0;

    double dn;

    if(p->geo_curv==1)
    {
        double nx,ny,len;
        sediment_curv::normal(p,dir,ii,jj,nx,ny,len);
        const int c0 = sediment_curv::q(p,ii,jj), c1 = sediment_curv::q(p,ia,ja);
        dn = (p->XC2D[c1]-p->XC2D[c0])*nx + (p->YC2D[c1]-p->YC2D[c0])*ny;
    }
    else
    dn = (dir==1) ? p->DXP[ii+marge] : p->DYP[jj+marge];

    if(fabs(dn)<1.0e-12)
    return 0.0;

    return -0.5*(bK(ii,jj) + bK(ia,ja))*(bz(ia,ja) - bz(ii,jj))/dn;
}

// the divergence of the slope fluxes of cell (ii,jj)
double sediment_exner::bend_div(lexer *p, sediment_fdm *s, int ii, int jj)
{
    if(p->S86==0)
    return 0.0;

    const double fx2 = bend_flux(p,s,1,ii,jj), fx1 = bend_flux(p,s,1,ii-1,jj);
    const double fy2 = bend_flux(p,s,2,ii,jj), fy1 = bend_flux(p,s,2,ii,jj-1);

    if(p->geo_curv==1)
    {
        const int c = sediment_curv::q(p,ii,jj);
        return (fx2*p->LEN1[c] - fx1*p->LEN1[sediment_curv::q(p,ii-1,jj)]
              + fy2*p->LEN2[c] - fy1*p->LEN2[sediment_curv::q(p,ii,jj-1)])/p->AREA2D[c];
    }

    return (fx2-fx1)/p->DXN[ii+marge] + (fy2-fy1)/p->DYN[jj+marge];
}
