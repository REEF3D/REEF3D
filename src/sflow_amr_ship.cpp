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

#include"sflow_amr_ship.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

sflow_amr_ship::sflow_amr_ship(lexer *p) : press(p), draft(p)
{
    for(int n=0; n<6; ++n)
    u[n]=0.0;
    c[0]=c[1]=c[2]=0.0;
    zref=0.0;

    // RK3 (A 210 3), as sixdof_obj
    alpha[0] = 1.0;
    alpha[1] = 0.25;
    alpha[2] = 2.0/3.0;
}

sflow_amr_ship::~sflow_amr_ship()
{
}

// solid heaviside of the body on the patch, as sixdof_obj::Hsolidface_2D
double sflow_amr_ship::Hsolid(lexer *p, slice &fs)
{
    double psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DYN[JP]);
    double phival_fb = fs(i,j);

    if(-phival_fb > psi)
    return 1.0;

    if(-phival_fb < -psi)
    return 0.0;

    return 0.5*(1.0 + -phival_fb/psi + (1.0/PI)*sin((PI*-phival_fb)/psi));
}

// X 10 2: direct forcing of the body velocity (the body is moved on level 0)
void sflow_amr_ship::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, int iter, slice &fs, slice &P, slice &Q, slice &w,
                                 slice &fx, slice &fy, slice &eta, bool finalize)
{
    if(p->X10!=2)
    return;

    double H, uf, vf;

    SLICELOOP4
    {
        uf = u[0] + u[4]*(zref - c[2]) - u[5]*(p->YP[JP] - c[1]);
        vf = u[1] + u[5]*(p->XP[IP] - c[0]) - u[3]*(zref - c[2]);

        H = Hsolid(p,fs);
        fx(i,j) += H*(uf - P(i,j))/(alpha[iter]*p->dt);
        fy(i,j) += H*(vf - Q(i,j))/(alpha[iter]*p->dt);
    }
}

// X 10 3: pressure of the moving body, as sixdof_sflow::isource2D/jsource2D
void sflow_amr_ship::isource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double dfdx;

    if(p->X10==3)
    SLICELOOP4
    {
    dfdx = (press(i+1,j)-press(i-1,j))/(p->DXP[IP]+p->DXP[IM1]);
    b->F(i,j) += b->WL(i,j)*dfdx/p->W1;
    }

    SLICELOOP4
    b->test(i,j) = press(i,j);
}

void sflow_amr_ship::jsource2D(lexer *p, fdm2D *b, ghostcell *pgc)
{
    double dfdy;

    if(p->X10==3)
    SLICELOOP4
    {
    dfdy = (press(i,j+1)-press(i,j-1))/(p->DYP[JP]+p->DYP[JM1]);
    b->G(i,j) += b->WL(i,j)*dfdy/p->W1;
    }
}

// =====================================================================================
//  sflow_amr: the moving body on the patches and the refinement zone around it
// =====================================================================================

#include"sflow_amr.h"
#include"6DOF_sflow.h"
#include"6DOF_obj.h"
#include<algorithm>
#include<mpi.h>
#include<cmath>

namespace
{
inline int fsh_s(int a, int l) { return (a>=0) ? (a>>l) : -(((-a)-1)>>l)-1; }
}

// the level-0 body level set at the centre of patch cell (ii,jj): bilinear from the level-0 cells
double sflow_amr::fs0_at(sflow_amr_patch &c, int ii, int jj)
{
    lexer *pp = c.pp;
    const int l = c.lev;
    const double x = pp->XP[ii+marge];
    const double y = pp->YP[jj+marge];

    const int ilo = p0->imin, ihi = p0->imin+p0->imax-1;
    const int jlo = p0->jmin, jhi = p0->jmin+p0->jmax-1;

    int ka = fsh_s(ii-EXT+c.I0,l) - O0i;
    int la = fsh_s(jj-EXT+c.J0,l) - O0j;
    if(x < p0->XP[ka+marge]) --ka;
    if(y < p0->YP[la+marge]) --la;
    int kb = ka+1, lb = la+1;
    ka = MAX(MIN(ka,ihi),ilo); kb = MAX(MIN(kb,ihi),ilo);
    la = MAX(MIN(la,jhi),jlo); lb = MAX(MIN(lb,jhi),jlo);

    double wx = (kb==ka) ? 0.0 : (x-p0->XP[ka+marge])/(p0->XP[kb+marge]-p0->XP[ka+marge]);
    double wy = (lb==la) ? 0.0 : (y-p0->YP[la+marge])/(p0->YP[lb+marge]-p0->YP[la+marge]);
    wx = MAX(MIN(wx,1.0),0.0);
    wy = MAX(MIN(wy,1.0),0.0);

    double fsv = 1.0e20;
    for(int nb=0; nb<ship6->objects(); ++nb)
    {
        slice &f = ship6->object(nb)->amr_fs();
        double v = (1.0-wx)*(1.0-wy)*f(ka,la) + wx*(1.0-wy)*f(kb,la) + (1.0-wx)*wy*f(ka,lb) + wx*wy*f(kb,lb);
        fsv = MIN(fsv,v);
    }
    return fsv;
}

// body fields on one patch: level set, kinematics and (X 10 3) the surface pressure
void sflow_amr::ship_fields(sflow_amr_patch &c, bool withpress)
{
    lexer *pp = c.pp;
    fdm2D *bp = c.b;
    sflow_amr_ship *s = c.pship;
    sixdof_obj *o0 = ship6->object(0);

    const int i0 = pp->imin, i1 = pp->imin+pp->imax;
    const int j0 = pp->jmin, j1 = pp->jmin+pp->jmax;

    // (G 7 1: linear in time within the level-0 step, ship_level)
    const double th = ship_th;
    for(int n=0; n<6; ++n)
    s->u[n] = (th<0.0) ? o0->amr_u(n) : (1.0-th)*sh_old[0].u[n] + th*o0->amr_u(n);
    for(int n=0; n<3; ++n)
    s->c[n] = (th<0.0) ? o0->amr_c(n) : (1.0-th)*sh_old[0].c[n] + th*o0->amr_c(n);
    s->zref = p0->ZP[marge];

    const double psi = 1.0e-8*p0->DXM;
    const double *xa = &pp->XP[i0+marge], *ya = &pp->YP[j0+marge];
    const int nx = i1-i0, ny = j1-j0;
    const double hx = (xa[nx-1]-xa[0])/double(nx-1), hy = (ya[ny-1]-ya[0])/double(ny-1);

    // draft: vertical rays through the patch cell centres, as sixdof_obj::ray_cast_2D_z
    {
        for(int ii=i0; ii<i1; ++ii)
        for(int jj=j0; jj<j1; ++jj)
        s->draft(ii,jj) = 0.0;

        for(int nb=0; nb<ship6->objects(); ++nb)
        {
            sixdof_obj *o = ship6->object(nb);
            double **tx = o->amr_tri(0), **ty = o->amr_tri(1), **tz = o->amr_tri(2);

            // patch away from the hull: no draft
            const reefamr_zone &z = zones[nb];
            if(z.bx1<xa[0]-psi || z.bx0>xa[nx-1]+psi || z.by1<ya[0]-psi || z.by0>ya[ny-1]+psi)
            continue;

            for(int n : ztri[nb])
            {
                double Ax=tx[n][0], Ay=ty[n][0], Az=tz[n][0];
                double Bx=tx[n][1], By=ty[n][1], Bz=tz[n][1];
                double Cx=tx[n][2], Cy=ty[n][2], Cz=tz[n][2];

                double xs = MIN3(Ax,Bx,Cx), xe = MAX3(Ax,Bx,Cx);
                double ys = MIN3(Ay,By,Cy), ye = MAX3(Ay,By,Cy);
                if(xe<xa[0]-psi || xs>xa[nx-1]+psi || ye<ya[0]-psi || ys>ya[ny-1]+psi)
                continue;

                double d = (By-Cy)*(Ax-Cx) + (Cx-Bx)*(Ay-Cy);
                if(fabs(d)<1.0e-30)
                continue;

                // first cell centre >= lo and first > hi (guess from the mean spacing, then corrected)
                auto first_ge = [](const double *c, int n, double h, double v)
                {
                    int k = MAX(MIN(int(ceil((v-c[0])/h)),n),0);
                    while(k>0 && c[k-1]>=v) --k;
                    while(k<n && c[k]<v) ++k;
                    return k;
                };
                auto first_gt = [](const double *c, int n, double h, double v)
                {
                    int k = MAX(MIN(int(floor((v-c[0])/h))+1,n),0);
                    while(k>0 && c[k-1]>v) --k;
                    while(k<n && c[k]<=v) ++k;
                    return k;
                };
                int a0 = first_ge(xa,nx,hx,xs);
                int a1 = first_gt(xa,nx,hx,xe+2.0*psi);
                int b0i = first_ge(ya,ny,hy,ys-2.0*psi);
                int b1i = first_gt(ya,ny,hy,ye);
                if(a0>=a1 || b0i>=b1i)
                continue;

                for(int a=a0; a<a1; ++a)
                for(int e=b0i; e<b1i; ++e)
                {
                    double Px = xa[a]-psi, Py = ya[e]+psi;
                    double u = ((By-Cy)*(Px-Cx) + (Cx-Bx)*(Py-Cy))/d;
                    double v = ((Cy-Ay)*(Px-Cx) + (Ax-Cx)*(Py-Cy))/d;
                    double w = 1.0-u-v;

                    if(u>0.0 && v>0.0 && w>0.0)
                    {
                        double Rz = u*Az + v*Bz + w*Cz;
                        if(p0->wd-Rz>0.0)
                        s->draft(a+i0,e+j0) = MAX(p0->wd-Rz,s->draft(a+i0,e+j0));
                    }
                }
            }
        }
    }

    // level set: from level 0 away from the hull.  Near the hull it is rebuilt from the
    // waterplane of the patch cells (draft > 0), with the interface half way between inside
    // and outside centres: interpolated from level 0 the waterline is blunted at a sharp stem
    // or transom and pulsates as the body moves through the coarse cells, which sends out
    // short waves on the fine grid
    for(int ii=i0; ii<i1; ++ii)
    for(int jj=j0; jj<j1; ++jj)
    bp->fs(ii,jj) = fs0_at(c,ii,jj);

    {
        const int R = 3;
        const double h = MIN(hx,hy);
        const double w = (R+1)*MAX(hx,hy);

        for(int nb=0; nb<(int)zones.size(); ++nb)
        {
            const reefamr_zone &z = zones[nb];
            if(z.bx1<xa[0]-w || z.bx0>xa[nx-1]+w || z.by1<ya[0]-w || z.by0>ya[ny-1]+w)
            continue;

            for(int a=0; a<nx; ++a)
            for(int e=0; e<ny; ++e)
            {
                if(xa[a]<z.bx0-w || xa[a]>z.bx1+w || ya[e]<z.by0-w || ya[e]>z.by1+w)
                continue;

                const bool in = s->draft(a+i0,e+j0)>0.0;
                double dmin = 1.0e20;
                for(int da=MAX(a-R,0); da<=MIN(a+R,nx-1); ++da)
                for(int de=MAX(e-R,0); de<=MIN(e+R,ny-1); ++de)
                if((s->draft(da+i0,de+j0)>0.0)!=in)
                dmin = MIN(dmin, sqrt(pow(xa[da]-xa[a],2.0) + pow(ya[de]-ya[e],2.0)));

                double val;
                if(dmin<1.0e19)
                val = dmin - 0.5*h;
                else
                val = MAX(fabs(bp->fs(a+i0,e+j0)),(R+0.5)*h);

                bp->fs(a+i0,e+j0) = in ? -val : val;
            }
        }
    }

    if(!withpress)
    return;

    const double ramp = o0->amr_ramp_draft(p0);

    // bilinear value of a patch field (clamped to the array)
    auto ipol = [&](slice &f, double x, double y)
    {
        const double *xa = &pp->XP[i0+marge], *ya = &pp->YP[j0+marge];
        const int nx = i1-i0, ny = j1-j0;
        int a = int(upper_bound(xa,xa+nx,x) - xa) - 1;
        int e = int(upper_bound(ya,ya+ny,y) - ya) - 1;
        a = MAX(MIN(a,nx-2),0);
        e = MAX(MIN(e,ny-2),0);
        double wx = MAX(MIN((x-xa[a])/(xa[a+1]-xa[a]),1.0),0.0);
        double wy = MAX(MIN((y-ya[e])/(ya[e+1]-ya[e]),1.0),0.0);
        return (1.0-wx)*(1.0-wy)*f(a+i0,e+j0) + wx*(1.0-wy)*f(a+1+i0,e+j0)
             + (1.0-wx)*wy*f(a+i0,e+1+j0) + wx*wy*f(a+1+i0,e+1+j0);
    };

    // surface pressure, as sixdof_obj::updateForcing_box/_oned/_stl
    for(int ii=i0; ii<i1; ++ii)
    for(int jj=j0; jj<j1; ++jj)
    {
        i=ii; j=jj;
        double H = s->Hsolid(pp,bp->fs);
        double xpos = pp->XP[IP] - p0->xg;
        double ypos = pp->YP[JP] - p0->yg;
        double pr = 0.0;

        if(p0->X400==2)
        {
            double Ls = p0->X110_xe[0] - p0->X110_xs[0];
            double Bs = p0->X110_ye[0] - p0->X110_ys[0];
            if(xpos<=Ls/2.0 && xpos>=-Ls/2.0 && ypos<=Bs/2.0 && ypos>=-Bs/2.0)
            pr = -H*p0->X401_p0*(1.0 - p0->X401_cl*pow(xpos/Ls,4.0))*(1.0 - p0->X401_cb*pow(ypos/Bs,2.0))
                   *exp(-p0->X401_a*pow(ypos/Bs,2.0))*ramp;
        }

        if(p0->X400==3)
        pr = p0->X401_p0*exp(-pow(xpos/p0->X401_a,2));

        if(p0->X400==10)
        {
            double etaval = 0.0;
            if(p0->X410==1 && ii>i0 && ii<i1-1 && jj>j0 && jj<j1-1)
            {
                double dfdx = (bp->fs(ii+1,jj) - bp->fs(ii-1,jj))/(pp->DXP[IP] + pp->DXP[IM1]);
                double dfdy = (bp->fs(ii,jj+1) - bp->fs(ii,jj-1))/(pp->DYP[JP] + pp->DYP[JM1]);
                double dnorm = sqrt(dfdx*dfdx + dfdy*dfdy);
                double nx = dfdx/(dnorm>1.0e-20?dnorm:1.0e20);
                double ny = dfdy/(dnorm>1.0e-20?dnorm:1.0e20);
                double xc = pp->XP[IP] + p0->X41*nx*pp->DXN[IP];
                double yc = pp->YP[JP] + p0->X41*ny*pp->DYN[JP];
                double fbval = ipol(bp->fs,xc,yc);
                if(fbval>-0.6*(1.0/2.0)*(pp->DXN[IP] + pp->DYN[JP]))
                etaval = ipol(bp->eta,xc,yc);
            }
            pr = -H*fabs(p0->W22)*p0->W1*(s->draft(ii,jj)+etaval)*ramp;
        }

        s->press(ii,jj) = pr;
    }
}

void sflow_amr::ship_patches(bool withpress)
{
    if(shipmode==0 || P.empty())
    return;

    double t0 = MPI_Wtime();

    zone_setup(p0);

    for(auto c : P)
    ship_fields(*SP(c),withpress);

    tm[8] += MPI_Wtime()-t0;
}

// G 7 1: the level-0 body now (the triangles, the level set, u, c)
void sflow_amr::ship_save(vector<shipsave> &S)
{
    const int nb = ship6->objects();
    S.resize(nb);
    for(int k=0; k<nb; ++k)
    {
        sixdof_obj *o = ship6->object(k);
        const int nt = o->amr_tricount();
        for(int d=0; d<3; ++d)
        {
            double **t = o->amr_tri(d);
            S[k].t[d].resize(3*(size_t)nt);
            for(int n=0; n<nt; ++n)
            for(int v=0; v<3; ++v)
            S[k].t[d][3*n+v] = t[n][v];
        }
        slice &f = ship6->object(k)->amr_fs();
        S[k].fs.assign(f.data(),f.data()+(size_t)p0->imax*p0->jmax);
        for(int n=0; n<6; ++n)
        S[k].u[n] = o->amr_u(n);
        for(int n=0; n<3; ++n)
        S[k].c[n] = o->amr_c(n);
    }
}

// into the level-0 objects: A, or (1-th) A + th B
void sflow_amr::ship_put(const vector<shipsave> &A, const vector<shipsave> *B, double th)
{
    for(int k=0; k<(int)A.size(); ++k)
    {
        sixdof_obj *o = ship6->object(k);
        const int nt = o->amr_tricount();
        for(int d=0; d<3; ++d)
        {
            double **t = o->amr_tri(d);
            for(int n=0; n<nt; ++n)
            for(int v=0; v<3; ++v)
            t[n][v] = (B==nullptr) ? A[k].t[d][3*n+v] : (1.0-th)*A[k].t[d][3*n+v] + th*(*B)[k].t[d][3*n+v];
        }
        slice &f = ship6->object(k)->amr_fs();
        const size_t nc = A[k].fs.size();
        for(size_t m=0; m<nc; ++m)
        f.data()[m] = (B==nullptr) ? A[k].fs[m] : (1.0-th)*A[k].fs[m] + th*(*B)[k].fs[m];
    }
}

// G 7 1: the body on the level-l patches at time t (in the level-0 step that started at sub_t0):
// the level-0 body linear in time between the start and the end of that step
void sflow_amr::ship_level(int l, double t, bool withpress)
{
    if(shipmode==0 || lev[l].empty())
    return;

    double t0 = MPI_Wtime();

    double th = (t - sub_t0)/sub_dt0;
    th = MAX(MIN(th,1.0),0.0);
    const bool ti = (th<1.0-1.0e-12) && !sh_old.empty();

    const double simt = p0->simtime;
    p0->simtime = t;

    if(ti)
    {
        ship_save(sh_new);
        ship_put(sh_old,&sh_new,th);
        ship_th = th;
    }

    zone_setup(p0);
    for(int n : lev[l])
    ship_fields(*SP(n),withpress);

    if(ti)
    {
        ship_put(sh_new,nullptr,0.0);
        ship_th = -1.0;
        zone_setup(p0);
    }

    p0->simtime = simt;

    tm[8] += MPI_Wtime()-t0;
}
