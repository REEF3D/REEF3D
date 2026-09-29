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

// triangles below the still water level and the refinement zone of every body (A 278, A 279)
void sflow_amr::ship_setup(lexer *p)
{
    const double psi = 1.0e-8*p->DXM;
    const int nbody = ship6->objects();

    shiptri.assign(nbody,vector<int>());
    zones.assign(nbody,shipzone());

    if((int)ship_x0.size()!=nbody)
    {
        ship_x0.resize(nbody);
        ship_y0.resize(nbody);
        for(int nb=0; nb<nbody; ++nb)
        {
            ship_x0[nb] = ship6->object(nb)->amr_c(0);
            ship_y0[nb] = ship6->object(nb)->amr_c(1);
        }
    }

    for(int nb=0; nb<nbody; ++nb)
    {
        sixdof_obj *o = ship6->object(nb);
        double **tx = o->amr_tri(0), **ty = o->amr_tri(1), **tz = o->amr_tri(2);

        shipzone &z = zones[nb];
        z.cx = o->amr_c(0);
        z.cy = o->amr_c(1);

        // heading: direction of motion, the yaw angle for a body at rest
        double ux = o->amr_u(0), uy = o->amr_u(1);
        double sp = sqrt(ux*ux + uy*uy);
        if(sp>1.0e-6)
        {
            z.ex = ux/sp;
            z.ey = uy/sp;
        }
        else
        {
            z.ex = cos(p->psi_fb);
            z.ey = sin(p->psi_fb);
        }

        z.smin = z.nmin = 1.0e20;
        z.smax = z.nmax = -1.0e20;
        z.bx0 = z.by0 = 1.0e20;
        z.bx1 = z.by1 = -1.0e20;

        for(int n=0; n<o->amr_tricount(); ++n)
        {
            if(tz[n][0]>p->wd+psi && tz[n][1]>p->wd+psi && tz[n][2]>p->wd+psi)
            continue;

            shiptri[nb].push_back(n);

            for(int q=0; q<3; ++q)
            {
                z.bx0 = MIN(z.bx0,tx[n][q]); z.bx1 = MAX(z.bx1,tx[n][q]);
                z.by0 = MIN(z.by0,ty[n][q]); z.by1 = MAX(z.by1,ty[n][q]);
            }

            for(int q=0; q<3; ++q)
            if(tz[n][q]<=p->wd+psi)
            {
                double dx = tx[n][q]-z.cx, dy = ty[n][q]-z.cy;
                double s = dx*z.ex + dy*z.ey;
                double r = -dx*z.ey + dy*z.ex;
                z.smin = MIN(z.smin,s); z.smax = MAX(z.smax,s);
                z.nmin = MIN(z.nmin,r); z.nmax = MAX(z.nmax,r);
            }
        }

        // the zone reaches ahead of the bow by the distance travelled until the next regrid
        z.sfront = z.smax + p->A278_r + 1.5*sp*p->dt*MAX(regrid_int,1);

        // the wake wedge reaches back to where the bow has been
        double trav = sqrt(pow(z.cx-ship_x0[nb],2.0) + pow(z.cy-ship_y0[nb],2.0));
        z.wake = MIN(p->A279_L, (z.smax-z.smin) + trav + p->A278_r);

        if(z.smin>z.smax)       // no part of the body in the water
        z.smin = z.smax = z.nmin = z.nmax = z.sfront = z.wake = 0.0;
    }
}

bool sflow_amr::ship_zone(double x, double y)
{
    const double r = p0->A278_r;
    const double ta = tan(p0->A279_a*PI/180.0);

    for(auto &z : zones)
    {
        if(z.smin>=z.smax)
        continue;

        double dx = x-z.cx, dy = y-z.cy;
        double s = dx*z.ex + dy*z.ey;
        double n = -dx*z.ey + dy*z.ex;

        // hull
        if(s>=z.smin-r && s<=z.sfront && n>=z.nmin-r && n<=z.nmax+r)
        return true;

        // wake: wedge from the bow with half angle A 279 a, length A 279 L at most (only where the bow has been)
        if(z.wake>0.0)
        {
            double d = z.smax - s;
            if(d>=0.0 && d<=z.wake)
            {
                double half = 0.5*(z.nmax-z.nmin) + r + d*ta;
                if(fabs(n - 0.5*(z.nmin+z.nmax))<=half)
                return true;
            }
        }
    }
    return false;
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

    for(int n=0; n<6; ++n)
    s->u[n] = o0->amr_u(n);
    for(int n=0; n<3; ++n)
    s->c[n] = o0->amr_c(n);
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
            const shipzone &z = zones[nb];
            if(z.bx1<xa[0]-psi || z.bx0>xa[nx-1]+psi || z.by1<ya[0]-psi || z.by0>ya[ny-1]+psi)
            continue;

            for(int n : shiptri[nb])
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
            const shipzone &z = zones[nb];
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

    ship_setup(p0);

    for(auto c : P)
    ship_fields(*c,withpress);

    tm[8] += MPI_Wtime()-t0;
}
