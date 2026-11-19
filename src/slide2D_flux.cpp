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

#include"slide2D.h"
#include"lexer.h"
#include"fdm2D.h"
#include<cmath>

// Fluxes of the slide layer: MUSCL (minmod) of h, eta_s = z_b + h, u, v; hydrostatic reconstruction at the faces
// (Audusse et al. 2004, second-order form with the centred bed term), HLL with the reduced gravity g' of the face.
// The x faces i+1/2 are stored at i, the y faces j+1/2 at j.

namespace
{
struct face_state
{
    double h, z, u, v;     // depth, bed, normal and tangential velocity at the face
};
}

double slide2D::minmod(double a, double b)
{
    if(a*b<=0.0)
    return 0.0;

    return fabs(a)<fabs(b) ? a : b;
}

double slide2D::vel(slice &H, slice &Q, int ii, int jj)
{
    return H(ii,jj)>hmin ? Q(ii,jj)/H(ii,jj) : 0.0;
}

void slide2D::rhs(lexer *p, slice &H, slice &QX, slice &QY)
{
    flux_x(p,H,QX,QY);

    if(p->j_dir==1)
    flux_y(p,H,QX,QY);

    const double ydir = p->j_dir==1 ? 1.0 : 0.0;

    SLICELOOP4
    {
        const double dx = p->DXN[IP];
        const double dy = p->DYN[JP];
        const double grc = gr(i,j);

        Rh(i,j) = -(FhX(i,j) - FhX(i-1,j))/dx;
        Rx(i,j) = -(FxLX(i,j) - FxRX(i-1,j))/dx - grc*0.5*(hE(i,j) + hW(i,j))*(zE(i,j) - zW(i,j))/dx;
        Ry(i,j) = -(FyX(i,j) - FyX(i-1,j))/dx;

        if(ydir>0.0)
        {
        Rh(i,j) += -(FhY(i,j) - FhY(i,j-1))/dy;
        Rx(i,j) += -(FxY(i,j) - FxY(i,j-1))/dy;
        Ry(i,j) += -(FyLY(i,j) - FyRY(i,j-1))/dy - grc*0.5*(hN(i,j) + hS(i,j))*(zN(i,j) - zS(i,j))/dy;
        }
    }

    sources(p,H,QX,QY);
}

void slide2D::flux_x(lexer *p, slice &H, slice &QX, slice &QY)
{
    auto fl = [p](int ii, int jj){ return p->flagslice4[(ii-p->imin)*p->jmax + (jj-p->jmin)]; };

    // reconstruction in x of the cell (ii,jj) at its east (side +1) or west (side -1) face
    auto recon = [&](int ii, int jj, int side)
    {
        face_state f;
        const double hc = H(ii,jj);
        const double zc = zb(ii,jj);
        const double uc = vel(H,QX,ii,jj);
        const double vc = vel(H,QY,ii,jj);

        f.h = hc;
        f.z = zc;
        f.u = uc;
        f.v = vc;

        if(hc<=hmin || p->A211==0)
        return f;

        const double dm = p->DXP[ii-1+marge];
        const double dp = p->DXP[ii+marge];
        const double hd = 0.5*side*p->DXN[ii+marge];

        const double hm = H(ii-1,jj), hp = H(ii+1,jj);
        const double em = hm + zb(ii-1,jj), ep = hp + zb(ii+1,jj), ec = hc + zc;

        const double sh = minmod((hp-hc)/dp, (hc-hm)/dm);
        const double se = minmod((ep-ec)/dp, (ec-em)/dm);

        f.h = MAX(hc + hd*sh, 0.0);
        f.z = ec + hd*se - f.h;

        // velocities: limited only between wet cells
        if(hm>hmin && hp>hmin)
        {
        f.u = uc + hd*minmod((vel(H,QX,ii+1,jj)-uc)/dp, (uc-vel(H,QX,ii-1,jj))/dm);
        f.v = vc + hd*minmod((vel(H,QY,ii+1,jj)-vc)/dp, (vc-vel(H,QY,ii-1,jj))/dm);
        }

        return f;
    };

    // the faces of each cell (for the centred bed term)
    for(i=-1; i<=p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    {
        const face_state e = recon(i,j,1);
        const face_state w = recon(i,j,-1);

        hE(i,j) = e.h;
        zE(i,j) = e.z;
        hW(i,j) = w.h;
        zW(i,j) = w.z;
    }

    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    {
        FhX(i,j) = FxLX(i,j) = FxRX(i,j) = FyX(i,j) = 0.0;

        const int fL = fl(i,j);
        const int fR = fl(i+1,j);

        if(fL<0 && fR<0)
        continue;

        face_state L = recon(i,j,1);
        face_state R = recon(i+1,j,-1);

        // walls: the mirror state
        if(fR<0)
        {
        R = L;
        R.u = -L.u;
        }

        if(fL<0)
        {
        L = R;
        L.u = -R.u;
        }

        if(L.h<=hmin && R.h<=hmin)
        continue;

        // hydrostatic reconstruction
        const double zf = MAX(L.z,R.z);
        const double hl = MAX(0.0, L.h + L.z - zf);
        const double hr = MAX(0.0, R.h + R.z - zf);

        const double gf = 0.5*((fL>0 ? gr(i,j) : gr(i+1,j)) + (fR>0 ? gr(i+1,j) : gr(i,j)));

        double fh=0.0, fx=0.0;

        // both cells at rest (held by the friction or dry): no transport across the face, the central pressure.
        // The HLL dissipation (hr - hl) moved material out of a deposit at rest and let its front creep.
        const bool rest = QX(i,j)==0.0 && QY(i,j)==0.0 && QX(i+1,j)==0.0 && QY(i+1,j)==0.0;

        if(rest)
        fx = 0.25*gf*(hl*hl + hr*hr);

        if(!rest && (hl>0.0 || hr>0.0))
        {
        const double cl = sqrt(gf*hl);
        const double cr = sqrt(gf*hr);
        const double sl = MIN(L.u - cl, R.u - cr);
        const double sr = MAX(L.u + cl, R.u + cr);

        const double FhL = hl*L.u,  FhR = hr*R.u;
        const double FxL = hl*L.u*L.u + 0.5*gf*hl*hl;
        const double FxR = hr*R.u*R.u + 0.5*gf*hr*hr;

            if(sl>=0.0)
            {
            fh = FhL;
            fx = FxL;
            }
            else if(sr<=0.0)
            {
            fh = FhR;
            fx = FxR;
            }
            else
            {
            fh = (sr*FhL - sl*FhR + sl*sr*(hr - hl))/(sr - sl);
            fx = (sr*FxL - sl*FxR + sl*sr*(hr*R.u - hl*L.u))/(sr - sl);
            }
        }

        FhX(i,j) = fh;
        FxLX(i,j) = fx + 0.5*gf*(L.h*L.h - hl*hl);
        FxRX(i,j) = fx + 0.5*gf*(R.h*R.h - hr*hr);
        FyX(i,j) = fh>0.0 ? fh*L.v : fh*R.v;
    }
}

void slide2D::flux_y(lexer *p, slice &H, slice &QX, slice &QY)
{
    auto fl = [p](int ii, int jj){ return p->flagslice4[(ii-p->imin)*p->jmax + (jj-p->jmin)]; };

    // reconstruction in y of the cell (ii,jj) at its north (side +1) or south (side -1) face: u normal, v tangential
    auto recon = [&](int ii, int jj, int side)
    {
        face_state f;
        const double hc = H(ii,jj);
        const double zc = zb(ii,jj);
        const double uc = vel(H,QY,ii,jj);
        const double vc = vel(H,QX,ii,jj);

        f.h = hc;
        f.z = zc;
        f.u = uc;
        f.v = vc;

        if(hc<=hmin || p->A211==0)
        return f;

        const double dm = p->DYP[jj-1+marge];
        const double dp = p->DYP[jj+marge];
        const double hd = 0.5*side*p->DYN[jj+marge];

        const double hm = H(ii,jj-1), hp = H(ii,jj+1);
        const double em = hm + zb(ii,jj-1), ep = hp + zb(ii,jj+1), ec = hc + zc;

        const double sh = minmod((hp-hc)/dp, (hc-hm)/dm);
        const double se = minmod((ep-ec)/dp, (ec-em)/dm);

        f.h = MAX(hc + hd*sh, 0.0);
        f.z = ec + hd*se - f.h;

        if(hm>hmin && hp>hmin)
        {
        f.u = uc + hd*minmod((vel(H,QY,ii,jj+1)-uc)/dp, (uc-vel(H,QY,ii,jj-1))/dm);
        f.v = vc + hd*minmod((vel(H,QX,ii,jj+1)-vc)/dp, (vc-vel(H,QX,ii,jj-1))/dm);
        }

        return f;
    };

    for(i=0; i<p->knox; ++i)
    for(j=-1; j<=p->knoy; ++j)
    {
        const face_state n = recon(i,j,1);
        const face_state s = recon(i,j,-1);

        hN(i,j) = n.h;
        zN(i,j) = n.z;
        hS(i,j) = s.h;
        zS(i,j) = s.z;
    }

    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    {
        FhY(i,j) = FxY(i,j) = FyLY(i,j) = FyRY(i,j) = 0.0;

        const int fL = fl(i,j);
        const int fR = fl(i,j+1);

        if(fL<0 && fR<0)
        continue;

        face_state L = recon(i,j,1);
        face_state R = recon(i,j+1,-1);

        if(fR<0)
        {
        R = L;
        R.u = -L.u;
        }

        if(fL<0)
        {
        L = R;
        L.u = -R.u;
        }

        if(L.h<=hmin && R.h<=hmin)
        continue;

        const double zf = MAX(L.z,R.z);
        const double hl = MAX(0.0, L.h + L.z - zf);
        const double hr = MAX(0.0, R.h + R.z - zf);

        const double gf = 0.5*((fL>0 ? gr(i,j) : gr(i,j+1)) + (fR>0 ? gr(i,j+1) : gr(i,j)));

        double fh=0.0, fy=0.0;

        const bool rest = QX(i,j)==0.0 && QY(i,j)==0.0 && QX(i,j+1)==0.0 && QY(i,j+1)==0.0;

        if(rest)
        fy = 0.25*gf*(hl*hl + hr*hr);

        if(!rest && (hl>0.0 || hr>0.0))
        {
        const double cl = sqrt(gf*hl);
        const double cr = sqrt(gf*hr);
        const double sl = MIN(L.u - cl, R.u - cr);
        const double sr = MAX(L.u + cl, R.u + cr);

        const double FhL = hl*L.u,  FhR = hr*R.u;
        const double FyL = hl*L.u*L.u + 0.5*gf*hl*hl;
        const double FyR = hr*R.u*R.u + 0.5*gf*hr*hr;

            if(sl>=0.0)
            {
            fh = FhL;
            fy = FyL;
            }
            else if(sr<=0.0)
            {
            fh = FhR;
            fy = FyR;
            }
            else
            {
            fh = (sr*FhL - sl*FhR + sl*sr*(hr - hl))/(sr - sl);
            fy = (sr*FyL - sl*FyR + sl*sr*(hr*R.u - hl*L.u))/(sr - sl);
            }
        }

        FhY(i,j) = fh;
        FyLY(i,j) = fy + 0.5*gf*(L.h*L.h - hl*hl);
        FyRY(i,j) = fy + 0.5*gf*(R.h*R.h - hr*hr);
        FxY(i,j) = fh>0.0 ? fh*L.v : fh*R.v;
    }
}

void slide2D::sources(lexer *p, slice &H, slice &QX, slice &QY)
{
    auto etaw = [&](int ii, int jj){ return (*weta)(ii,jj) + H(ii,jj) - hstep(ii,jj); };

    auto wet = [p,this](int ii, int jj){ return p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]==1 && (*wWL)(ii,jj)>wdry
                                          && p->flagslice4[(ii-p->imin)*p->jmax + (jj-p->jmin)]>0; };

    SLICELOOP4
    {
        const double hh = H(i,j);

        if(hh<=hmin)
        continue;

        // along-slope gravity (inclined reference plane, the slide alone)
        Rx(i,j) += p->W20*hh;
        Ry(i,j) += p->W21*hh*p->y_dir;

        // the water pressure on the slide surface: - r sub g h grad(eta_w); one-sided next to dry cells. The water column
        // is held over the sub-steps (the water moves in its own step): eta_w = eta^n + h_s - h_s^n, so a slide that thickens
        // lifts the water above it within the step (otherwise the frozen surface slope drove a two-layer mode unstable in
        // the long semi-implicit steps)
        const double sb = sub_cell(p,i,j);

        if(sb<=0.0)
        continue;

        const double ec = etaw(i,j);
        double gx = 0.0, gy = 0.0;

        const bool we = wet(i+1,j), ww = wet(i-1,j);

        if(we && ww)
        gx = (etaw(i+1,j) - etaw(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
        else if(we)
        gx = (etaw(i+1,j) - ec)/p->DXP[IP];
        else if(ww)
        gx = (ec - etaw(i-1,j))/p->DXP[IM1];

        if(p->j_dir==1)
        {
        const bool wn = wet(i,j+1), ws = wet(i,j-1);

        if(wn && ws)
        gy = (etaw(i,j+1) - etaw(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
        else if(wn)
        gy = (etaw(i,j+1) - ec)/p->DYP[JP];
        else if(ws)
        gy = (ec - etaw(i,j-1))/p->DYP[JM1];
        }

        Rx(i,j) -= r*sb*g*hh*gx;
        Ry(i,j) -= r*sb*g*hh*gy*p->y_dir;
    }
}
