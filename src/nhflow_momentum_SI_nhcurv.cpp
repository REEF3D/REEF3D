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

#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<cmath>

// ---------------------------------------------------------------------------------------------------------------------
// NHFLOW semi-implicit step with the non-hydrostatic pressure (A 520 > 0) on one curvilinear grid (G 50 2).
// The coupled eta-phi system of nhflow_momentum_SI_nh.cpp with the geometry of the curvilinear grid:
//  - the face terms of phi (nh_faces, nh_dfaces) along the face normal, over the normal distance of the two centres
//    (fdn); the layer discharges of the faces are the normal ones, so face_update and the divergence are unchanged.
//    The non-orthogonal part of the face gradient of phi is not included (as that of eta without A 539 1): exact on
//    orthogonal grids (also turned ones), first order in the skewness otherwise.
//  - the divergence of the rows with the face lengths over the cell area (nh_assemble_col: adddiv, as div)
//  - [u . grad z]_k and the bed condition: the gradient of the interface heights from the centres of the neighbours
//    along the two index directions (as cgrad), grad z = a Z_1 + b Z_2 with Z_1 = z_E - z_W, Z_2 = z_N - z_S. Then
//    u . grad z = (u . a) Z_1 + (u . b) Z_2: the constraint keeps its form with the differences along the index
//    directions over the centre distances L_1, L_2 as slopes (NSX, NSY, with the change of eta of the solve) and the
//    velocities (u . a) L_1, (u . b) L_2 in place of UH, VH (nh_curv_metric, nh_curv_vel: nUZx, nUZy)
//  - the phi face terms of a cell (UH, VH Cartesian): from the mean normals and the mean face terms of the open
//    faces of the two index directions, as the gravity term of cell_update (nh_cell_faces_curv)
// Second order in space (fourth order on a turned uniform grid, nhflow_momentum_SI_ho.cpp), the nonlinear terms with the
// second-order time levels of the Cartesian step (tl2; nh_curv_vel, nh_cell_faces_curv).
// ---------------------------------------------------------------------------------------------------------------------

void nhflow_momentum_SI::nh_curv_metric(lexer *p, int ie, int iw, bool oe, bool ow, int jn, int js, bool on, bool os)
{
    // cell (i,j): the vectors between the centres of its neighbours along the two index directions (ie, iw, jn, js:
    // the cells of the slopes of nh_prepare; an open boundary side: the ghost centre at the normal distance of the
    // face along its normal), their lengths L_1, L_2 and the velocity transform T
    auto lij = [&](int a, int b) { return (a-p->imin)*p->jmax + (b-p->jmin); };
    const double x0 = p->XC2D[lij(i,j)], y0 = p->YC2D[lij(i,j)];

    // the vector from the centre of (i,j) to the neighbour across face (fi,fj) of direction dir, sign s
    auto nb = [&](int dir, int a, int b, bool open, int fi, int fj, double s, double &vx, double &vy)
    {
        vx = vy = 0.0;

        if(a==i && b==j)
        return;

        if(!open)
        {
        vx = p->XC2D[lij(a,b)] - x0;
        vy = p->YC2D[lij(a,b)] - y0;
        return;
        }

        double nx,ny;
        fnorm(p,dir,fi,fj,nx,ny);
        const double dn = fdn(p,dir,fi,fj);
        vx = s*nx*dn;
        vy = s*ny*dn;
    };

    double ex,ey,wx_,wy_,nx_,ny_,sx,sy;
    nb(1,ie,j,oe,i,j,1.0,ex,ey);
    nb(1,iw,j,ow,i-1,j,-1.0,wx_,wy_);
    nb(2,i,jn,on,i,j,1.0,nx_,ny_);
    nb(2,i,js,os,i,j-1,-1.0,sx,sy);

    double e1x = ex - wx_, e1y = ey - wy_, e2x = nx_ - sx, e2y = ny_ - sy;
    double L1 = sqrt(e1x*e1x + e1y*e1y), L2 = sqrt(e2x*e2x + e2y*e2y);
    const bool s1 = L1>0.0, s2 = L2>0.0;

    // no neighbour along a direction: no slope along it (L = 0); the direction of its face normal for the transform
    if(!s1)
    fnorm(p,1,i,j,e1x,e1y);

    if(!s2)
    fnorm(p,2,i,j,e2x,e2y);

    const double det = e1x*e2y - e1y*e2x;
    const size_t q = (size_t)lij(i,j);

    nCL1[q] = s1 ? L1 : 0.0;
    nCL2[q] = s2 ? L2 : 0.0;

    if(fabs(det)>0.0)
    {
    // grad f = (a1, c1) F_1 + (b1, d1) F_2:  u . grad f = (u a1 + v c1) F_1 + (u b1 + v d1) F_2
    nCT[4*q]   =  e2y/det*nCL1[q];
    nCT[4*q+1] = -e2x/det*nCL1[q];
    nCT[4*q+2] = -e1y/det*nCL2[q];
    nCT[4*q+3] =  e1x/det*nCL2[q];
    }

    else
    nCT[4*q] = nCT[4*q+1] = nCT[4*q+2] = nCT[4*q+3] = 0.0;
}

void nhflow_momentum_SI::nh_curv_vel(lexer *p)
{
    // the velocities of [u . grad z] and of the bed condition along the index directions, times the centre distances
    if(!tl2)
    {
    nUZx.assign(lN,0.0);
    nUZy.assign(lN,0.0);
    }

    SLICELOOP4
    if(p->wet[IJ]==1)
    {
    const size_t q = (size_t)((i-p->imin)*p->jmax + (j-p->jmin));

        for(k=0;k<p->knoz;++k)
        {
        // (tl2: the velocities of t^{n+1}, the predictor plus the acceleration of the step before, as nh_prepare)
        const double u = tl2 ? nUZx[PIJK] : UHD[PIJK], v = tl2 ? nUZy[PIJK] : VHD[PIJK];
        nUZx[PIJK] = nCT[4*q]*u + nCT[4*q+1]*v;
        nUZy[PIJK] = nCT[4*q+2]*u + nCT[4*q+3]*v;
        }
    }
}

void nhflow_momentum_SI::nh_cell_faces_curv(lexer *p, fdm_nhf *d, double *cx, double *cy)
{
    // the phi face terms of cell (i,j) as a Cartesian vector per layer: per index direction the mean unit normal and
    // the mean face term of the open faces, the vector from the two (with open faces in one direction only its
    // component along that normal), as the gravity term of cell_update on a curvilinear grid
    const int K = p->knoz;
    auto wn = [&](int a, int b) { return etaI(a,b) + d->depth(a,b) > p->A544; };

    const bool oe = Hx(i,j)>0.0 && wn(i+1,j), ow = Hx(i-1,j)>0.0 && wn(i-1,j);
    const bool on = Hy(i,j)>0.0 && wn(i,j+1), os = Hy(i,j-1)>0.0 && wn(i,j-1);

    double ax=0.0, ay=0.0, bx=0.0, by=0.0, nx, ny;
    int na=0, nb=0;

    if(oe) { fnorm(p,1,i,j,nx,ny);   ax += nx; ay += ny; ++na; }
    if(ow) { fnorm(p,1,i-1,j,nx,ny); ax += nx; ay += ny; ++na; }
    if(on) { fnorm(p,2,i,j,nx,ny);   bx += nx; by += ny; ++nb; }
    if(os) { fnorm(p,2,i,j-1,nx,ny); bx += nx; by += ny; ++nb; }

    const double la = sqrt(ax*ax + ay*ay), lb = sqrt(bx*bx + by*by);

    // tl2: a wall side counts in the mean of the face terms of its direction, with zero (as nh_cell_faces)

    if(tl2 && na>0)
    na += int(wall(p,i,j,1) + wall(p,i,j,2) + 0.5);

    if(tl2 && nb>0)
    nb += int(wall(p,i,j,3) + wall(p,i,j,4) + 0.5);

    if(la>0.0) { ax /= la; ay /= la; }
    if(lb>0.0) { bx /= lb; by /= lb; }

    const double det = ax*by - ay*bx;

    for(k=0;k<K;++k)
    {
    const double ag = na>0 ? ((oe ? NQx[PIJK] : 0.0) + (ow ? NQx[PIm1JK] : 0.0))/double(na) : 0.0;
    const double bg = nb>0 ? ((on ? NQy[PIJK] : 0.0) + (os ? NQy[PIJm1K] : 0.0))/double(nb) : 0.0;

    cx[k] = cy[k] = 0.0;

        if(na>0 && nb>0 && fabs(det)>1.0e-3)
        {
        cx[k] = (ag*by - ay*bg)/det;
        cy[k] = (ax*bg - ag*bx)/det;
        }

        else
        if(na>0)
        {
        cx[k] = ag*ax;
        cy[k] = ag*ay;
        }

        else
        if(nb>0)
        {
        cx[k] = bg*bx;
        cy[k] = bg*by;
        }
    }

    // turned uniform grid: the parts along the two face normals with (-1, 9, 9, -1)/16 of four faces
    if(ho && hoT==1)
    {
    const bool h4x = ho_cell(p,d,1,i,j), h4y = ho_cell(p,d,2,i,j);

        if(h4x || h4y)
        for(k=0;k<K;++k)
        {
        double a1 = cx[k]*hon1x + cy[k]*hon1y, a2 = cx[k]*hon2x + cy[k]*hon2y;

        if(h4x)
        a1 = (9.0*(NQx[PIm1JK] + NQx[PIJK]) - NQx[PIm2JK] - NQx[PIp1JK])/16.0;

        if(h4y)
        a2 = (9.0*(NQy[PIJm1K] + NQy[PIJK]) - NQy[PIJm2K] - NQy[PIJp1K])/16.0;

        cx[k] = a1*hon1x + a2*hon2x;
        cy[k] = a1*hon1y + a2*hon2y;
        }
    }
}

void nhflow_momentum_SI::nh_curv_phi_corr(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // A 539 1: the non-orthogonal part of the phi face terms (NQcx, NQcy; added to the face terms of the step by
    // ho_phi_add, as the explicit part of the fourth-order operators). The face term of layer k is
    //     NQ_k = dt (H_f/h_f)/dsig_k ( -d_n(dsig_k h phi_k^c) + phi_{k+1} d_n z_{k+1} - phi_k d_n z_k ),
    // phi_k^c the mean of phi on the two interfaces of the layer, phi_k, phi_{k+1} the face means; nh_faces takes the
    // normal derivatives d_n from the two cells of the face over their normal distance. Here the rest of them on a
    // non-orthogonal grid, k . grad(f)/L with the mean of the two cell gradients (fcorr, as for eta), explicit in the
    // solve: with phi extrapolated from the two steps before in the first solve of the step, with phi of that solve
    // in a second one (step_single, nhco_it 1). One solve with the extrapolated phi went unstable on a grid sheared
    // by 22 deg (sloshing kh 1.6: the pressure grew by a factor 2.6 per period from 8 periods on), with phi of the step
    // before the wave grew by 2 % per period once the eta part was centred in time; with the second solve it loses
    // 0.06 % per period
    const double dt = p->dt;
    const int K = p->knoz;
    const int N = lN;

    if(NQcx==nullptr)
    {
    p->Darray_lazy(NQcx,N);
    p->Darray_lazy(NQcy,N);
    }

    hoPS.resize(N);

    const double r = ho_phi ? p->dt/MAX(p->dt_old,1.0e-20) : 0.0;

    // a further solve of the step (nhco_it > 0): phi of the solve before
    if(nhco_it>0)
    for(int q=0; q<N; ++q)
    {
    hoPS[q] = PHI[q];
    NQcx[q] = 0.0;
    NQcy[q] = 0.0;
    }

    else
    for(int q=0; q<N; ++q)
    {
    hoPS[q] = PHI[q] + r*(PHI[q] - PHIo[q]);
    PHIo[q] = PHI[q];
    NQcx[q] = 0.0;
    NQcy[q] = 0.0;
    }

    ho_phi = true;

    auto in = [&](int a, int c) { return p->flagslice4[(a-p->imin)*p->jmax + (c-p->jmin)]>0; };
    auto hh = [&](int a, int c) { return nWLh[(a-p->imin)*p->jmax + (c-p->jmin)]; };
    auto ps = [&](int a, int c, int kk) { return kk<K ? hoPS[PIX(a,c,kk)] : 0.0; };

    const size_t NS = (size_t)p->imax*p->jmax;
    auto sq = [&](int a, int c) { return (size_t)((a-p->imin)*p->jmax + (c-p->jmin)); };

    // the non-orthogonal part of the face derivative of a cell field (in tS, all cells of the grid with the halo)
    std::vector<double> cbx(NS,0.0), cby(NS,0.0), chx(NS,0.0), chy(NS,0.0), cax(NS,0.0), cay(NS,0.0);

    auto corr = [&](std::vector<double> &fx, std::vector<double> &fy)
    {
        exchange4(p,pgc,tS);

        SLICELOOP1
        fx[sq(i,j)] = (Hx(i,j)>0.0 && in(i,j) && in(i+1,j)) ? fcorr(p,tS,1,i,j) : 0.0;

        SLICELOOP2
        fy[sq(i,j)] = (Hy(i,j)>0.0 && in(i,j) && in(i,j+1)) ? fcorr(p,tS,2,i,j) : 0.0;
    };

    auto fill = [&](const std::function<double(int,int)> &f)
    {
        for(int a=p->imin; a<p->imin+p->imax; ++a)
        for(int c=p->jmin; c<p->jmin+p->jmax; ++c)
        tS(a,c) = f(a,c);
    };

    // the bed and the layer thickness of the nh terms (z_m = bed + sig_m h)
    fill([&](int a, int c) { return d->bed(a,c); });
    corr(cbx,cby);
    fill([&](int a, int c) { return hh(a,c); });
    corr(chx,chy);

    for(int kk=0; kk<K; ++kk)
    {
    const double s = p->DZN[kk+marge];
    const double sl = p->ZN[kk+marge], su = p->ZN[kk+1+marge];

    fill([&](int a, int c) { return s*hh(a,c)*0.5*(ps(a,c,kk) + ps(a,c,kk+1)); });
    corr(cax,cay);

    k = kk;

        SLICELOOP1
        if(Hx(i,j)>0.0 && in(i,j) && in(i+1,j))
        {
        const double hf = 0.5*(hh(i,j) + hh(i+1,j));
        const double fm = 0.5*(ps(i,j,kk) + ps(i+1,j,kk)), um = 0.5*(ps(i,j,kk+1) + ps(i+1,j,kk+1));
        const size_t q = sq(i,j);

        if(hf>1.0e-20)
        NQcx[PIJK] = dt*(Hx(i,j)/hf)/s*(-cax[q] + um*(cbx[q] + su*chx[q]) - fm*(cbx[q] + sl*chx[q]));
        }

        SLICELOOP2
        if(Hy(i,j)>0.0 && in(i,j) && in(i,j+1))
        {
        const double hf = 0.5*(hh(i,j) + hh(i,j+1));
        const double fm = 0.5*(ps(i,j,kk) + ps(i,j+1,kk)), um = 0.5*(ps(i,j,kk+1) + ps(i,j+1,kk+1));
        const size_t q = sq(i,j);

        if(hf>1.0e-20)
        NQcy[PIJK] = dt*(Hy(i,j)/hf)/s*(-cay[q] + um*(cby[q] + su*chy[q]) - fm*(cby[q] + sl*chy[q]));
        }
    }

    // open boundary faces: no pressure term
    for(const bc_face &f : bcf)
    {
    i=f.fi;
    j=f.fj;
    double *NQ = (f.cs==1 || f.cs==4) ? NQcx : NQcy;
        for(k=0;k<K;++k)
        NQ[PIJK] = 0.0;
    }

    pexchF(p,pgc,NQcx);
    pexchF(p,pgc,NQcy);
}


// ---------------------------------------------------------------------------------------------------------------------
// The wall sides of a cell in the mean of its face terms on a curvilinear grid with the second-order time levels
// (cell_update; the Cartesian branch there): a wall side of an index direction with an open face counts as a face,
// its term is the one that keeps the normal velocity at the wall at zero, g H d_n eta = the column mean of the sources
// of the predictor along the mean normal of the open faces of that direction (zero in free waves). Without the wall
// sides a wall cell took the full gradient of its one open face, as the Cartesian step did before them.
// (ax, ay), (bx, by): the sums of the normals of the open faces of the two index directions, na, nb their numbers,
// ag, bg the sums of their face terms (per dt g). wa, wb: the number of wall sides counted
// ---------------------------------------------------------------------------------------------------------------------
void nhflow_momentum_SI::curv_wall_sides(lexer *p, double ax, double ay, double bx, double by, int na, int nb,
                                         int &wa, int &wb, double &ag, double &bg)
{
    wa = na>0 ? int(wall(p,i,j,1) + wall(p,i,j,2) + 0.5) : 0;
    wb = nb>0 ? int(wall(p,i,j,3) + wall(p,i,j,4) + 0.5) : 0;

    if(wa==0 && wb==0)
    return;

    const double dt = p->dt;
    double sx=0.0, sy=0.0;

    for(k=0;k<p->knoz;++k)
    {
    sx += p->DZN[KP]*dt*(FSx[PIJK] - (Bx!=nullptr ? WL0(i,j)*Bx[IJK] : 0.0));
    sy += p->DZN[KP]*dt*(FSy[PIJK] - (By!=nullptr ? WL0(i,j)*By[IJK] : 0.0));
    }

    if(wa>0)
    {
    const double l = sqrt(ax*ax + ay*ay);
    ag += double(wa)*(l>0.0 ? (sx*ax + sy*ay)/l : 0.0)/(dt*g);
    }

    if(wb>0)
    {
    const double l = sqrt(bx*bx + by*by);
    bg += double(wb)*(l>0.0 ? (sx*bx + sy*by)/l : 0.0)/(dt*g);
    }
}
