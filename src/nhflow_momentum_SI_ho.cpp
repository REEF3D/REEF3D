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
// Fourth-order horizontal operators of the semi-implicit step, on one uniform Cartesian grid and on one curvilinear
// grid (G 50 2) that is a uniform grid of rectangles turned by any angle (ho_turned). On the turned grid the stencils
// run along the two index directions with the normal distances of the centres: eta and phi as on the Cartesian grid,
// the momentum of the cells through its components along the two face normals (ho_predictor, ho_cell_gravity,
// nh_cell_faces_curv). A flume of 15 cells per wavelength turned by 30 deg then gives the waves of its Cartesian twin
// to 0.03 mm of 10.5 mm amplitude after 10 wavelengths (second order: 6.7 mm off).
//
// The second-order step carries linear waves too slowly, with an error of (k dx)^2/6 in the wave speed of a long wave:
// (k dx)^2/8 from the two means between the cells and the faces (the cell predictor to the faces, the face terms of
// gravity and of the non-hydrostatic pressure back to the cells) and (k dx)^2/24 from the compact gradient and
// divergence. Here
//  - the cell predictor goes to the faces with (-1, 9, 9, -1)/16 of four cells (ho_predictor),
//  - the face terms go to the cells with (-1, 9, 9, -1)/16 of four faces (ho_cell_gravity, nh_cell),
//  - the face gradient of eta and the face term of phi get the part -delta^2/12 of the compact one (ho_gradient,
//    ho_phi_corr): with the compact divergence of the faces, which stays as it is (the volume is exact), the Laplacian
//    is fourth order. This part is explicit, so the matrices of the eta system and of the coupled eta-phi system keep
//    their compact stencil: for eta^{n+1} it takes eta of the estimate pass (eta^n in the estimate pass), for phi the
//    pressure of the two steps before extrapolated to this one.
// A face or a cell whose stencil leaves the wet deep fluid cells (walls, open boundaries, the shoreline) keeps the
// second-order terms. The advection of the predictor and the nonlinear face depths are unchanged.
// ---------------------------------------------------------------------------------------------------------------------

void nhflow_momentum_SI::ho_begin(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // curvilinear grid (G 50 2): only a uniform grid of rectangles, turned (ho_turned: checked once)
    if(p->geo_curv==1)
    {
        if(hoT<0)
        ho_turned(p,pgc);

    ho = (hoT==1);

        if(!ho)
        return;

    ho_ux = ho_uy = true;
    }

    else
    {
    ho = true;

    // uniform spacing of the direction, on all ranks
    double ux=0.0, uy=0.0;

    for(int q=0; q<p->knox; ++q)
    ux = MAX(ux,fabs(p->DXN[q+marge] - p->DXN[marge])/p->DXN[marge]);

    for(int q=0; q<p->knoy; ++q)
    uy = MAX(uy,fabs(p->DYN[q+marge] - p->DYN[marge])/p->DYN[marge]);

    ho_ux = pgc->globalmax(ux)<1.0e-8;
    ho_uy = pgc->globalmax(uy)<1.0e-8 && p->j_dir==1;
    }

    // wet deep fluid cells, with those of the neighbouring ranks (the ghost cells of walls and open boundaries: 0)
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
    hoC(ii,jj) = (p->flagslice4[q]>0 && p->wet[q]==1 && p->deep[q]==1 && !df_col(p,ii,jj)) ? 1.0 : 0.0;

    // VRANS: no fourth-order stencil over a column with a porous cell
    if(vr && vrPc[q])
    hoC(ii,jj) = 0.0;
    }

    // and the cells of subcritical flow, Froude number of the layers below 0.8 (t^n): at the critical point of a dam break
    // (the tail of the rarefaction, u = c) the cubic overshoot of the (-1, 9, 9, -1)/16 means grew a two-cell wave at
    // gravity Courant numbers near 1 (Stoker dam break, theta 1, dx 0.05 m: EMERGENCY STOP at 3.4 s with dt 0.02 s, at
    // 1.5 s with dt 0.03 s). Waves (the Froude number of the orbital velocity well below 0.5) and rivers keep them. A
    // limiter of the means instead (v24: the harmonic mean of the two curvatures at kinks) stopped the time convergence
    // of the standing waves near 2e-7 (its switches are not smooth in time)
    const double frl = 0.8;

    LOOP
    if(hoC(i,j)>0.5 && WL0(i,j)>p->A544)
    {
    const double v = p->j_dir==1 ? d->V[IJK] : 0.0;

        if(d->U[IJK]*d->U[IJK] + v*v > frl*frl*g*WL0(i,j))
        hoC(i,j) = 0.0;
    }

    exchange4(p,pgc,hoC);
}

// eta of the estimate pass, for the explicit part of the gradient of eta^{n+1} in the second pass
void nhflow_momentum_SI::ho_estimate(lexer *p, ghostcell *pgc)
{
    exchange4(p,pgc,etaI);

    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    hoE(ii,jj) = etaI(ii,jj);
}

// face (ii,jj) of direction dir (1: x, 2: y): the four cells of its stencil are wet deep fluid cells
bool nhflow_momentum_SI::ho_face(lexer *p, int dir, int ii, int jj)
{
    if(dir==1)
    return ho_ux && hoC(ii-1,jj)>0.5 && hoC(ii,jj)>0.5 && hoC(ii+1,jj)>0.5 && hoC(ii+2,jj)>0.5;

    return ho_uy && hoC(ii,jj-1)>0.5 && hoC(ii,jj)>0.5 && hoC(ii,jj+1)>0.5 && hoC(ii,jj+2)>0.5;
}

// cell (ii,jj): its four faces of direction dir lie between wet deep fluid cells that are wet at the end of the step
bool nhflow_momentum_SI::ho_cell(lexer *p, fdm_nhf *d, int dir, int ii, int jj)
{
    if((dir==1 && !ho_ux) || (dir==2 && !ho_uy))
    return false;

    const int di = dir==1 ? 1 : 0, dj = dir==2 ? 1 : 0;

    for(int m=-2; m<=2; ++m)
    {
    const int a = ii+m*di, b = jj+m*dj;

    if(hoC(a,b)<0.5 || etaI(a,b) + d->depth(a,b) <= p->A544)
    return false;
    }

    return true;
}

// the face gradient of the face (i,j): g0 that of eta^n, c0 the explicit part of that of eta^{n+1}; the part beyond
// (f_R - f_L)/dx is -delta^2 (f_R - f_L)/(12 dx)
void nhflow_momentum_SI::ho_gradient(lexer *p, fdm_nhf *d, int dir, double &g0, double &c0)
{
    if(!ho_face(p,dir,i,j))
    return;

    slice &e1 = si_pass>0 ? (slice&)hoE : (slice&)d->eta;

    if(dir==1)
    {
    const double c = 1.0/(12.0*(hoT==1 ? hod1 : p->DXP[IP]));
    g0 -= c*(d->eta(i+2,j) - 3.0*d->eta(i+1,j) + 3.0*d->eta(i,j) - d->eta(i-1,j));
    c0 -= c*(e1(i+2,j) - 3.0*e1(i+1,j) + 3.0*e1(i,j) - e1(i-1,j));
    }

    else
    {
    const double c = 1.0/(12.0*(hoT==1 ? hod2 : p->DYP[JP]));
    g0 -= c*(d->eta(i,j+2) - 3.0*d->eta(i,j+1) + 3.0*d->eta(i,j) - d->eta(i,j-1));
    c0 -= c*(e1(i,j+2) - 3.0*e1(i,j+1) + 3.0*e1(i,j) - e1(i,j-1));
    }
}

// the column w of the face (i,j) of the cell array UHc (the momentum normal to the face): (-1, 9, 9, -1)/16 of the cell
// transports, or of the cell velocities times the face depth of the mass flux (hMx, hMy >= 0, nh_begin), as the
// second-order mean of si_assemble
void nhflow_momentum_SI::ho_predictor(lexer *p, int dir, const double *UHc, double *w)
{
    if(!ho_face(p,dir,i,j))
    return;

    const int di = dir==1 ? 1 : 0, dj = dir==2 ? 1 : 0;
    const double hm = dir==1 ? hMx[IJ] : hMy[IJ];
    const double c0 = hm<0.0 ? 1.0 : hm/WL0(i-di,j-dj);
    const double c1 = hm<0.0 ? 1.0 : hm/WL0(i,j);
    const double c2 = hm<0.0 ? 1.0 : hm/WL0(i+di,j+dj);
    const double c3 = hm<0.0 ? 1.0 : hm/WL0(i+2*di,j+2*dj);

    // turned grid: the momentum normal to the face from the two Cartesian components of the cells
    if(hoT==1)
    {
    const bool old = (UHc==UH0 || UHc==VH0);
    const double *A = old ? UH0 : UHD, *B = old ? VH0 : VHD;
    const double nx = dir==1 ? hon1x : hon2x, ny = dir==1 ? hon1y : hon2y;
    const int s = dir==1 ? lsi : lsj;

        for(k=0;k<p->knoz;++k)
        {
        const int q = PIJK;
        w[k] = (9.0*(c1*(nx*A[q] + ny*B[q]) + c2*(nx*A[q+s] + ny*B[q+s]))
                - c0*(nx*A[q-s] + ny*B[q-s]) - c3*(nx*A[q+2*s] + ny*B[q+2*s]))/16.0;
        }

    return;
    }

    if(dir==1)
    for(k=0;k<p->knoz;++k)
    w[k] = (9.0*(c1*UHc[PIJK] + c2*UHc[PIp1JK]) - c0*UHc[PIm1JK] - c3*UHc[PIp2JK])/16.0;

    else
    for(k=0;k<p->knoz;++k)
    w[k] = (9.0*(c1*UHc[PIJK] + c2*UHc[PIJp1K]) - c0*UHc[PIJm1K] - c3*UHc[PIJp2K])/16.0;
}

// the gravity term of the cell (i,j) (cell_update): dt g times (-1, 9, 9, -1)/16 of H_f times the face gradient
// (theta weighted, with the explicit part) of its four faces of a direction
void nhflow_momentum_SI::ho_cell_gravity(lexer *p, fdm_nhf *d, double &gu, double &gv)
{
    const double dt = p->dt;

    // turned grid: the parts of the cell vector along the two face normals (orthogonal), each from its four faces
    if(hoT==1)
    {
    double a1 = gu*hon1x + gv*hon1y, a2 = gu*hon2x + gv*hon2y;

        if(ho_cell(p,d,1,i,j))
        {
        auto aF = [&](int f)
        {
            return Hx(f,j)*(theta*((etaI(f+1,j) - etaI(f,j))/hod1 + cx0(f,j)) + (1.0-theta)*gx0(f,j));
        };

        a1 = dt*g*(9.0*(aF(i-1) + aF(i)) - aF(i-2) - aF(i+1))/16.0;
        }

        if(ho_cell(p,d,2,i,j))
        {
        auto aF = [&](int f)
        {
            return Hy(i,f)*(theta*((etaI(i,f+1) - etaI(i,f))/hod2 + cy0(i,f)) + (1.0-theta)*gy0(i,f));
        };

        a2 = dt*g*(9.0*(aF(j-1) + aF(j)) - aF(j-2) - aF(j+1))/16.0;
        }

    gu = a1*hon1x + a2*hon2x;
    gv = a1*hon1y + a2*hon2y;
    return;
    }

    if(ho_cell(p,d,1,i,j))
    {
    auto aF = [&](int f)
    {
        return Hx(f,j)*(theta*((etaI(f+1,j) - etaI(f,j))/p->DXP[f+marge] + cx0(f,j)) + (1.0-theta)*gx0(f,j)
                        + (shipon ? shx(f,j) : 0.0));
    };

    gu = dt*g*(9.0*(aF(i-1) + aF(i)) - aF(i-2) - aF(i+1))/16.0;
    }

    if(p->j_dir==1 && ho_cell(p,d,2,i,j))
    {
    auto aF = [&](int f)
    {
        return Hy(i,f)*(theta*((etaI(i,f+1) - etaI(i,f))/p->DYP[f+marge] + cy0(i,f)) + (1.0-theta)*gy0(i,f)
                        + (shipon ? shy(i,f) : 0.0));
    };

    gv = dt*g*(9.0*(aF(j-1) + aF(j)) - aF(j-2) - aF(j+1))/16.0;
    }
}

// curvilinear grid (G 50 2): hoT 1 if it is a uniform grid of rectangles (turned by any angle): on all ranks the
// same unit normals of the faces of a direction, orthogonal to each other, the same normal distance of the centres
// and the same face lengths. The fourth-order operators then act along the two index directions as on the Cartesian
// grid: scalars (eta, phi) as they are, the momentum of the cells through its components along the two normals
void nhflow_momentum_SI::ho_turned(lexer *p, ghostcell *pgc)
{
    double n1x=0.0, n1y=0.0, n2x=0.0, n2y=0.0, d1=0.0, d2=0.0, l1=0.0, l2=0.0;
    double v[8] = {-1.0e30,-1.0e30,-1.0e30,-1.0e30,-1.0e30,-1.0e30,-1.0e30,-1.0e30};
    double w[8] = {-1.0e30,-1.0e30,-1.0e30,-1.0e30,-1.0e30,-1.0e30,-1.0e30,-1.0e30};

    // the largest and the smallest of every quantity over the interior faces of all ranks
    auto take = [&](int q, double x) { v[q] = MAX(v[q],x); w[q] = MAX(w[q],-x); };
    auto in = [&](int a, int b) { return p->flagslice4[(a-p->imin)*p->jmax + (b-p->jmin)]>0; };

    SLICELOOP4
    {
        if(in(i+1,j))
        {
        fnorm(p,1,i,j,n1x,n1y);
        take(0,n1x); take(1,n1y); take(2,fdn(p,1,i,j)); take(3,flen(p,1,i,j));
        }

        if(in(i,j+1))
        {
        fnorm(p,2,i,j,n2x,n2y);
        take(4,n2x); take(5,n2y); take(6,fdn(p,2,i,j)); take(7,flen(p,2,i,j));
        }
    }

    // (1e-6: node coordinates of the order of 1e6 m with cells of 0.1 m leave 1e-8 of the cell size)
    const double tol = 1.0e-6;
    bool ok = true;

    for(int q=0;q<8;++q)
    {
    v[q] = pgc->globalmax(v[q]);
    w[q] = -pgc->globalmax(w[q]);

    const double sc = (q==2 || q==3 || q==6 || q==7) ? MAX(fabs(v[q]),1.0e-20) : 1.0;

    if(v[q]<-1.0e29 || fabs(v[q]-w[q])>tol*sc)
    ok = false;
    }

    hon1x = v[0]; hon1y = v[1]; hod1 = v[2];
    hon2x = v[4]; hon2y = v[5]; hod2 = v[6];

    // orthogonal normals, rectangles (the face length of a direction is the centre distance of the other one)
    if(ok && (fabs(hon1x*hon2x + hon1y*hon2y)>tol || fabs(v[3]-hod2)>tol*hod2 || fabs(v[7]-hod1)>tol*hod1))
    ok = false;

    hoT = ok ? 1 : 0;

    if(p->mpirank==0 && !patch)
    {
        if(hoT==1)
        cout<<"NHFLOW semi-implicit: the curvilinear grid is a turned uniform grid ("<<hod1<<" x "<<hod2
            <<" m), fourth-order horizontal operators"<<endl;

        else
        cout<<"NHFLOW semi-implicit: curvilinear grid, second-order horizontal operators (fourth order on a uniform "
            <<"grid of rectangles only)"<<endl;
    }
}
