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
// Fourth-order horizontal operators of the semi-implicit step (A 508 1), on one uniform Cartesian grid
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
//    their compact stencil: for eta^{n+1} it takes eta of the estimate pass (A 505 1; eta^n in the estimate pass and
//    with A 505 0), for phi the pressure of the two steps before extrapolated to this one.
// A face or a cell whose stencil leaves the wet deep fluid cells (walls, open boundaries, the shoreline) keeps the
// second-order terms. The advection of the predictor and the nonlinear face depths are unchanged.
// ---------------------------------------------------------------------------------------------------------------------

void nhflow_momentum_SI::ho_begin(lexer *p, fdm_nhf *d, ghostcell *pgc)
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

    // wet deep fluid cells, with those of the neighbouring ranks (the ghost cells of walls and open boundaries: 0)
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
    hoC(ii,jj) = (p->flagslice4[q]>0 && p->wet[q]==1 && p->deep[q]==1) ? 1.0 : 0.0;
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
    const double c = 1.0/(12.0*p->DXP[IP]);
    g0 -= c*(d->eta(i+2,j) - 3.0*d->eta(i+1,j) + 3.0*d->eta(i,j) - d->eta(i-1,j));
    c0 -= c*(e1(i+2,j) - 3.0*e1(i+1,j) + 3.0*e1(i,j) - e1(i-1,j));
    }

    else
    {
    const double c = 1.0/(12.0*p->DYP[JP]);
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

    if(dir==1)
    for(k=0;k<p->knoz;++k)
    w[k] = (9.0*(c1*UHc[IJK] + c2*UHc[Ip1JK]) - c0*UHc[Im1JK] - c3*UHc[Ip2JK])/16.0;

    else
    for(k=0;k<p->knoz;++k)
    w[k] = (9.0*(c1*UHc[IJK] + c2*UHc[IJp1K]) - c0*UHc[IJm1K] - c3*UHc[IJp2K])/16.0;
}

// the gravity term of the cell (i,j) (cell_update): dt g times (-1, 9, 9, -1)/16 of H_f times the face gradient
// (theta weighted, with the explicit part) of its four faces of a direction
void nhflow_momentum_SI::ho_cell_gravity(lexer *p, fdm_nhf *d, double &gu, double &gv)
{
    const double dt = p->dt;

    if(ho_cell(p,d,1,i,j))
    {
    auto aF = [&](int f)
    {
        return Hx(f,j)*(theta*((etaI(f+1,j) - etaI(f,j))/p->DXP[f+marge] + cx0(f,j)) + (1.0-theta)*gx0(f,j));
    };

    gu = dt*g*(9.0*(aF(i-1) + aF(i)) - aF(i-2) - aF(i+1))/16.0;
    }

    if(p->j_dir==1 && ho_cell(p,d,2,i,j))
    {
    auto aF = [&](int f)
    {
        return Hy(i,f)*(theta*((etaI(i,f+1) - etaI(i,f))/p->DYP[f+marge] + cy0(i,f)) + (1.0-theta)*gy0(i,f));
    };

    gv = dt*g*(9.0*(aF(j-1) + aF(j)) - aF(j-2) - aF(j+1))/16.0;
    }
}
