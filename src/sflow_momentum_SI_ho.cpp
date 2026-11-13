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

#include"sflow_momentum_SI.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"hydro_sflow.h"
#include<cmath>

// ---------------------------------------------------------------------------------------------------------------------
// Fourth-order horizontal operators of the semi-implicit SFLOW, on one uniform Cartesian grid (as NHFLOW v15 / 0004)
//
// The second-order step carries linear waves too slowly, with an error of (k dx)^2/6 in the wave speed of a long wave:
// (k dx)^2/8 from the two means between the cells and the faces (the cell predictor to the faces, the face terms of
// gravity and of the non-hydrostatic pressure back to the cells) and (k dx)^2/24 from the compact gradient and
// divergence. Here
//  - the cell predictor goes to the faces with (-1, 9, 9, -1)/16 of four cells (ho_pred),
//  - the face terms go to the cells with (-1, 9, 9, -1)/16 of four faces (cell_update, nh_cell),
//  - the face gradient of eta and the pressure terms of the faces get the part -delta^2/12 of the compact ones
//    (ho_gradient, ho_phi_corr): with the compact divergence of the faces, which stays as it is (the volume is exact),
//    the Laplacian is fourth order. This part is explicit, so the eta system and the coupled eta-P system keep their
//    compact matrices: for eta^{n+1} it takes eta extrapolated from t^n and t^{n-1}, for P the pressure of the two
//    steps before extrapolated to this one.
// A face or a cell whose stencil leaves the wet deep fluid cells (walls, open boundaries, the shoreline, dry cells)
// keeps the second-order terms. Second order stays on curvilinear, multi-block and refined grids, with REEF3D::HYDRO
// (rain and point sources) and with the moving bodies X 10 2/3. The A 220 3 dispersion term h^3 lap(eta) and the
// advection of the predictor are unchanged.
// ---------------------------------------------------------------------------------------------------------------------

void sflow_momentum_SI::ho_begin(lexer *p, fdm2D *b, ghostcell *pgc)
{
    // the grid: one Cartesian grid, uniform spacing of a direction on all ranks
    if(!ho_grid)
    {
    ho_grid = true;

    double ux=0.0, uy=0.0;

    for(int q=0; q<p->knox; ++q)
    ux = MAX(ux,fabs(p->DXN[q+marge] - p->DXN[marge])/p->DXN[marge]);

    for(int q=0; q<p->knoy; ++q)
    uy = MAX(uy,fabs(p->DYN[q+marge] - p->DYN[marge])/p->DYN[marge]);

    const bool one = p->geo_curv==0 && prun==nullptr && !patch && !amrgrid && !hydro_sflow::wetting(p)
                     && p->X10!=2 && p->X10!=3;

    ho_ux = one && pgc->globalmax(ux)<1.0e-8;
    ho_uy = one && pgc->globalmax(uy)<1.0e-8 && p->j_dir==1;
    }

    ho = ho_ux || ho_uy;

    if(!ho)
    return;

    // wet deep fluid cells (WL above 10 A 244, not next to the shoreline), with those of the neighbouring ranks (the
    // ghost cells of walls and open boundaries: 0), where the flow is smooth: a Froude number below 1/2 and a change of
    // the depth to the next cells below 10 % of the depth (a wave of a/h 0.1 at 16 cells per wave length: 4 %). Bores,
    // steep fronts and supercritical flow keep the second-order operators: there is no limiter, and across the bore of
    // the Stoker dam break the fourth-order operators raised the error from 0.15 to 0.8 % (0.21 % with the Froude
    // number alone: the foot of the bore). And a Courant number of the shortest wave below 1 (below)
    auto jump = [&](int a, int c) { return p->flagslice4[(a-p->imin)*p->jmax + (c-p->jmin)]>0 ? fabs(WL0(a,c) - WL0(i,j)) : 0.0; };

    SLICELOOP4
    {
    const double h = WL0(i,j);
    bool f = p->flagslice4[IJ]>0 && p->wet[IJ]==1 && p->deep[IJ]==1 && h>10.0*p->A244;

    f = f && (UH0(i,j)*UH0(i,j) + VH0(i,j)*VH0(i,j)) < 0.25*g*h*h*h;
    f = f && MAX(jump(i-1,j),jump(i+1,j)) < 0.1*h;

    // the explicit part is stable up to a Courant number of about 1 of the shortest wave of the grid (2 dx), with the
    // phase speed of the model: sqrt(g h) hydrostatic, sqrt(g h (1 + B (kh)^2)/(1 + (B + 1/3)(kh)^2)) with the
    // non-hydrostatic pressure, kh = pi h/dx (B of A 220 3, 0 with A 220 2). Hydrostatic standing waves at 16 cells per
    // wave length: stable at c dt/dx 1.0, growing from 1.45 on (sflow_2d_hydro_soil_si, c dt/dx 2.8: 0.1 m/s in a
    // basin at rest); larger steps keep the second-order operators there
    if(f)
    {
    const double dxm = p->j_dir==1 ? MIN(p->DXN[IP],p->DYN[JP]) : p->DXN[IP];
    const double kh = PI*h/dxm;
    const double D = nh ? (1.0 + Bdisp*kh*kh)/(1.0 + (Bdisp + 1.0/3.0)*kh*kh) : 1.0;
    f = sqrt(g*h*D)*p->dt <= dxm;
    }

    if(p->j_dir==1)
    f = f && MAX(jump(i,j-1),jump(i,j+1)) < 0.1*h;

    hoC(i,j) = f ? 1.0 : 0.0;
    }

    exchange4(p,pgc,hoC);

    // eta^{n+1} extrapolated from t^n and t^{n-1} (the explicit part of the fourth-order gradient of eta^{n+1})
    const double r = ho_eo ? p->dt/MAX(p->dt_old,1.0e-20) : 0.0;

    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
    hoEx(ii,jj) = (1.0+r)*eta0(ii,jj) - r*hoEo(ii,jj);
    hoEo(ii,jj) = eta0(ii,jj);
    }

    ho_eo = true;
}

// face (ii,jj) of direction dir (1: x, 2: y): the four cells of its stencil are wet deep fluid cells
bool sflow_momentum_SI::ho_face(lexer *p, int dir, int ii, int jj)
{
    if(!ho)
    return false;

    if(dir==1)
    return ho_ux && hoC(ii-1,jj)>0.5 && hoC(ii,jj)>0.5 && hoC(ii+1,jj)>0.5 && hoC(ii+2,jj)>0.5;

    return ho_uy && hoC(ii,jj-1)>0.5 && hoC(ii,jj)>0.5 && hoC(ii,jj+1)>0.5 && hoC(ii,jj+2)>0.5;
}

// cell (ii,jj): its four faces of direction dir lie between wet deep fluid cells that are wet at the end of the step
bool sflow_momentum_SI::ho_cell(lexer *p, fdm2D *b, int dir, int ii, int jj)
{
    if(!ho || (dir==1 && !ho_ux) || (dir==2 && !ho_uy))
    return false;

    const int di = dir==1 ? 1 : 0, dj = dir==2 ? 1 : 0;

    for(int m=-2; m<=2; ++m)
    {
    const int a = ii+m*di, c = jj+m*dj;

    if(hoC(a,c)<0.5 || etaI(a,c) + b->depth(a,c) <= p->A244)
    return false;
    }

    return true;
}

// the face gradient of the face (i,j): g0 that of eta^n, c0 the explicit part of that of eta^{n+1}; the part beyond
// (f_R - f_L)/dx is -delta^2 (f_R - f_L)/(12 dx)
void sflow_momentum_SI::ho_gradient(lexer *p, fdm2D *b, int dir, double &g0, double &c0)
{
    c0 = 0.0;

    if(!ho_face(p,dir,i,j))
    return;

    if(dir==1)
    {
    const double c = 1.0/(12.0*p->DXP[IP]);
    g0 -= c*(eta0(i+2,j) - 3.0*eta0(i+1,j) + 3.0*eta0(i,j) - eta0(i-1,j));
    c0  = -c*(hoEx(i+2,j) - 3.0*hoEx(i+1,j) + 3.0*hoEx(i,j) - hoEx(i-1,j));
    }

    else
    {
    const double c = 1.0/(12.0*p->DYP[JP]);
    g0 -= c*(eta0(i,j+2) - 3.0*eta0(i,j+1) + 3.0*eta0(i,j) - eta0(i,j-1));
    c0  = -c*(hoEx(i,j+2) - 3.0*hoEx(i,j+1) + 3.0*hoEx(i,j) - hoEx(i,j-1));
    }
}

// the face predictor of the face (i,j) from the cell array Q (the momentum normal to the face): (-1, 9, 9, -1)/16 of
// four cells instead of the mean of two
double sflow_momentum_SI::ho_pred(lexer *p, int dir, slice &Q)
{
    if(dir==1)
    return (9.0*(Q(i,j) + Q(i+1,j)) - Q(i-1,j) - Q(i+2,j))/16.0;

    return (9.0*(Q(i,j) + Q(i,j+1)) - Q(i,j-1) - Q(i,j+2))/16.0;
}

// the non-hydrostatic pressure terms of the faces (nh_assemble, after the face coefficients): q_NH = -R P_R + L P_L
// with -delta^2 P/12 along the face normal in place of P as the explicit part (the constant c0 of the face), P
// extrapolated to this step from the two steps before; the four cells of the stencil non-hydrostatic
void sflow_momentum_SI::ho_phi_corr(lexer *p, fdm2D *b, ghostcell *pgc)
{
    if(!ho)
    return;

    const double r = nh_x0 ? p->dt/MAX(p->dt_old,1.0e-20) : 0.0;

    // P extrapolated (nPo: P^{n-1}, nh_solve), with the cells of the neighbouring ranks (nZ: scratch of nh_solve)
    SLICELOOP4
    nZ(i,j) = (1.0+r)*b->press(i,j)/p->W1 - r*nPo(i,j);

    exchange4(p,pgc,nZ);

    auto pe = [&](int a, int c) { return nZ(a,c); };
    auto act = [&](int a, int c) { return nA(a,c)>0.5; };

    SLICELOOP1
    if(Hx(i,j)>0.0 && ho_face(p,1,i,j) && act(i-1,j) && act(i,j) && act(i+1,j) && act(i+2,j))
    {
    const double dL = pe(i+1,j) - 2.0*pe(i,j) + pe(i-1,j);
    const double dR = pe(i+2,j) - 2.0*pe(i+1,j) + pe(i,j);
    nx0(i,j) += (nxR(i,j)*dR - nxL(i,j)*dL)/12.0;
    }

    if(p->j_dir==1)
    SLICELOOP2
    if(Hy(i,j)>0.0 && ho_face(p,2,i,j) && act(i,j-1) && act(i,j) && act(i,j+1) && act(i,j+2))
    {
    const double dL = pe(i,j+1) - 2.0*pe(i,j) + pe(i,j-1);
    const double dR = pe(i,j+2) - 2.0*pe(i,j+1) + pe(i,j);
    ny0(i,j) += (nyR(i,j)*dR - nyL(i,j)*dL)/12.0;
    }

    // The A 220 3 dispersion term g B grad(h^3 lap(eta)) keeps its compact operators: its fourth-order part as an
    // explicit term would remove the remaining phase error of A 220 3, but it is stiff (unstable at practical steps)
}
