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

#ifndef NHFLOW_CURV_H_
#define NHFLOW_CURV_H_

#include"increment.h"
#include"lexer.h"

class fdm_nhf;
class slice;

using namespace std;

// NHFLOW on a curvilinear horizontal grid (G 50, p->geo_curv==1), the sigma grid in the vertical.
//
// As SFLOW (sflow_curv.h): the state stays cell-centred with Cartesian velocity components (U,V,W)
// and conserved variables (UH,VH,WH = velocity times the water depth WL); the faces of the logical
// grid carry the normal vectors of the geometry layer, n1 = (SX1,SY1)/LEN1 on face i+1/2 and
// n2 = (SX2,SY2)/LEN2 on face j+1/2.  The reconstruction stays in index space, the Riemann problem of
// a face is solved in its normal direction:
//   - the face velocities Us/Un (i-faces) and Ve/Vw (j-faces) hold the normal velocity u.n
//   - momentum flux of UH: UH (u.n) + P n_x, of VH: VH (u.n) + P n_y, of WH: WH (u.n),
//     P = 0.5 g eta^2 + g eta d_face (pre-balanced)
//   - continuity flux: UH n_x + VH n_y
//   - horizontal divergence: sum of the face fluxes times the face lengths over the cell area
//   - bed-slope source: g eta sum(d_face S_face)/A (a lake at rest stays at rest on any grid)
// Horizontal derivatives (sigma metrics, pressure gradient, divergence of the Poisson right-hand
// side) by the chain rule with the cell centres (sflow_curv::grad5); the Laplacians compact with the
// face lengths and the centre distances along the face normals, the non-orthogonal part deferred.
// The Cartesian kernels are not touched: every curvilinear path is a branch on p->geo_curv.

class nhflow_curv : public increment
{
public:
    // unsupported options are refused
    static void check(lexer*);

    // unit normal (towards +i / +j) and length of face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (i,j)
    static inline void normal(lexer *p, int dir, int ii, int jj, double &nx, double &ny, double &len)
    {
        const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
        len = (dir==1) ? p->LEN1[q] : p->LEN2[q];
        if(len>0.0)
        {
        nx = ((dir==1) ? p->SX1[q] : p->SX2[q])/len;
        ny = ((dir==1) ? p->SY1[q] : p->SY2[q])/len;
        }
        else
        {
        nx = (dir==1) ? 1.0 : 0.0;
        ny = (dir==1) ? 0.0 : 1.0;
        }
    }

    // after the reconstruction: Us/Un = u.n1 from (Us,Vs)/(Un,Vn), Ve/Vw = u.n2 from (Ue,Ve)/(Uw,Vw)
    static void normal_velocities(lexer*, fdm_nhf*);

    // horizontal divergence of face fluxes per unit length (i-faces Fx, j-faces Fy) of cell (i,j,k)
    static inline double divh(lexer *p, const double *Fx, const double *Fy, int i, int j, int k)
    {
        return (Fx[IJK]*p->LEN1[IJ] - Fx[Im1JK]*p->LEN1[Im1J]
              + Fy[IJK]*p->LEN2[IJ] - Fy[IJm1K]*p->LEN2[IJm1])/p->AREA2D[IJ];
    }

    // bed-slope source g eta sum(d_face S_face)/A of cell (i,j): comp 1 x, comp 2 y
    static double bedslope(lexer*, fdm_nhf*, int comp, int i, int j);

    // gradient of a cell-centred 3D field at cell (i,j,k) along the sigma layer (index space chain
    // rule), and of the four horizontal neighbour values
    static void grad5(lexer*, int i, int j, double fE, double fW, double fN, double fS, double &fx, double &fy);

    // gradient from limited index derivatives f_xi, f_eta (per index) of cell (i,j)
    static void grad_idx(lexer*, int i, int j, double fxi, double feta, double &fx, double &fy);

    // limited gradient (the limiter of nhflow_gradient, A 524 2): the limited derivatives along the
    // i- and j-lines of the cell centres, the gradient that has them as directional derivatives
    // (identical to the Cartesian one on Cartesian nodes); fE..fS the neighbours of fC
    static void grad_lim(lexer*, int i, int j, double fC, double fE, double fW, double fN, double fS, double &fx, double &fy);
    static void grad_lim(lexer*, slice&, int i, int j, double &fx, double &fy);
    
    // compact Laplacian of a slice at cell (i,j) (orthogonal part + non-orthogonal correction)
    static double lap(lexer*, slice&, int i, int j);

    // ghost cells of a velocity pair (u, v: U, V or UH, VH), after ghostcell::start4V has set them
    // component by component: the side walls (j-faces of solid cells, not the open y edges) mirror
    // the normal velocity u.n2 of the wall face; the river outflow (B 60) acts on the normal velocity
    static void pair_bc(lexer*, double *u, double *v);
    // a grid of the mesh refinement is deleted: its entry of the pairing of start4V
    static void pair_forget(const lexer*);
    // G 52 1 (multi-block grids): the i-walls act as the j-walls, the velocity pair mirrored at the face
    // normal (pair_bc), the other fields of start4V symmetric (iwall_sym); iwall: side 0 (-i) / 1 (+i) a wall
    static bool iwall(lexer*, int side);
    static void iwall_sym(lexer*, double *f);

    // fluxes through the boundary faces (after ghostcell::start1V / start2V, which set the Cartesian
    // ones): walls the pressure of the ghost cell along the face normal, in- and outflow faces and
    // open edges the flux of the ghost cell's state through the face; comp 1 UH, 2 VH, 3 WH, 4 the
    // continuity flux
    static void bc_fluxes(lexer*, fdm_nhf*, int comp, double *Fx, double *Fy);

    // ---- non-hydrostatic pressure (P at the sigma levels, FIJK)
    // horizontal gradient of the layer pressure 0.5 (P(k) + P(k+1)) of cell (i,j,k)
    static void pgrad_layer(lexer*, const double *P, int i, int j, int k, double &px, double &py);
    // horizontal divergence of (U,V) at sigma level k (between cells k-1 and k) of column (i,j), as
    // the Poisson right-hand sides: the layers interpolated to the level, next to the bed (k==0) the
    // lower cell weighted down with the bed slope along the index direction (A 522)
    static double divh_level(lexer*, fdm_nhf*, const double *U, const double *V, int i, int j, int k);
    // 2 (sigx d/dx + sigy d/dy) dP/dsigma at level k (the cross terms of the sigma Laplacian)
    static double sigma_cross(lexer*, const double *P, int i, int j, int k);
    // compact Laplacian coefficients L/(A dn) of the faces to i+1, i-1, j+1, j-1 (Poisson rows)
    static void lapcoef(lexer*, int i, int j, double &cIp, double &cIm, double &cJp, double &cJm);
    // non-orthogonal flux of grad P at level k through the faces of column (i,j) over the area
    // (deferred correction of the compact rows), faces to neighbours with ok(i,j) only
    template<class OK>
    static double kflux_level(lexer *p, const double *P, int ii, int jj, int kk, OK ok);

    // ---- momentum diffusion (A 512 1, explicit; nhflow_ediff): the stress form of the Cartesian
    // kernels, 2 nu u_xx + nu u_yy + nu v_xy = nu Lap_h(u) + nu d/dx(u_x + v_y), with the compact
    // Laplacian (non-orthogonal part included) and the chain rule; vertical and sigma terms as there.
    // comp 1: UH into F, 2: VH, 3: WH
    static void ediff(lexer*, fdm_nhf*, int comp, double *F, double *UH, double *VH, double *WH);
    // the parts of it kept explicit by the implicit diffusion (A 512 2, nhflow_idiff): the non-orthogonal
    // part of the Laplacian, d/dx (d/dy) of the horizontal divergence (diff_div first, for u and v), the
    // vertical transpose stress and the sigma cross terms, per unit viscosity
    static void diff_div(lexer*, const double *UH, const double *VH);
    static double diff_explicit(lexer*, int comp, int i, int j, int k, const double *UH, const double *VH, const double *WH);
    // the non-orthogonal part of the horizontal Laplacian of a cell field at layer k (faces to fluid cells)
    static double kflux3(lexer*, const double *f, int i, int j, int k);
    // horizontal gradient of the sigma derivative (f(k+1) - f(k-1))/(dsig) of a cell field at layer k
    static void grad_dsig(lexer*, const double *f, int i, int j, int k, double &qx, double &qy);
    
    // width of cell (i,j) across the i-faces and across the j-faces (time step)
    static double width1(lexer*, int i, int j);
    static double width2(lexer*, int i, int j);

    // the non-orthogonal part k = S - d |S|^2/(S.d) of face i+1/2 (dir 1) or j+1/2 (dir 2)
    static void kvec(lexer*, int dir, int i, int j, double &kx, double &ky);
};

template<class OK>
double nhflow_curv::kflux_level(lexer *p, const double *P, int ii, int jj, int kk, OK ok)
{
    auto F = [&](int a, int b) { return P[(a-p->imin)*p->jmax*p->kmaxF + (b-p->jmin)*p->kmaxF + kk-p->kmin]; };
    auto g = [&](int a, int b, double &gx, double &gy) { grad5(p,a,b,F(a+1,b),F(a-1,b),F(a,b+1),F(a,b-1),gx,gy); };
    double s=0.0, gx0, gy0, gx, gy;
    g(ii,jj,gx0,gy0);
    const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};
    for(int dd=0; dd<4; ++dd)
    {
        const int a=ii+di[dd], b=jj+dj[dd];
        if(!ok(a,b))
        continue;
        const int dir = (dd<2) ? 1 : 2;
        const int fi = (dd==1) ? ii-1 : ii, fj = (dd==3) ? jj-1 : jj;
        double kx,ky;
        kvec(p,dir,fi,fj,kx,ky);
        g(a,b,gx,gy);
        const double sg = (dd==0 || dd==2) ? 1.0 : -1.0;
        s += sg*0.5*((gx0+gx)*kx + (gy0+gy)*ky);
    }
    return s/p->AREA2D[(ii-p->imin)*p->jmax + (jj-p->jmin)];
}

#endif
