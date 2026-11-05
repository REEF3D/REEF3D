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

#include"sflow_ediff.h"
#include"sflow_curv.h"
#include"lexer.h"
#include"fdm2D.h"
#include<cmath>

sflow_ediff::sflow_ediff(lexer* p)
{
}

sflow_ediff::~sflow_ediff()
{
}

void sflow_ediff::diff_u(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &UHdiff, slice &UH, slice &U, slice &V, slice &WL, double alpha)
{
    if(own_eddyv==1)
    algebraic_eddyv(p,b,U,V,WL);

    if(p->geo_curv==1)
    {
        SLICELOOP4
        UHdiff(i,j) = UH(i,j);
        curv_stress(p,b,U,V,WL,1);
        return;
    }

	SLICELOOP4
    {
    UHdiff(i,j) = UH(i,j);
    
    WETDRY
	b->F(i,j) += WL(i,j)*viscosity(p,b)*(2.0*dxx(p,U) + dyy(p,U) + dxy(p,V));
    }
}

void sflow_ediff::diff_v(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &VHdiff, slice &VH, slice &U, slice &V, slice &WL, double alpha)
{
    if(p->geo_curv==1)
    {
        SLICELOOP4
        VHdiff(i,j) = VH(i,j);
        curv_stress(p,b,U,V,WL,2);
        return;
    }

	SLICELOOP4
    {
    VHdiff(i,j) = VH(i,j);
    
    WETDRY
	b->G(i,j) += WL(i,j)*viscosity(p,b)*(dxx(p,V) + 2.0*dyy(p,V) + dxy(p,U))*p->y_dir;
    }
}

void sflow_ediff::diff_w(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &WHdiff, slice &WH, slice &U, slice &V, slice &W, slice &WL, double alpha)
{
	SLICELOOP4
    {
    WHdiff(i,j) = WH(i,j);
    
    WETDRY
	b->H(i,j) += WL(i,j)*viscosity(p,b)*(dxx(p,W) + dyy(p,W));
    }
}

void sflow_ediff::diff_scalar(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &f, double sig, double alpha)
{
}

double sflow_ediff::viscosity(lexer *p, fdm2D *b)
{
    double visc = p->W2 + b->eddyv(i,j);
    
    if(p->A246==2 && b->breaking(i,j)==1)
    visc += p->A250;
    
    return visc;
}

// dry neighbours (inside the domain, wet 0) are free-slip, zero normal gradient: the cell's own value is
// used instead of the dry cell's (U = 0 there would act as a no-slip wall at the shoreline); walls and
// open boundaries keep their ghost values
#define SFLOW_DRY(X) (p->flagslice4[X]>0 && p->wet[X]==0)

double sflow_ediff::dxx(lexer *p, slice &f)
{
    const double fp = SFLOW_DRY(Ip1J) ? f(i,j) : f(i+1,j);
    const double fm = SFLOW_DRY(Im1J) ? f(i,j) : f(i-1,j);
    
    return ((fp - f(i,j))/p->DXP[IP] - (f(i,j) - fm)/p->DXP[IM1])/p->DXN[IP];
}

double sflow_ediff::dyy(lexer *p, slice &f)
{
    if(p->j_dir==0)
    return 0.0;
    
    const double fp = SFLOW_DRY(IJp1) ? f(i,j) : f(i,j+1);
    const double fm = SFLOW_DRY(IJm1) ? f(i,j) : f(i,j-1);
    
    return ((fp - f(i,j))/p->DYP[JP] - (f(i,j) - fm)/p->DYP[JM1])/p->DYN[JP];
}

double sflow_ediff::dxy(lexer *p, slice &f)
{
    if(p->j_dir==0)
    return 0.0;
    
    if(p->flagslice4[Ip1Jp1]<0 || p->flagslice4[Ip1Jm1]<0 || p->flagslice4[Im1Jp1]<0 || p->flagslice4[Im1Jm1]<0)
    return 0.0;
    
    if(p->wet[Ip1Jp1]==0 || p->wet[Ip1Jm1]==0 || p->wet[Im1Jp1]==0 || p->wet[Im1Jm1]==0)   // no cross term next to dry cells
    return 0.0;
    
    return (f(i+1,j+1) - f(i+1,j-1) - f(i-1,j+1) + f(i-1,j-1))/((p->DXP[IP]+p->DXP[IM1])*(p->DYP[JP]+p->DYP[JM1]));
}

// ---------------------------------------------------------------------------------------------------------------------
// curvilinear grids (G 50): the depth-integrated stress tau = nu h (grad u + grad u^T) through the faces of a cell over
// its area (finite volume), nu h of a face the mean of its two cells.  Per face with the face vector S (length L,
// sflow_curv::normal) and the cell centres P, N:
//   x-momentum  nu h ( du/dn L  +  k.grad u  +  u_x S_x + v_x S_y )
// du/dn L = (u_N - u_P) L/dn (the orthogonal part, sflow_curv::lapcoef), k the non-orthogonal part of S
// (sflow_curv::kvec) and the transposed part from the mean of the two cell gradients.  A face to a wall or a dry cell
// carries no stress (free slip), as the dry neighbours of the Cartesian operator; the cell gradients use the cell's own
// value for such neighbours.  On a Cartesian grid: the Laplacian of the 5-point stencil and the transposed part
// central, nu h (2 u_xx + u_yy + v_xy) as the Cartesian operator to second order.
// Both cells of a face compute the same flux: the stress is conserved across the faces of the blocks and of the
// patches of one level (the ghost cells hold the neighbour's state).
void sflow_ediff::curv_stress(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL, int comp)
{
    auto fluid = [&](int a, int d) { const int q=(a-p->imin)*p->jmax+(d-p->jmin); return p->flagslice4[q]>0 && p->wet[q]==1; };

    // cell gradient, non-fluid neighbours replaced by the cell's value
    auto cgrad = [&](slice &f, int a, int d, double &gx, double &gy)
    {
        const double fc = f(a,d);
        sflow_curv::grad5(p,a,d, fluid(a+1,d) ? f(a+1,d) : fc, fluid(a-1,d) ? f(a-1,d) : fc,
                                  fluid(a,d+1) ? f(a,d+1) : fc, fluid(a,d-1) ? f(a,d-1) : fc, gx, gy);
    };

    auto nuh = [&](int a, int d) { return (p->W2 + b->eddyv(a,d) + ((p->A246==2 && b->breaking(a,d)==1) ? p->A250 : 0.0))*WL(a,d); };

    slice &f = (comp==1) ? U : V;
    const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};

    SLICELOOP4
    WETDRY
    {
        double s = 0.0, uxP, uyP, vxP, vyP;
        cgrad(U,i,j,uxP,uyP);
        cgrad(V,i,j,vxP,vyP);

        for(int d=0; d<4; ++d)
        {
            const int a=i+di[d], e=j+dj[d];

            if(!fluid(a,e))
            continue;

            const int dir = (d<2) ? 1 : 2;
            const int fi = (d==1) ? i-1 : i, fj = (d==3) ? j-1 : j;
            const double sg = (d==0 || d==2) ? 1.0 : -1.0;   // the face vector points from P to N

            double nx, ny, len, kx, ky, uxN, uyN, vxN, vyN;
            sflow_curv::normal(p,dir,fi,fj,nx,ny,len);
            sflow_curv::kvec(p,dir,fi,fj,kx,ky);
            const double dn = sflow_curv::dnorm(p,dir,fi,fj);
            cgrad(U,a,e,uxN,uyN);
            cgrad(V,a,e,vxN,vyN);

            const double Sx = sg*nx*len, Sy = sg*ny*len;
            kx *= sg; ky *= sg;

            const double gx = 0.5*((comp==1 ? uxP : vxP) + (comp==1 ? uxN : vxN));
            const double gy = 0.5*((comp==1 ? uyP : vyP) + (comp==1 ? uyN : vyN));

            // transposed part: x-momentum u_x S_x + v_x S_y, y-momentum u_y S_x + v_y S_y
            const double tr = (comp==1) ? 0.5*(uxP+uxN)*Sx + 0.5*(vxP+vxN)*Sy
                                        : 0.5*(uyP+uyN)*Sx + 0.5*(vyP+vyN)*Sy;

            s += 0.5*(nuh(i,j) + nuh(a,e))*((f(a,e) - f(i,j))*len/dn + gx*kx + gy*ky + tr);
        }

        const double A = sflow_curv::area(p,i,j);

        if(comp==1)
        b->F(i,j) += s/A;
        else
        b->G(i,j) += s/A;
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// the algebraic eddy viscosity on a grid that computes it itself (own_eddyv): the cells and their first ghost layer,
// so that both cells of a face of two grids see the same viscosity
//   A 260 3: Prandtl mixing length  nu_t = (0.4 A261 h)^2 |S|,  |S| = sqrt(2 u_x^2 + 2 v_y^2 + (u_y + v_x)^2)
//   A 260 4: parabolic              nu_t = A262 u* h,  u* = sqrt(cf) |u|, cf = g n^2 / h^(1/3), n = ks^(1/6)/20
// (sflow_turb_prandtl, sflow_turb_parabolic; the gradients central, sflow_curv::grad5 on curvilinear grids)
void sflow_ediff::algebraic_eddyv(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    if(p->A260!=3 && p->A260!=4)
    return;

    auto fluid = [&](int a, int d) { const int q=(a-p->imin)*p->jmax+(d-p->jmin); return p->flagslice4[q]>0 && p->wet[q]==1; };

    auto cgrad = [&](slice &f, int a, int d, double &gx, double &gy)
    {
        const double fc = f(a,d);
        const double fE = fluid(a+1,d) ? f(a+1,d) : fc, fW = fluid(a-1,d) ? f(a-1,d) : fc;
        const double fN = fluid(a,d+1) ? f(a,d+1) : fc, fS = fluid(a,d-1) ? f(a,d-1) : fc;

        if(p->geo_curv==1)
        sflow_curv::grad5(p,a,d,fE,fW,fN,fS,gx,gy);

        else
        {
        gx = (fE - fW)/(p->DXP[a+marge] + p->DXP[a-1+marge]);
        gy = (p->j_dir==1) ? (fN - fS)/(p->DYP[d+marge] + p->DYP[d-1+marge]) : 0.0;
        }
    };

    for(int a=-1; a<p->knox+1; ++a)
    for(int d=-1; d<p->knoy+1; ++d)
    {
        const int q=(a-p->imin)*p->jmax+(d-p->jmin);

        if(p->flagslice4[q]<=0)
        continue;

        if(p->wet[q]==0)
        {
        b->eddyv(a,d) = 0.0;
        continue;
        }

        const double h = MAX(WL(a,d),0.0);

        if(p->A260==3)
        {
            double ux, uy, vx, vy;
            cgrad(U,a,d,ux,uy);
            cgrad(V,a,d,vx,vy);
            const double l = 0.4*p->A261*h;
            b->eddyv(a,d) = l*l*sqrt(2.0*ux*ux + 2.0*vy*vy + (uy+vx)*(uy+vx));
        }

        if(p->A260==4)
        {
            const double n = pow(b->ks(a,d),1.0/6.0)/20.0;
            const double cf = n*n*fabs(p->W22)/pow(MAX(h,1.0e-6),1.0/3.0);
            b->eddyv(a,d) = p->A262*sqrt(cf*(U(a,d)*U(a,d) + V(a,d)*V(a,d)))*h;
        }
    }
}
