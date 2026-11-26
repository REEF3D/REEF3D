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

#include"sflow_turb_io.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_curv.h"

sflow_turb_io::sflow_turb_io(lexer* p) : kin(p), eps(p)
{

}

sflow_turb_io::~sflow_turb_io()
{
}

void sflow_turb_io::print_2D(lexer *p, fdm2D *b, ghostcell *pgc, std::ostream &result)
{
    iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(kin));
	result.write((char*)&ffn, sizeof (float));
	}
    
    
    iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(eps));
	result.write((char*)&ffn, sizeof (float));
	}
    
}

void sflow_turb_io::kinget(int ii, int jj, double val)
{
    kin(ii,jj)=val;
}

void sflow_turb_io::epsget(int ii, int jj, double val)
{
    eps(ii,jj)=val;
}
    
double sflow_turb_io::kinval(int ii, int jj)
{
    val=kin(ii,jj);

    return val;
}
    
double sflow_turb_io::epsval(int ii, int jj)
{
    val=eps(ii,jj);

    return val;
}
    
void sflow_turb_io::name_pvtp(lexer *p, fdm2D *b, ghostcell *pgc,ofstream &result)
{
    result<<"<PDataArray type=\"Float32\" Name=\"kin\"/>\n";
	
	if(p->A260==1)
	result<<"<PDataArray type=\"Float32\" Name=\"epsilon\"/>\n";
	if(p->A260==2)
    result<<"<PDataArray type=\"Float32\" Name=\"omega\"/>\n";
}

void sflow_turb_io::name_vtp(lexer *p, fdm2D *b, ghostcell *pgc,std::ostream &result, int *offset, int &n)
{
    result<<"<DataArray type=\"Float32\" Name=\"kin\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
	if(p->A260==1)
	result<<"<DataArray type=\"Float32\" Name=\"epsilon\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
	if(p->A260==2)
    result<<"<DataArray type=\"Float32\" Name=\"omega\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
}
    
void sflow_turb_io::offset_ParaView_2D(lexer *p, int *offset, int &n)
{
    offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
}

void sflow_turb_io::timesource(lexer* p, fdm2D *b, slice &fn)
{
    int count=0;
    SLICELOOP4
    {
        b->M.p[count] += 1.0/p->dt;

        b->rhsvec.V[count] += b->L(i,j) + fn(i,j)/p->dt;

	++count;
    }
}

void sflow_turb_io::wall_bc(lexer* p, fdm2D *b, slice &f)
{
    // SFLOW side walls are free-slip in the momentum equations (sflow_momentum_func): zero gradient of k and
    // epsilon / omega there and at dry neighbours, without wall functions
    n=0;
    SLICELOOP4
    {
        if(p->wet[IJ]==1)
        {
        if(p->flagslice4[Im1J]<0 || p->wet[Im1J]==0)
        {
        b->rhsvec.V[n] -= b->M.s[n]*f(i,j);
        b->M.s[n] = 0.0;
        }

        if(p->flagslice4[Ip1J]<0 || p->wet[Ip1J]==0)
        {
        b->rhsvec.V[n] -= b->M.n[n]*f(i,j);
        b->M.n[n] = 0.0;
        }

        if(p->flagslice4[IJm1]<0 || p->wet[IJm1]==0)
        {
        b->rhsvec.V[n] -= b->M.e[n]*f(i,j);
        b->M.e[n] = 0.0;
        }

        if(p->flagslice4[IJp1]<0 || p->wet[IJp1]==0)
        {
        b->rhsvec.V[n] -= b->M.w[n]*f(i,j);
        b->M.w[n] = 0.0;
        }
        }

    ++n;
    }

    // dry cells: identity rows, f = 0
    n=0;
    SLICELOOP4
    {
        if(p->wet[IJ]==0)
        {
        b->M.p[n] = 1.0;
        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;
        b->M.w[n] = 0.0;
        b->M.e[n] = 0.0;
        b->rhsvec.V[n] = 0.0;
        }
    ++n;
    }
}

// Curvilinear grid (G 50 1/2): k-epsilon and k-omega (A 260 1/2) on the faces of the logical grid.
// Convection: implicit first-order upwind in advective form, as sflow_ifou on the Cartesian grid, with the face
// fluxes (u.n) L of the normal velocity (mean of the two cells, the own cell at walls and dry neighbours) over the
// cell area A:
//   M.p = (F_w+ - F_e- + F_s+ - F_n-)/A,  M.s = -F_w+/A,  M.n = F_e-/A,  M.e = -F_s+/A,  M.w = F_n-/A
// (F_w, F_s the flux in through the faces i-1/2, j-1/2, F_e, F_n out through i+1/2, j+1/2; + and - the positive and
// negative parts). The rows sum to zero and the off-diagonals are not positive: k and epsilon / omega stay positive
// at any Courant number of the semi-implicit time steps. On a Cartesian grid (L/A = 1/dx) these are the rows of
// sflow_ifou with the face velocities averaged from the cells.
void sflow_turb_io::curv_convection(lexer *p, fdm2D *b)
{
    auto fluid = [&](int a, int d) { const int q=(a-p->imin)*p->jmax+(d-p->jmin); return p->flagslice4[q]>0 && p->wet[q]==1; };

    n=0;
    SLICELOOP4
    {
        auto flux = [&](int dir, int fi, int fj, int a, int d)
        {
            double nx,ny,len;
            sflow_curv::normal(p,dir,fi,fj,nx,ny,len);

            double u = b->U(i,j), v = b->V(i,j);

            if(fluid(a,d))
            {
            u = 0.5*(u + b->U(a,d));
            v = 0.5*(v + b->V(a,d));
            }

            return (u*nx + v*ny)*len;
        };

        const double A = sflow_curv::area(p,i,j);
        const double fe = flux(1,i,j,i+1,j), fw = flux(1,i-1,j,i-1,j);
        const double fn = flux(2,i,j,i,j+1), fs = flux(2,i,j-1,i,j-1);

        const double wp = MAX(fw,0.0), em = MIN(fe,0.0), sp = MAX(fs,0.0), nm = MIN(fn,0.0);

        b->M.p[n] = (wp - em + sp - nm)/A;
        b->M.s[n] = -wp/A;
        b->M.n[n] =  em/A;
        b->M.e[n] = -sp/A;
        b->M.w[n] =  nm/A;
    ++n;
    }
}

// diffusion nu + nu_t/sig with the face mean of the two cells (as sflow_idiff::diff_scalar), the orthogonal part of
// the face gradients (sflow_curv::lapcoef: L/(A dn)); walls and dry neighbours by wall_bc
void sflow_turb_io::curv_diffusion(lexer *p, fdm2D *b, double sig)
{
    n=0;
    SLICELOOP4
    {
        double cE,cW,cN,cS;
        sflow_curv::lapcoef(p,i,j,cE,cW,cN,cS);

        const double vc = p->W2 + b->eddyv(i,j)/sig;
        const double dE = 0.5*(vc + p->W2 + b->eddyv(i+1,j)/sig);
        const double dW = 0.5*(vc + p->W2 + b->eddyv(i-1,j)/sig);
        const double dN = 0.5*(vc + p->W2 + b->eddyv(i,j+1)/sig);
        const double dS = 0.5*(vc + p->W2 + b->eddyv(i,j-1)/sig);

        b->M.p[n] += cE*dE + cW*dW + cN*dN + cS*dS;
        b->M.n[n] -= cE*dE;
        b->M.s[n] -= cW*dW;
        b->M.w[n] -= cN*dN;
        b->M.e[n] -= cS*dS;
    ++n;
    }
}

// cell-centred velocity gradients by the chain rule with the cell centres (sflow_curv::grad5); walls (free-slip in
// SFLOW), open boundaries and dry neighbours zero-gradient (mirror), as the Cartesian production of k
void sflow_turb_io::curv_velgrad(lexer *p, fdm2D *b, int ii, int jj, double &ux, double &uy, double &vx, double &vy)
{
    auto fluid = [&](int a, int d) { const int q=(a-p->imin)*p->jmax+(d-p->jmin); return p->flagslice4[q]>0 && p->wet[q]==1; };

    auto g = [&](slice &f, double &fx, double &fy)
    {
        const double fc = f(ii,jj);
        const double fE = fluid(ii+1,jj) ? f(ii+1,jj) : fc, fW = fluid(ii-1,jj) ? f(ii-1,jj) : fc;
        const double fN = fluid(ii,jj+1) ? f(ii,jj+1) : fc, fS = fluid(ii,jj-1) ? f(ii,jj-1) : fc;
        sflow_curv::grad5(p,ii,jj,fE,fW,fN,fS,fx,fy);
    };

    g(b->U,ux,uy);
    g(b->V,vx,vy);
}
