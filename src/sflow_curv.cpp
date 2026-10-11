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

#include"sflow_curv.h"
#include"lexer.h"
#include"fdm2D.h"
#include"slice.h"
#include<cmath>

void sflow_curv::normal(lexer *p, int dir, int ii, int jj, double &nx, double &ny, double &len)
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

void sflow_curv::check(lexer *p)
{
    int bad = 0;

    auto refuse = [&](bool cond, const char *what)
    {
        if(cond)
        {
            if(p->mpirank==0)
            cout<<"!!! curvilinear grid (G 50): "<<what<<" is not available yet !!!"<<endl;
            ++bad;
        }
    };

    refuse(p->j_dir!=1, "a 1D grid (j_dir 0)");
    refuse(p->A220!=0, "non-hydrostatic pressure and Boussinesq (A 220 must be 0)");
    refuse(p->A212>0 || p->A260>0, "horizontal diffusion and turbulence (A 212, A 260)");
    refuse(p->G1>0, "mesh refinement (G 1)");
    refuse(p->B98>=1 || p->B90>=1, "wave generation and absorption (B 90, B 98)");
    refuse(p->B440>0 || p->B441>0 || p->B442>0, "patch boundary conditions (B 440-442)");
    refuse(p->S10>0, "sediment transport (S 10)");
    refuse(p->X10>0, "6DOF bodies (X 10)");
    refuse(p->A750==1, "the SEASTATE coupling (A 750)");
    refuse(p->I11==1, "the potential flow initialisation (I 11)");

    if(bad>0)
    exit(1);
}

// normal velocity at the faces: Us/Un from (Us,Vs)/(Un,Vn), Ve/Vw from (Ue,Ve)/(Uw,Vw)
void sflow_curv::normal_velocities(lexer *p, fdm2D *b)
{
    double nx,ny,len;

    SLICELOOP1
    {
    normal(p,1,i,j,nx,ny,len);
    b->Us(i,j) = b->Us(i,j)*nx + b->Vs(i,j)*ny;
    b->Un(i,j) = b->Un(i,j)*nx + b->Vn(i,j)*ny;
    }

    SLICELOOP2
    {
    normal(p,2,i,j,nx,ny,len);
    b->Ve(i,j) = b->Ue(i,j)*nx + b->Ve(i,j)*ny;
    b->Vw(i,j) = b->Uw(i,j)*nx + b->Vw(i,j)*ny;
    }
}

void sflow_curv::flux(lexer *p, fdm2D *b, int ipol)
{
    const double g = fabs(p->W22);
    double nx,ny,len,Ps,Pn;

    // i-faces (normal velocity in Us/Un)
    SLICELOOP1
    {
    normal(p,1,i,j,nx,ny,len);

    Ps = 0.5*g*b->ETAs(i,j)*b->ETAs(i,j) + g*b->ETAs(i,j)*b->dfx(i,j);
    Pn = 0.5*g*b->ETAn(i,j)*b->ETAn(i,j) + g*b->ETAn(i,j)*b->dfx(i,j);

        if(ipol==1)
        {
        b->Fs(i,j) = b->UHs(i,j)*b->Us(i,j) + Ps*nx;
        b->Fn(i,j) = b->UHn(i,j)*b->Un(i,j) + Pn*nx;
        }

        if(ipol==2)
        {
        b->Fs(i,j) = b->VHs(i,j)*b->Us(i,j) + Ps*ny;
        b->Fn(i,j) = b->VHn(i,j)*b->Un(i,j) + Pn*ny;
        }

        if(ipol==3)
        {
        b->Fs(i,j) = b->WHs(i,j)*b->Us(i,j);
        b->Fn(i,j) = b->WHn(i,j)*b->Un(i,j);
        }

        if(ipol==4)
        {
        b->Fs(i,j) = b->UHs(i,j)*nx + b->VHs(i,j)*ny;
        b->Fn(i,j) = b->UHn(i,j)*nx + b->VHn(i,j)*ny;
        }
    }

    // j-faces (normal velocity in Ve/Vw)
    SLICELOOP2
    {
    normal(p,2,i,j,nx,ny,len);

    Ps = 0.5*g*b->ETAe(i,j)*b->ETAe(i,j) + g*b->ETAe(i,j)*b->dfy(i,j);
    Pn = 0.5*g*b->ETAw(i,j)*b->ETAw(i,j) + g*b->ETAw(i,j)*b->dfy(i,j);

        if(ipol==1)
        {
        b->Fe(i,j) = b->UHe(i,j)*b->Ve(i,j) + Ps*nx;
        b->Fw(i,j) = b->UHw(i,j)*b->Vw(i,j) + Pn*nx;
        }

        if(ipol==2)
        {
        b->Fe(i,j) = b->VHe(i,j)*b->Ve(i,j) + Ps*ny;
        b->Fw(i,j) = b->VHw(i,j)*b->Vw(i,j) + Pn*ny;
        }

        if(ipol==3)
        {
        b->Fe(i,j) = b->WHe(i,j)*b->Ve(i,j);
        b->Fw(i,j) = b->WHw(i,j)*b->Vw(i,j);
        }

        if(ipol==4)
        {
        b->Fe(i,j) = b->UHe(i,j)*nx + b->VHe(i,j)*ny;
        b->Fw(i,j) = b->UHw(i,j)*nx + b->VHw(i,j)*ny;
        }
    }
}

void sflow_curv::flux_bc(lexer *p, fdm2D *b, int ipol, slice &Fx, slice &Fy, int inflow, int outflow)
{
    const double g = fabs(p->W22);
    double wc,sw,un,nx,ny,len,P,Fnrm;

    // normal-momentum flux of a wall face, as Cartesian components (mass, WH: none)
    //   side -1: the wall is on the -n side of the cell (interior right state), +1: on the +n side
    auto wall = [&](slice &F, int fi, int fj, int dir, double dface, double side)
    {
        normal(p,dir,fi,fj,nx,ny,len);
        un = b->U(i,j)*nx + b->V(i,j)*ny;
        sw = fabs(un) + sqrt(g*wc);
        P = 0.5*g*b->eta(i,j)*b->eta(i,j) + g*b->eta(i,j)*dface;
        Fnrm = P + wc*un*un + side*sw*wc*un;

        if(ipol==1)
        F(fi,fj) = Fnrm*nx;

        else
        if(ipol==2)
        F(fi,fj) = Fnrm*ny;

        else
        F(fi,fj) = 0.0;
    };

    SLICELOOP4
    {
    wc = MAX(b->eta(i,j) + b->depth(i,j), 0.0);

    if(p->flagslice4[Im1J]<0)
    wall(Fx,i-1,j,1,b->dfx(i-1,j),-1.0);

    if(p->flagslice4[Ip1J]<0)
    wall(Fx,i,j,1,b->dfx(i,j),1.0);

    if(p->flagslice4[IJm1]<0)
    wall(Fy,i,j-1,2,b->dfy(i,j-1),-1.0);

    if(p->flagslice4[IJp1]<0)
    wall(Fy,i,j,2,b->dfy(i,j),1.0);
    }

    // in- and outflow: physical flux of the ghost state through the i-face
    auto ghost = [&](int gi, int gj, int fi, int fj)
    {
        const double eg = b->eta(gi,gj);
        const double wl = MAX(eg + b->depth(gi,gj), 0.0);

        normal(p,1,fi,fj,nx,ny,len);
        un = b->U(gi,gj)*nx + b->V(gi,gj)*ny;
        P = 0.5*g*eg*eg + g*eg*b->dfx(fi,fj);

        if(ipol==1)
        Fx(fi,fj) = wl*b->U(gi,gj)*un + P*nx;

        if(ipol==2)
        Fx(fi,fj) = wl*b->V(gi,gj)*un + P*ny;

        if(ipol==3)
        Fx(fi,fj) = wl*b->W(gi,gj)*un;

        if(ipol==4)
        Fx(fi,fj) = wl*un;
    };

    if(inflow==1)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];

    ghost(i-1,j,i-1,j);
    }

    if(outflow==1)
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];

    ghost(i+1,j,i,j);
    }
}

void sflow_curv::divergence(lexer *p, fdm2D *b, slice &f, slice &Fx, slice &Fy)
{
    SLICELOOP4
    WETDRY
    f(i,j) -= (Fx(i,j)*p->LEN1[IJ] - Fx(i-1,j)*p->LEN1[Im1J]
             + Fy(i,j)*p->LEN2[IJ] - Fy(i,j-1)*p->LEN2[IJm1])/p->AREA2D[IJ];
}

void sflow_curv::pgrad(lexer *p, fdm2D *b, slice &eta, slice &f, int comp)
{
    const double g = fabs(p->W22);
    const double *S1 = (comp==1) ? p->SX1 : p->SY1;
    const double *S2 = (comp==1) ? p->SX2 : p->SY2;

    SLICELOOP4
    WETDRY
    f(i,j) += g*eta(i,j)*(b->dfx(i,j)*S1[IJ] - b->dfx(i-1,j)*S1[Im1J]
                        + b->dfy(i,j)*S2[IJ] - b->dfy(i,j-1)*S2[IJm1])/p->AREA2D[IJ];
}

double sflow_curv::sigma(lexer *p, fdm2D *b, int ii, int jj, double c)
{
    const int q = (ii-p->imin)*p->jmax + (jj-p->jmin);
    const int qw = q - p->jmax, qs = q - 1;
    const double u = b->U(ii,jj), v = b->V(ii,jj);

    // mean normal vectors of the opposite faces of the cell
    const double s1x = 0.5*(p->SX1[q]+p->SX1[qw]), s1y = 0.5*(p->SY1[q]+p->SY1[qw]);
    const double s2x = 0.5*(p->SX2[q]+p->SX2[qs]), s2y = 0.5*(p->SY2[q]+p->SY2[qs]);
    const double l1 = sqrt(s1x*s1x + s1y*s1y), l2 = sqrt(s2x*s2x + s2y*s2y);
    const double A = p->AREA2D[q];

    double s = 0.0;

    if(l1>0.0)
    s += (fabs(u*s1x + v*s1y)/l1 + c)*l1/A;

    if(l2>0.0)
    s += (fabs(u*s2x + v*s2y)/l2 + c)*l2/A;

    return s;
}

void sflow_curv::wall_mirror(lexer *p, slice &u, slice &v, int ii, int jj, int cs)
{
    double nx,ny,len,un;

    // the wall face of the boundary cell
    if(cs==1) normal(p,1,ii-1,jj,nx,ny,len);
    if(cs==4) normal(p,1,ii,jj,nx,ny,len);
    if(cs==2) normal(p,2,ii,jj,nx,ny,len);
    if(cs==3) normal(p,2,ii,jj-1,nx,ny,len);

    for(int qq=0; qq<3; ++qq)
    {
        int ai=ii, aj=jj, gi=ii, gj=jj;

        if(cs==1) { ai=ii+qq; gi=ii-qq-1; }
        if(cs==4) { ai=ii-qq; gi=ii+qq+1; }
        if(cs==2) { aj=jj-qq; gj=jj+qq+1; }
        if(cs==3) { aj=jj+qq; gj=jj-qq-1; }

        un = u(ai,aj)*nx + v(ai,aj)*ny;

        u(gi,gj) = u(ai,aj) - 2.0*un*nx;
        v(gi,gj) = v(ai,aj) - 2.0*un*ny;
    }
}
