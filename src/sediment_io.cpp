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

#include"sediment_f.h"
#include"sediment_curv.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"bedshear.h"
#include"suspended.h"

double sediment_f::bedshear_point(lexer *p, ghostcell *pgc)
{
	return s->tau_eff(i,j);
}

void sediment_f::fill_bedk(lexer *p, fdm *a,ghostcell *pgc)
{
    SLICELOOP4
    s->bedk(i,j)=0;
    
    LOOP
    if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
    s->bedk(i,j)=k+1;
}

void sediment_f::fill_PQ_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
    double zval,xip,yip;
    
    SLICELOOP1
    {
    k=s->bedk(i,j);
    
    xip= p->XN[IP1];
	yip= p->YP[JP];
    zval = 0.5*(s->bedzh(i,j)+s->bedzh(i+1,j)) + 1.6*p->DZN[KP];
    
    s->P(i,j) = a->P(i,j) = p->ccipol1_a(a->u,xip,yip,zval);
    }
    
    SLICELOOP2
    {
    k=s->bedk(i,j);
    
    xip= p->XP[IP];
	yip= p->YN[JP1];
    zval = 0.5*(s->bedzh(i,j)+s->bedzh(i,j+1)) + 1.6*p->DZN[KP];
    
    s->Q(i,j) = a->Q(i,j)  = p->ccipol2_a(a->v,xip,yip,zval);
    }
    
    SLICELOOP4
    {
    k=s->bedk(i,j);
    
    s->ro(i,j) = a->ro(i,j,k);
    s->waterlevel(i,j) = a->WL(i,j);
    }
    
    pgc->gcsl_start1(p,s->P,10);
	pgc->gcsl_start2(p,s->Q,11);
    pgc->gcsl_start4(p,s->ro,1);
    pgc->gcsl_start4(p,s->waterlevel,1);
}

void sediment_f::fill_PQ_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    k=0;
    
    // curvilinear grid (G 50): the cell velocities of the bottom layer (with their ghost cells) and
    // the normal velocities of the faces (sediment_curv.h)
    if(p->geo_curv==1)
    {
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
        const int c = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax - p->kmin;
        s->ub(ii,jj) = d->U[c];
        s->vb(ii,jj) = d->V[c];
    }
    
    sediment_curv::face_normals(p,pgc,s,0);
    
    SLICELOOP4
    s->waterlevel(i,j) = d->WL(i,j);
    
    pgc->gcsl_start4(p,s->waterlevel,1);
    return;
    }
    
    SLICELOOP1
    s->P(i,j) = 0.5*(d->U[IJK] + d->U[Ip1JK]);
    
    SLICELOOP2
    s->Q(i,j) = 0.5*(d->V[IJK] + d->V[IJp1K]);
    
    SLICELOOP4
    s->waterlevel(i,j) = d->WL(i,j);
    
    pgc->gcsl_start1(p,s->P,10);
	pgc->gcsl_start2(p,s->Q,11);  
    pgc->gcsl_start4(p,s->waterlevel,1);
    
    // open domain edges (in- and outflow, IO 1, 2): the face velocities of the bottom layer from the boundary
    // ghost cells, so that the Exner fluxes cross them (the ghost faces were left at zero: the bed load could
    // neither enter nor leave the domain); closed edges stay closed
    k=0;
    
    if(p->origin_i==0)
    for(j=0; j<p->knoy; ++j)
    {
    i=0;
    if(p->flag4[IJK]>0 && p->IO[Im1JK]>0)
    s->P(-1,j) = 0.5*(d->U[Im1JK] + d->U[IJK]);
    }
    
    if(p->origin_i+p->knox==p->gknox)
    for(j=0; j<p->knoy; ++j)
    {
    i=p->knox-1;
    if(p->flag4[IJK]>0 && p->IO[Ip1JK]>0)
    s->P(i,j) = 0.5*(d->U[IJK] + d->U[Ip1JK]);
    }
    
    if(p->j_dir==1 && p->origin_j==0)
    for(i=0; i<p->knox; ++i)
    {
    j=0;
    if(p->flag4[IJK]>0 && p->IO[IJm1K]>0)
    s->Q(i,-1) = 0.5*(d->V[IJm1K] + d->V[IJK]);
    }
    
    if(p->j_dir==1 && p->origin_j+p->knoy==p->gknoy)
    for(i=0; i<p->knox; ++i)
    {
    j=p->knoy-1;
    if(p->flag4[IJK]>0 && p->IO[IJp1K]>0)
    s->Q(i,j) = 0.5*(d->V[IJK] + d->V[IJp1K]);
    }
}

void sediment_f::fill_PQ_sflow(lexer *p, fdm2D *b,ghostcell *pgc,slice &P, slice &Q)
{
    // curvilinear grid (G 50): the cell velocities and the normal velocities of the faces (sediment_curv.h)
    if(p->geo_curv==1)
    {
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    {
        s->ub(ii,jj) = b->U(ii,jj);
        s->vb(ii,jj) = b->V(ii,jj);
    }
    
    sediment_curv::face_normals(p,pgc,s,1);
    inflow_PQ(p);
    
    SLICELOOP4
    s->waterlevel(i,j) = b->hp(i,j);
    
    pgc->gcsl_start4(p,s->waterlevel,1);
    return;
    }
    
    SLICELOOP1
    s->P(i,j) = P(i,j);
    
    SLICELOOP2
    s->Q(i,j) = Q(i,j);
    
    // water depth (non-equilibrium S33 4, van Rijn reference concentration, Engelund-Hansen)
    SLICELOOP4
    s->waterlevel(i,j) = b->hp(i,j);
    
    pgc->gcsl_start1(p,s->P,10);
	pgc->gcsl_start2(p,s->Q,11);
    pgc->gcsl_start4(p,s->waterlevel,1);
    
    // domain edge faces: the face velocities of SFLOW (in- and outflow; zero at walls), the ghost faces of
    // gcsl_start1/2 were zero, so the bed load could neither enter nor leave the domain and the first cell
    // had half its velocity in the bed shear stress
    if(p->origin_i==0)
    for(j=0; j<p->knoy; ++j)
    s->P(-1,j) = P(-1,j);
    
    if(p->origin_i+p->knox==p->gknox)
    for(j=0; j<p->knoy; ++j)
    s->P(p->knox-1,j) = P(p->knox-1,j);
    
    if(p->j_dir==1 && p->origin_j==0)
    for(i=0; i<p->knox; ++i)
    s->Q(i,-1) = Q(i,-1);
    
    if(p->j_dir==1 && p->origin_j+p->knoy==p->gknoy)
    for(i=0; i<p->knox; ++i)
    s->Q(i,p->knoy-1) = Q(i,p->knoy-1);
}

void sediment_f::inflow_PQ(lexer *p)
{
    // SFLOW on a curvilinear grid (the Cartesian grid takes the face velocities of the flow at the domain edges,
    // fill_PQ_sflow): the face velocities of the sediment at the inflow faces (gcslin; the ghost faces are not filled
    // by the boundary conditions of P, Q and stayed 0): the velocity of the first cell, so that its bed shear, bed load and
    // transport direction are those of the flow (with 0 the first row had half the velocity and a fraction of the
    // bed load, the second row eroded twice as fast).  The faces stay closed for the sediment (dfbed_inflow).
    if(p->geo_curv!=1)
    return;

    double nx,ny,len;

    for(int n=0; n<p->gcslin_count; ++n)
    {
        const int ii = p->gcslin[n][0], jj = p->gcslin[n][1], cs = p->gcslin[n][3];
        const int dir = (cs==1 || cs==4) ? 1 : 2;
        const int fi = (cs==1) ? ii-1 : ii;
        const int fj = (cs==3) ? jj-1 : jj;

        {
            // the cells beyond the face take the velocity of the first cell: the tangential part of the face and the
            // curvature of the streamlines there (the ghost velocities of the native inflow and of the open face of a
            // block differ)
            const int di = (cs==1) ? -1 : (cs==4 ? 1 : 0);
            const int dj = (cs==3) ? -1 : (cs==2 ? 1 : 0);
            for(int q=1; q<=p->margin; ++q)
            {
            const int ia = ii + q*di, ja = jj + q*dj;
                if(ia>=p->imin && ia<p->imin+p->imax && ja>=p->jmin && ja<p->jmin+p->jmax)
                {
                s->ub(ia,ja) = s->ub(ii,jj);
                s->vb(ia,ja) = s->vb(ii,jj);
                }
            }

            sediment_curv::normal(p,dir,fi,fj,nx,ny,len);
            const double un = (p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]==0) ? 0.0 : s->ub(ii,jj)*nx + s->vb(ii,jj)*ny;

            if(dir==1)
            s->P(fi,fj) = un;
            else
            s->Q(fi,fj) = un;
        }
    }
}

double sediment_f::qbeval(int ii, int jj)
{
    double val;

    val=s->qbe(ii,jj);

    return val;
}

void sediment_f::qbeget(int ii, int jj, double val)
{
    s->qbe(ii,jj)=val;
}

double sediment_f::bedzhval(int ii, int jj)
{
    double val;

    val=s->bedzh(ii,jj);

    return val;
}

void sediment_f::ctimesave(lexer *p, fdm* a)
{
    psusp->ctimesave(p,a);
}
