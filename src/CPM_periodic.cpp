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


#include"CPM.h"
#include"lexer.h"
#include"ghostcell.h"

/*--------------------------------------------------------------------
periodic boundaries for the parcels (DIVEMesh C 21 x, C 22 y)

serial (periodic = 1): the kernel deposits into the ghost cells beyond the periodic sides,
pfold adds them to the cells on the opposite side before start4a_sum (which copies the
periodic ghost cells and would overwrite the deposits). Parcels crossing a periodic side are
moved by the period after the grid-limited step, so the limiter sees the move into the ghost
cell, which is the copy of the cell on the opposite side.

parallel (periodic = 2): the ranks at the periodic sides are neighbours, start4a_sum sums the
deposits and the parcels go to the neighbour by the regular exchange; the receiving rank moves
them by the period (part::xchange_fillback_flag).
--------------------------------------------------------------------*/

void CPM::pfold(lexer *p, field &f)
{
    if(perx==1)
    {
        for(j=-1;j<p->knoy+1;++j)
        for(k=-1;k<p->knoz+1;++k)
        {
            f(p->knox-1,j,k) += f(-1,j,k);
            f(0,j,k) += f(p->knox,j,k);
            f(-1,j,k) = 0.0;
            f(p->knox,j,k) = 0.0;
        }
    }
    
    if(pery==1)
    {
        for(i=-1;i<p->knox+1;++i)
        for(k=-1;k<p->knoz+1;++k)
        {
            f(i,p->knoy-1,k) += f(i,-1,k);
            f(i,0,k) += f(i,p->knoy,k);
            f(i,-1,k) = 0.0;
            f(i,p->knoy,k) = 0.0;
        }
    }
}

/*--------------------------------------------------------------------
deposits of the parcels in the ghost cells (kernel): solid volume, momentum, coupling terms

The kernel of a parcel within half a cell of a side of the subdomain deposits into the ghost cells
beyond that side. start4a_sum adds the ghost cells of a face to the cells of the neighbour subdomain.
A parcel next to two sides with neighbours (an edge of the subdomain, decompositions with two or three
directions split) also deposits into the ghost cells beyond the edge, which belong to the subdomain
across the edge: start4a_sum does not reach them, and this share of the parcel was lost. The cells along
such an edge (the lines where four subdomains meet) then had 1/16 too little solid volume, 6 % of the
bed (2 x 2 subdomains in x and y: 0.06 % of the solid volume of the grid, the same share of the
momentum and of the drag of the two-way coupling).

edge_fold moves the deposits of the edge and corner ghost cells into the ghost cells of a face of a
neighbour, from where start4a_sum takes them on:
  x: the ghost planes i = -1 and i = knox, their cells with j or k in the ghost layer, go to the
     x neighbour, into its cells i = 0 and i = knox-1 of the same j and k (ghost cells of its y or z faces)
  y: the ghost planes j = -1 and j = knoy, their cells with k in the ghost layer (i within the
     subdomain), go to the y neighbour, into its ghost cells of the z faces
The subdomains are a tensor product (one neighbour per face), so the faces of two neighbours match.
--------------------------------------------------------------------*/
void CPM::deposit_sum(lexer *p, ghostcell *pgc, field &f)
{
    pfold(p,f);
    
    edge_fold(p,pgc,f);
    
    pgc->start4a_sum(p,f,1);
}

void CPM::edge_fold(lexer *p, ghostcell *pgc, field &f)
{
    const int nx = p->knox, ny = p->knoy, nz = p->knoz;
    int sn[6], rn[6];
    double *sb[6], *rb[6];
    
    auto exchange = [&](int qa, int qb, int count)
    {
        for(int q=0;q<6;++q)
        {
            const int c = (q==qa || q==qb) ? count : 0;
            
            sn[q] = rn[q] = c;
            
            if(ef_send[q].size()<size_t(c)+1)
            {
                ef_send[q].resize(size_t(c)+1);
                ef_recv[q].resize(size_t(c)+1);
            }
            
            sb[q] = ef_send[q].data();
            rb[q] = ef_recv[q].data();
        }
    };
    
    // x: sides 1 (-x, buffer 0) and 4 (+x, buffer 3)
    if(xsplit==1 && (ysplit==1 || zsplit==1))
    {
        const int count = (ny+2)*(nz+2) - ny*nz;
        
        exchange(0,3,count);
        
        // the frame of the ghost plane ig, added into the plane ic of the neighbour
        auto frame = [&](int ig, int ic, double *buf, bool send)
        {
            int c=0;
            
            for(int jj=-1;jj<=ny;++jj)
            for(int kk=-1;kk<=nz;++kk)
            if(jj<0 || jj>=ny || kk<0 || kk>=nz)
            {
                if(send)
                {
                    buf[c] = f(ig,jj,kk);
                    f(ig,jj,kk) = 0.0;
                }
                else
                f(ic,jj,kk) += buf[c];
                
                ++c;
            }
        };
        
        frame(-1,0,sb[0],true);
        frame(nx,0,sb[3],true);
        
        for(int c=0;c<count;++c)
        rb[0][c] = rb[3][c] = 0.0;
        
        pgc->gcpartx(sn,rn,sb,rb);
        
        if(p->nb1>=0)
        frame(0,0,rb[0],false);
        
        if(p->nb4>=0)
        frame(0,nx-1,rb[3],false);
    }
    
    // y: sides 3 (-y, buffer 2) and 2 (+y, buffer 1)
    if(ysplit==1 && zsplit==1)
    {
        const int count = 2*nx;
        
        exchange(2,1,count);
        
        auto frame = [&](int jg, int jc, double *buf, bool send)
        {
            int c=0;
            
            for(int ii=0;ii<nx;++ii)
            for(int kk=-1;kk<=nz;kk+=nz+1)
            {
                if(send)
                {
                    buf[c] = f(ii,jg,kk);
                    f(ii,jg,kk) = 0.0;
                }
                else
                f(ii,jc,kk) += buf[c];
                
                ++c;
            }
        };
        
        frame(-1,0,sb[2],true);
        frame(ny,0,sb[1],true);
        
        for(int c=0;c<count;++c)
        rb[2][c] = rb[1][c] = 0.0;
        
        pgc->gcpartx(sn,rn,sb,rb);
        
        if(p->nb3>=0)
        frame(0,0,rb[2],false);
        
        if(p->nb2>=0)
        frame(0,ny-1,rb[1],false);
    }
}

// serial periodic: parcels whose new position PX1 lies beyond a periodic side are moved by the
// period, together with their old position PX0 (the RK2 stages combine both)
void CPM::periodic_wrap(lexer *p, double *PX0, double *PY0, double *PX1, double *PY1)
{
    double Lx = p->global_xmax - p->global_xmin;
    double Ly = p->global_ymax - p->global_ymin;
    
    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE)
    {
        if(perx==1)
        {
            if(PX1[n]>=p->global_xmax)
            {
                PX1[n] -= Lx;
                PX0[n] -= Lx;
            }
            else if(PX1[n]<p->global_xmin)
            {
                PX1[n] += Lx;
                PX0[n] += Lx;
            }
        }
        
        if(pery==1)
        {
            if(PY1[n]>=p->global_ymax)
            {
                PY1[n] -= Ly;
                PY0[n] -= Ly;
            }
            else if(PY1[n]<p->global_ymin)
            {
                PY1[n] += Ly;
                PY0[n] += Ly;
            }
        }
    }
}

// serial periodic: the ghost cells beyond a periodic side carry the flags of the cells on the
// opposite side, so the interpolation of the fluid velocities (which skips solid faces) uses the
// periodic ghost values next to the periodic side
void CPM::periodic_flags(lexer *p)
{
    auto id = [&](int ii, int jj, int kk) {return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;};
    int *fl[4] = {p->flag1,p->flag2,p->flag3,p->flag4};
    
    if(perx==1)
    for(int jj=p->jmin;jj<p->jmin+p->jmax;++jj)
    for(int kk=p->kmin;kk<p->kmin+p->kmax;++kk)
    for(int q=1;q<=-p->imin;++q)
    for(int f=0;f<4;++f)
    {
        fl[f][id(-q,jj,kk)] = fl[f][id(p->knox-q,jj,kk)];
        fl[f][id(p->knox-1+q,jj,kk)] = fl[f][id(q-1,jj,kk)];
    }
    
    if(pery==1)
    for(int ii=p->imin;ii<p->imin+p->imax;++ii)
    for(int kk=p->kmin;kk<p->kmin+p->kmax;++kk)
    for(int q=1;q<=-p->imin;++q)
    for(int f=0;f<4;++f)
    {
        fl[f][id(ii,-q,kk)] = fl[f][id(ii,p->knoy-q,kk)];
        fl[f][id(ii,p->knoy-1+q,kk)] = fl[f][id(ii,q-1,kk)];
    }
}
