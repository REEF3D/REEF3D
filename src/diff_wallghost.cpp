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

#include"diff_wallghost.h"
#include"lexer.h"
#include"ghostcell.h"
#include<cmath>

diff_wallghost::diff_wallghost(lexer *p) : P(p), gcb(nullptr), gcb_count(0)
{
}

diff_wallghost::~diff_wallghost()
{
}

// index of the first ghost cell of face cs of cell (ii,jj,kk); -1 if outside the local arrays
int diff_wallghost::ghost_index(lexer *p, int ii, int jj, int kk, int cs)
{
    if(cs==1) --ii;
    if(cs==4) ++ii;
    if(cs==3) --jj;
    if(cs==2) ++jj;
    if(cs==5) --kk;
    if(cs==6) ++kk;

    if(ii<p->imin || ii>=p->imin+p->imax || jj<p->jmin || jj>=p->jmin+p->jmax || kk<p->kmin || kk>=p->kmin+p->kmax)
    return -1;

    return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin;
}

// fill P on the cells of the component (flag>0) with 0, 1, m or m^2 (m = global i+j+k), 0 elsewhere,
// run the ghost-cell routines and read the first ghost cell of every boundary face
void diff_wallghost::probe(lexer *p, ghostcell *pgc, int c, int *flag, int mode, std::vector<double> &res)
{
    for(i=p->imin; i<p->imin+p->imax; ++i)
    for(j=p->jmin; j<p->jmin+p->jmax; ++j)
    for(k=p->kmin; k<p->kmin+p->kmax; ++k)
    {
        int q = (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin;
        double m = double(i+p->origin_i + j+p->origin_j + k+p->origin_k);
        double val = 0.0;

        if(flag[q]>0)
        {
        if(mode==1)
        val = 1.0;

        if(mode==2)
        val = m;

        if(mode==3)
        val = m*m;
        }

        P.V[q] = val;
    }

    if(c==0)
    pgc->start1(p,P,gcv_probe);

    if(c==1)
    pgc->start2(p,P,gcv_probe);

    if(c==2)
    pgc->start3(p,P,gcv_probe);

    res.assign(gcb_count,0.0);

    for(int q=0; q<gcb_count; ++q)
    {
    int g = ghost_index(p,gcb[q][0],gcb[q][1],gcb[q][2],gcb[q][3]);

        if(g>=0)
        res[q] = P.V[g];
    }
}

void diff_wallghost::update(lexer *p, ghostcell *pgc, int c, int gcv)
{
    int *flag;

    if(c==0)
    {
    gcb = p->gcb1;
    gcb_count = p->gcb1_count;
    flag = p->flag1;
    }

    if(c==1)
    {
    gcb = p->gcb2;
    gcb_count = p->gcb2_count;
    flag = p->flag2;
    }

    if(c==2)
    {
    gcb = p->gcb3;
    gcb_count = p->gcb3_count;
    flag = p->flag3;
    }

    gcv_probe = gcv;

    // which face writes which ghost cell
    owner.assign(size_t(p->imax)*size_t(p->jmax)*size_t(p->kmax),-1);

    for(int q=0; q<gcb_count; ++q)
    {
    int cs = gcb[q][3];

        if(cs<1 || cs>6 || ((cs==2 || cs==3) && p->j_dir==0))
        continue;

    int g = ghost_index(p,gcb[q][0],gcb[q][1],gcb[q][2],cs);

        if(g<0)
        continue;

        owner[g] = (owner[g]==-1) ? q : -2;
    }

    std::vector<double> z, a, b, d;

    probe(p,pgc,c,flag,0,z);
    probe(p,pgc,c,flag,1,a);
    probe(p,pgc,c,flag,2,b);
    probe(p,pgc,c,flag,3,d);

    w1.assign(gcb_count,0.0);
    w2.assign(gcb_count,0.0);
    valid.assign(gcb_count,0);

    for(int q=0; q<gcb_count; ++q)
    {
    int cs = gcb[q][3];

        if(cs<1 || cs>6)
        continue;

    int g = ghost_index(p,gcb[q][0],gcb[q][1],gcb[q][2],cs);

        if(g<0 || owner[g]!=q)
        continue;

    // the next interior cell is one step further away from the wall
    const double s = (cs==1 || cs==3 || cs==5) ? 1.0 : -1.0;
    const double m = double(gcb[q][0]+p->origin_i + gcb[q][1]+p->origin_j + gcb[q][2]+p->origin_k);

    const double A = a[q]-z[q];
    const double B = b[q]-z[q];
    const double C = d[q]-z[q];

    const double c2 = s*(B - A*m);
    const double c1 = A - c2;

    // homogeneous, and only the two cells on the normal line
    const double check = c1*m*m + c2*(m+s)*(m+s);

        if(fabs(z[q])<1.0e-12 && fabs(C-check)<=1.0e-9*(1.0+m*m) && fabs(c1)<10.0 && fabs(c2)<10.0)
        {
        w1[q] = c1;
        w2[q] = c2;
        valid[q] = 1;
        }

    }
}

bool diff_wallghost::coef(lexer *p, int ii, int jj, int kk, int cs, double &c1, double &c2)
{
    int g = ghost_index(p,ii,jj,kk,cs);

    if(g<0)
    return false;

    int q = owner[g];

    if(q<0 || !valid[q] || gcb[q][0]!=ii || gcb[q][1]!=jj || gcb[q][2]!=kk || gcb[q][3]!=cs)
    return false;

    c1 = w1[q];
    c2 = w2[q];

    return true;
}
