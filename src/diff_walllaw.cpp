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

#include"diff_walllaw.h"
#include"lexer.h"
#include"bc_noflux.h"
#include"matrix_diag.h"

diff_walllaw::diff_walllaw() : active(false)
{
}

long long diff_walllaw::key(lexer *p, int ii, int jj, int kk, int cs) const
{
    return ((long long)((ii-p->imin)*p->jmax + (jj-p->jmin))*p->kmax + (kk-p->kmin))*8 + cs;
}

static inline bool walllaw_type(int bc)
{
    return bc==5 || bc==21 || bc==22 || bc==41 || bc==42 || bc==43;
}

void diff_walllaw::update(lexer *p, int c)
{
    S.clear();
    active = p->B10!=0 && (c!=1 || p->j_dir==1);

    if(!active)
    return;

    int **gcb = c==0 ? p->gcb1 : (c==1 ? p->gcb2 : p->gcb3);
    const int gcb_count = c==0 ? p->gcb1_count : (c==1 ? p->gcb2_count : p->gcb3_count);
    int **gcdf = c==0 ? p->gcdf1 : (c==1 ? p->gcdf2 : p->gcdf3);
    const int gcdf_count = c==0 ? p->gcdf1_count : (c==1 ? p->gcdf2_count : p->gcdf3_count);

    // the faces normal to the component carry no wall law (bcmom: tangential faces only)
    const int n1 = c==0 ? 1 : (c==1 ? 2 : 5);
    const int n2 = c==0 ? 4 : (c==1 ? 3 : 6);

    for(int q=0; q<gcb_count; ++q)
    if(gcb[q][3]>0 && walllaw_type(gcb[q][4]) && gcb[q][3]!=n1 && gcb[q][3]!=n2 && !bc_periodic_face(p,gcb[q][0],gcb[q][1],gcb[q][2],gcb[q][3]))
    S.insert(key(p,gcb[q][0],gcb[q][1],gcb[q][2],gcb[q][3]));

    for(int q=0; q<gcdf_count; ++q)
    if(gcdf[q][3]>0 && gcdf[q][3]!=n1 && gcdf[q][3]!=n2)
    S.insert(key(p,gcdf[q][0],gcdf[q][1],gcdf[q][2],gcdf[q][3]));
}

bool diff_walllaw::face(lexer *p, int ii, int jj, int kk, int cs) const
{
    if(!active)
    return false;

    return S.count(key(p,ii,jj,kk,cs))>0;
}

// zero gradient through the wall-law faces in the rows of the implicit diffusion matrix of component c
// (rows in the order of ULOOP / VLOOP / WLOOP, as assembled): M.p += M.x, M.x = 0
void diff_walllaw::rows(lexer *p, matrix_diag &M, int c)
{
    if(!active)
    return;

    int n=0;

    auto apply = [&]()
    {
        for(int cs=1; cs<=6; ++cs)
        {
            if((cs==2 || cs==3) && p->j_dir==0)
            continue;

            if(!face(p,i,j,k,cs))
            continue;

            double &mx = cs==1 ? M.s[n] : cs==4 ? M.n[n] : cs==3 ? M.e[n] : cs==2 ? M.w[n] : cs==5 ? M.b[n] : M.t[n];
            M.p[n] += mx;
            mx = 0.0;
        }
        ++n;
    };

    if(c==0)
    ULOOP
    apply();

    if(c==1)
    VLOOP
    apply();

    if(c==2)
    WLOOP
    apply();
}
