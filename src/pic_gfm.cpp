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

#include"pic_gfm.h"
#include"lexer.h"
#include"fdm.h"
#include"field.h"
#include<algorithm>
#include<map>

namespace
{
inline int gid(lexer *p, int i, int j, int k)
{
    return (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin;
}

inline double theta(double fin, double fout)
{
    return std::max(PIC_GFM_THETA, fin/(fin-fout));
}

const std::function<bool(lexer*,int,int,int,int)> *gfm_keep = nullptr;

std::map<lexer*,field*> gfm_air;
}

void pic_gfm_air_field(lexer *p, field *f)
{
    if(p!=nullptr)
    {
        gfm_air[p] = f;
        return;
    }
    for(auto it=gfm_air.begin(); it!=gfm_air.end();)
    if(it->second==f)
    it = gfm_air.erase(it);
    else
    ++it;
}

double pic_gfm_air(lexer *p, int i, int j, int k)
{
    if(gfm_air.empty())
    return 0.0;
    const auto it = gfm_air.find(p);
    if(it==gfm_air.end())
    return 0.0;
    return (*it->second)(i,j,k);
}

double pic_gfm_theta(double fin, double fout)
{
    return theta(fin,fout);
}

void pic_gfm_keep(const std::function<bool(lexer*,int,int,int,int)> *keep)
{
    gfm_keep = keep;
}

bool gfm_single(lexer *p)
{
    return p->F81==1 || ((p->F41==1 || p->F41==2) && p->F30>0 && p->F81==0);
}

void pic_gfm_matrix(lexer *p, fdm *a)
{
    int n=0;

    for(int i=0; i<p->knox; ++i)
    for(int j=0; j<p->knoy; ++j)
    for(int k=0; k<p->knoz; ++k)
    {
        if(p->flag4[gid(p,i,j,k)]<=0)    // the cells of LOOP (PCHECK)
        continue;

        const double f0 = a->phi(i,j,k);

        // a cell in the air: p = 0
        if(f0<0.0)
        {
        a->M.p[n] = 1.0;
        a->M.n[n] = a->M.s[n] = a->M.e[n] = a->M.w[n] = a->M.t[n] = a->M.b[n] = 0.0;
        a->rhsvec.V[n] = 0.0;
        ++n;
        continue;
        }

        // the faces of the water to the air: p = 0 at the surface between
        double *coef[6] = {&a->M.n[n],&a->M.s[n],&a->M.w[n],&a->M.e[n],&a->M.t[n],&a->M.b[n]};
        const int di[6] = {1,-1,0,0,0,0};
        const int dj[6] = {0,0,1,-1,0,0};
        const int dk[6] = {0,0,0,0,1,-1};

        for(int s=0; s<6; ++s)
        {
            if((s==2 || s==3) && p->j_dir==0)
            continue;

            const int ni = i+di[s], nj = j+dj[s], nk = k+dk[s];

            if(p->flag4[gid(p,ni,nj,nk)]<0 || *coef[s]==0.0)
            continue;

            const double f1 = a->phi(ni,nj,nk);

            if(f1>=0.0)
            continue;

            // a face at a patch edge: the composite operator (cfd_amr::pr_matrix, pr_apply)
            if(gfm_keep!=nullptr && (*gfm_keep)(p,i,j,k,s))
            continue;

            const double c = -(*coef[s]);
            a->M.p[n] += c*(1.0/theta(f0,f1) - 1.0);
            *coef[s] = 0.0;
        }

        ++n;
    }
}

double pic_gfm_grad(lexer *p, fdm *a, field &pr, int i, int j, int k, int di, int dj, int dk, double dx)
{
    const double f0 = a->phi(i,j,k);
    const double f1 = a->phi(i+di,j+dj,k+dk);

    if(f0>=0.0 && f1>=0.0)
    return (pr(i+di,j+dj,k+dk) - pr(i,j,k))/dx;

    // the pressure of the surface: 0, or that of the air pocket (F 41 2) for the pressure itself
    const bool air = &pr==&a->press;

    if(f0>=0.0)
    return ((air ? pic_gfm_air(p,i+di,j+dj,k+dk) : 0.0) - pr(i,j,k))/(theta(f0,f1)*dx);

    if(f1>=0.0)
    return (pr(i+di,j+dj,k+dk) - (air ? pic_gfm_air(p,i,j,k) : 0.0))/(theta(f1,f0)*dx);

    return 0.0;
}
