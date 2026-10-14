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

#include"sediment_xfilter.h"
#include"sediment_curv.h"
#include"lexer.h"
#include"ghostcell.h"

sediment_xfilter::sediment_xfilter(lexer *p) : own(p), m(p), h(p), dh(p), a(p), b(p), wxm(p), wxp(p), wym(p), wyp(p)
{
    beta = 0.25*(1.0-MAX(0.0,MIN(1.0,p->S102)));
    ydir = (p->j_dir==1 && p->gknoy>1)?1:0;
}

void sediment_xfilter::mask(lexer *p, slice *amr)
{
    SLICEBASELOOP
    own(i,j) = 0.0;

    SEDSLICELOOP
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    own(i,j) = 1.0;

    if(amr!=nullptr)
    SLICEBASELOOP
    if((*amr)(i,j)<0.5)
    own(i,j) = 0.0;

    SLICEBASELOOP
    m(i,j) = own(i,j);
}

void sediment_xfilter::weights(lexer *p)
{
    auto area = [&](int ii, int jj)
    {
        if(p->geo_curv==1)
        return sediment_curv::area(p,ii,jj);
        return p->DXN[ii+marge]*p->DYN[jj+marge];
    };

    SLICEBASELOOP
    {
    wxm(i,j)=wxp(i,j)=wym(i,j)=wyp(i,j)=0.0;

        if(own(i,j)>0.5)
        {
        const double Ai = area(i,j);

        if(i+p->origin_i>0 && m(i-1,j)>0.5)
        wxm(i,j) = beta*MIN(Ai,area(i-1,j))/Ai;

        if(i+p->origin_i<p->gknox-1 && m(i+1,j)>0.5)
        wxp(i,j) = beta*MIN(Ai,area(i+1,j))/Ai;

            if(ydir==1)
            {
            if(j+p->origin_j>0 && m(i,j-1)>0.5)
            wym(i,j) = beta*MIN(Ai,area(i,j-1))/Ai;

            if(j+p->origin_j<p->gknoy-1 && m(i,j+1)>0.5)
            wyp(i,j) = beta*MIN(Ai,area(i,j+1))/Ai;
            }
        }
    }
}

void sediment_xfilter::copy(lexer *p, slice &f, int what)
{
    if(what==0)
    SLICEBASELOOP
    h(i,j) = f(i,j);

    if(what==1)
    SLICEBASELOOP
    dh(i,j) = h(i,j) - f(i,j);
}

void sediment_xfilter::sweep1(lexer *p, int src)
{
    slice4 &in = (src==1) ? h : dh;

    SLICEBASELOOP
    {
    a(i,j) = in(i,j) + wxm(i,j)*(in(i-1,j)-in(i,j)) + wxp(i,j)*(in(i+1,j)-in(i,j));
    b(i,j) = in(i,j) + wym(i,j)*(in(i,j-1)-in(i,j)) + wyp(i,j)*(in(i,j+1)-in(i,j));
    }
}

void sediment_xfilter::sweep2(lexer *p, slice &f, int set)
{
    double out;

    SLICEBASELOOP
    if(own(i,j)>0.5)
    {
        if(ydir==1)
        out = 0.5*(a(i,j) + wym(i,j)*(a(i,j-1)-a(i,j)) + wyp(i,j)*(a(i,j+1)-a(i,j))
                 + b(i,j) + wxm(i,j)*(b(i-1,j)-b(i,j)) + wxp(i,j)*(b(i+1,j)-b(i,j)));
        else
        out = a(i,j);

        if(set==1)
        f(i,j) = out;
        else
        f(i,j) += out;
    }
}

void sediment_xfilter::run(lexer *p, ghostcell *pgc, slice &f, int outer_iter, int inner_iter)
{
    mask(p,nullptr);
    pgc->gcsl_start4(p,m,1);
    weights(p);

    for(int qn=0; qn<outer_iter; ++qn)
    {
        copy(p,f,0);
        pgc->gcsl_start4(p,h,1);
        sweep1(p,1);
        pgc->gcsl_start4(p,a,1);
        pgc->gcsl_start4(p,b,1);
        sweep2(p,f,1);

        for(int qqn=0; qqn<inner_iter; ++qqn)
        {
            copy(p,f,1);
            pgc->gcsl_start4(p,dh,1);
            sweep1(p,2);
            pgc->gcsl_start4(p,a,1);
            pgc->gcsl_start4(p,b,1);
            sweep2(p,f,0);
        }
    }
}
