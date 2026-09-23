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

#include "ddweno_f_nug.h"
#include "lexer.h"
#include "field.h"
#include "slice.h"
#include "ghostcell.h"

ddweno_f_nug::ddweno_f_nug(lexer* pp):weno_nug_func(pp)
{
    p=pp;

    // only the uniform-flux coefficient set is used
    uf = vf = wf = 0;
}

ddweno_f_nug::~ddweno_f_nug()
{
}

double ddweno_f_nug::ddwenox(field &f, double uw)
{
    if(uw>=0.0)
    {
        iqmin(f);
        return weno_min_x();
    }
    else
    {
        iqmax(f);
        return weno_max_x();
    }
}

double ddweno_f_nug::ddwenoy(field &f, double uw)
{
    if(uw>=0.0)
    {
        jqmin(f);
        return weno_min_y();
    }
    else
    {
        jqmax(f);
        return weno_max_y();
    }
}

double ddweno_f_nug::ddwenoz(field &f, double uw)
{
    if(uw>=0.0)
    {
        kqmin(f);
        return weno_min_z();
    }
    else
    {
        kqmax(f);
        return weno_max_z();
    }
}

double ddweno_f_nug::dswenox(slice &f, double uw)
{
    if(uw>=0.0)
    {
        isqmin(f);
        return weno_min_x();
    }
    else
    {
        isqmax(f);
        return weno_max_x();
    }
}

double ddweno_f_nug::dswenoy(slice &f, double uw)
{
    if(uw>=0.0)
    {
        jsqmin(f);
        return weno_min_y();
    }
    else
    {
        jsqmax(f);
        return weno_max_y();
    }
}

void ddweno_f_nug::dsdiffx(slice &f, slice &dq)
{
    // faces i-3 .. i+2 of the cells 0 .. knox-1
    for(int ii=-3; ii<p->knox+2; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    dq(ii,jj) = (f(ii+1,jj)-f(ii,jj))/p->DXP[ii+marge];
}

void ddweno_f_nug::dsdiffy(slice &f, slice &dq)
{
    // faces j-3 .. j+2 of the cells 0 .. knoy-1
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=-3; jj<p->knoy+2; ++jj)
    dq(ii,jj) = (f(ii,jj+1)-f(ii,jj))/p->DYP[jj+marge];
}
