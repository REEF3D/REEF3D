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

#include"weno_nug_func.h"
#include"lexer.h"

weno_nug_func::weno_nug_func(lexer *p)
{
    // the shared tables belong to the first lexer (the rank grid); an instance on another
    // grid (a mesh refinement patch) builds its own
    if(iniflag==1 && p!=s_lexer)
    own_ini(p);
    else
    ini(p);

    weno_nug_func::p = p;
}

weno_nug_func::weno_nug_func(lexer *p, int)
{
    own_ini(p);

    weno_nug_func::p = p;
}

void weno_nug_func::own_ini(lexer* p)
{
    own_qfx.resize(p->knox+2*marge);
    own_qfy.resize(p->knoy+2*marge);
    own_qfz.resize(p->knoz+2*marge);

    own_cfx.resize(p->knox+2*marge);
    own_cfy.resize(p->knoy+2*marge);
    own_cfz.resize(p->knoz+2*marge);

    own_isfx.resize(p->knox+2*marge);
    own_isfy.resize(p->knoy+2*marge);
    own_isfz.resize(p->knoz+2*marge);

    qfx=own_qfx.data(); qfy=own_qfy.data(); qfz=own_qfz.data();
    cfx=own_cfx.data(); cfy=own_cfy.data(); cfz=own_cfz.data();
    isfx=own_isfx.data(); isfy=own_isfy.data(); isfz=own_isfz.data();

    precalc_qf(p);
    precalc_cf(p);
    precalc_isf(p);
}

void weno_nug_func::ini(lexer* p)
{
    if(!iniflag)
    {
        s_qfx.resize(p->knox+2*marge);
        s_qfy.resize(p->knoy+2*marge);
        s_qfz.resize(p->knoz+2*marge);

        s_cfx.resize(p->knox+2*marge);
        s_cfy.resize(p->knoy+2*marge);
        s_cfz.resize(p->knoz+2*marge);

        s_isfx.resize(p->knox+2*marge);
        s_isfy.resize(p->knoy+2*marge);
        s_isfz.resize(p->knoz+2*marge);
    }

    qfx=s_qfx.data(); qfy=s_qfy.data(); qfz=s_qfz.data();
    cfx=s_cfx.data(); cfy=s_cfy.data(); cfz=s_cfz.data();
    isfx=s_isfx.data(); isfy=s_isfy.data(); isfz=s_isfz.data();

    if(!iniflag)
    {
        precalc_qf(p);
        precalc_cf(p);
        precalc_isf(p);

        iniflag = true;
        s_lexer = p;
    }
}

void weno_nug_func::dsdiffx(slice &f, slice &dq)
{
    // faces i-3 .. i+2 of the cells 0 .. knox-1
    for(int ii=-3; ii<p->knox+2; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    dq(ii,jj) = (f(ii+1,jj)-f(ii,jj))/p->DXP[ii+marge];
}

void weno_nug_func::dsdiffy(slice &f, slice &dq)
{
    // faces j-3 .. j+2 of the cells 0 .. knoy-1
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=-3; jj<p->knoy+2; ++jj)
    dq(ii,jj) = (f(ii,jj+1)-f(ii,jj))/p->DYP[jj+marge];
}
