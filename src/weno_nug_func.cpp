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
    own_tables = true;
    own_nx = p->knox+2*marge;
    own_ny = p->knoy+2*marge;
    own_nz = p->knoz+2*marge;

    p->Darray(qfx,own_nx,2,6,2);
    p->Darray(qfy,own_ny,2,6,2);
    p->Darray(qfz,own_nz,2,6,2);

    p->Darray(cfx,own_nx,2,6);
    p->Darray(cfy,own_ny,2,6);
    p->Darray(cfz,own_nz,2,6);

    p->Darray(isfx,own_nx,2,6,3);
    p->Darray(isfy,own_ny,2,6,3);
    p->Darray(isfz,own_nz,2,6,3);

    precalc_qf(p);
    precalc_cf(p);
    precalc_isf(p);
}

weno_nug_func::~weno_nug_func()
{
    if(own_tables==1)
    {
        p->del_Darray(qfx,own_nx,2,6,2);
        p->del_Darray(qfy,own_ny,2,6,2);
        p->del_Darray(qfz,own_nz,2,6,2);

        p->del_Darray(cfx,own_nx,2,6);
        p->del_Darray(cfy,own_ny,2,6);
        p->del_Darray(cfz,own_nz,2,6);

        p->del_Darray(isfx,own_nx,2,6,3);
        p->del_Darray(isfy,own_ny,2,6,3);
        p->del_Darray(isfz,own_nz,2,6,3);
    }
}

void weno_nug_func::ini(lexer* p)
{
    if(!iniflag)
    {
        p->Darray(s_qfx,p->knox+2*marge,2,6,2);
        p->Darray(s_qfy,p->knoy+2*marge,2,6,2);
        p->Darray(s_qfz,p->knoz+2*marge,2,6,2);

        p->Darray(s_cfx,p->knox+2*marge,2,6);
        p->Darray(s_cfy,p->knoy+2*marge,2,6);
        p->Darray(s_cfz,p->knoz+2*marge,2,6);

        p->Darray(s_isfx,p->knox+2*marge,2,6,3);
        p->Darray(s_isfy,p->knoy+2*marge,2,6,3);
        p->Darray(s_isfz,p->knoz+2*marge,2,6,3);
    }

    qfx=s_qfx; qfy=s_qfy; qfz=s_qfz;
    cfx=s_cfx; cfy=s_cfy; cfz=s_cfz;
    isfx=s_isfx; isfy=s_isfy; isfz=s_isfz;

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

double ****weno_nug_func::s_qfx,****weno_nug_func::s_qfy,****weno_nug_func::s_qfz;
double ***weno_nug_func::s_cfx,***weno_nug_func::s_cfy,***weno_nug_func::s_cfz;
double ****weno_nug_func::s_isfx,****weno_nug_func::s_isfy,****weno_nug_func::s_isfz;
