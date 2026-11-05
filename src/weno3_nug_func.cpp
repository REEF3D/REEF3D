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

#include"weno3_nug_func.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"

weno3_nug_func::weno3_nug_func(lexer* p):epsilon(0.0),psi(1.0e-6)
{
    // the shared tables belong to the first lexer (the rank grid); an instance on another grid
    // (a mesh refinement patch) builds its own
    if(iniflag==1 && p!=s_lexer)
    own_ini(p);
    else
    ini(p);
}

weno3_nug_func::~weno3_nug_func()
{
    if(own_tables==1)
    {
    own_p->del_Darray(qfx,own_nx,2,4,2);
    own_p->del_Darray(qfy,own_ny,2,4,2);
    own_p->del_Darray(qfz,own_nz,2,4,2);

    own_p->del_Darray(cfx,own_nx,2,4);
    own_p->del_Darray(cfy,own_ny,2,4);
    own_p->del_Darray(cfz,own_nz,2,4);

    own_p->del_Darray(isfx,own_nx,2,4);
    own_p->del_Darray(isfy,own_ny,2,4);
    own_p->del_Darray(isfz,own_nz,2,4);
    }
}

void weno3_nug_func::own_ini(lexer* p)
{
    own_tables=1;
    own_p=p;
    own_nx=p->knox+8;
    own_ny=p->knoy+8;
    own_nz=p->knoz+8;

    p->Darray(qfx,p->knox+8,2,4,2);
    p->Darray(qfy,p->knoy+8,2,4,2);
    p->Darray(qfz,p->knoz+8,2,4,2);

    p->Darray(cfx,p->knox+8,2,4);
    p->Darray(cfy,p->knoy+8,2,4);
    p->Darray(cfz,p->knoz+8,2,4);

    p->Darray(isfx,p->knox+8,2,4);
    p->Darray(isfy,p->knoy+8,2,4);
    p->Darray(isfz,p->knoz+8,2,4);

    precalc_qf(p);
    precalc_cf(p);
    precalc_isf(p);
}

void weno3_nug_func::ini(lexer* p)
{
    if(iniflag==0)
    {
    p->Darray(s_qfx,p->knox+8,2,4,2);
    p->Darray(s_qfy,p->knoy+8,2,4,2);
    p->Darray(s_qfz,p->knoz+8,2,4,2);

    p->Darray(s_cfx,p->knox+8,2,4);
    p->Darray(s_cfy,p->knoy+8,2,4);
    p->Darray(s_cfz,p->knoz+8,2,4);

    p->Darray(s_isfx,p->knox+8,2,4);
    p->Darray(s_isfy,p->knoy+8,2,4);
    p->Darray(s_isfz,p->knoz+8,2,4);
    }

    qfx=s_qfx; qfy=s_qfy; qfz=s_qfz;
    cfx=s_cfx; cfy=s_cfy; cfz=s_cfz;
    isfx=s_isfx; isfy=s_isfy; isfz=s_isfz;

    if(iniflag==0)
    {
    precalc_qf(p);
    precalc_cf(p);
    precalc_isf(p);

    iniflag=1;
    s_lexer=p;
    }
}

double ****weno3_nug_func::s_qfx,****weno3_nug_func::s_qfy,****weno3_nug_func::s_qfz;
double ***weno3_nug_func::s_cfx,***weno3_nug_func::s_cfy,***weno3_nug_func::s_cfz;
double ***weno3_nug_func::s_isfx,***weno3_nug_func::s_isfy,***weno3_nug_func::s_isfz;
int weno3_nug_func::iniflag(0);
lexer *weno3_nug_func::s_lexer(nullptr);
