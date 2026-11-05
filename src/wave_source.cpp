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

#include"wave_source.h"
#include"wave_lib.h"
#include"lexer.h"
#include<cmath>

void wave_lexer_context::save(const lexer *p)
{
    wT=p->wT; wC=p->wC; wA=p->wA; wk=p->wk; wL=p->wL; wwp=p->wwp; ww=p->ww; wd=p->wd;
    wTp=p->wTp; wHs=p->wHs; wH=p->wH; wAs=p->wAs; phiin=p->phiin; wts=p->wts; wte=p->wte;
    wLp=p->wLp; ww_s=p->ww_s; ww_e=p->ww_e;
    wN=p->wN;

    B84=p->B84; B85=p->B85; B86=p->B86; B87=p->B87; B91=p->B91; B92=p->B92; B93=p->B93;
    B94=p->B94; B130=p->B130; B133=p->B133; B136=p->B136; B138=p->B138; B138_1=p->B138_1;
    B138_2=p->B138_2; B139=p->B139;
    B87_1=p->B87_1; B87_2=p->B87_2; B88=p->B88; B91_1=p->B91_1; B91_2=p->B91_2;
    B93_1=p->B93_1; B93_2=p->B93_2; B94_wdt=p->B94_wdt; B105_1=p->B105_1; B131=p->B131;
    B132_s=p->B132_s; B132_e=p->B132_e; B134=p->B134;
}

void wave_lexer_context::load(lexer *p) const
{
    p->wT=wT; p->wC=wC; p->wA=wA; p->wk=wk; p->wL=wL; p->wwp=wwp; p->ww=ww; p->wd=wd;
    p->wTp=wTp; p->wHs=wHs; p->wH=wH; p->wAs=wAs; p->phiin=phiin; p->wts=wts; p->wte=wte;
    p->wLp=wLp; p->ww_s=ww_s; p->ww_e=ww_e;
    p->wN=wN;

    p->B84=B84; p->B85=B85; p->B86=B86; p->B87=B87; p->B91=B91; p->B92=B92; p->B93=B93;
    p->B94=B94; p->B130=B130; p->B133=B133; p->B136=B136; p->B138=B138; p->B138_1=B138_1;
    p->B138_2=B138_2; p->B139=B139;
    p->B87_1=B87_1; p->B87_2=B87_2; p->B88=B88; p->B91_1=B91_1; p->B91_2=B91_2;
    p->B93_1=B93_1; p->B93_2=B93_2; p->B94_wdt=B94_wdt; p->B105_1=B105_1; p->B131=B131;
    p->B132_s=B132_s; p->B132_e=B132_e; p->B134=B134;
}

wave_source::wave_source(int iid, int itype) : id(iid), type(itype), H(0.0), T(0.0), rot(0.0), phase(0.0),
                                               ts(0.0), te(1.0e20), t_ramp(0.0), x0(0.0), y0(0.0), seed(0),
                                               tshift(0.0), cr(1.0), sr(0.0), lib(nullptr)
{
}

wave_source::~wave_source()
{
    delete lib;
}

bool wave_source::active(const lexer *p) const
{
    return p->simtime>=ts && p->simtime<=te && p->simtime>=ctx.wts && p->simtime<=ctx.wte;
}

double wave_source::ramp(const lexer *p) const
{
    if(t_ramp<=0.0)
    return 1.0;

    const double pi = 3.14159265358979323846;
    double r = 1.0;

    if(p->simtime < ts+t_ramp)
    r = 0.5*(1.0 - cos(pi*(p->simtime-ts)/t_ramp));

    if(te<1.0e19 && p->simtime > te-t_ramp)
    r = fmin(r, 0.5*(1.0 - cos(pi*(te-p->simtime)/t_ramp)));

    return fmax(0.0,fmin(1.0,r));
}

void wave_source::local(double x, double y, double &xs, double &ys) const
{
    xs =  (x-x0)*cr + (y-y0)*sr;
    ys = -(x-x0)*sr + (y-y0)*cr;
}
