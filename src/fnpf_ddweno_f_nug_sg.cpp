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

#include "fnpf_ddweno_f_nug.h"
#include "lexer.h"
#include "field.h"
#include "slice.h"

fnpf_ddweno_f_nug::fnpf_ddweno_f_nug(lexer* pp) : weno_nug_func(pp)
{
    p=pp;

    // only the uniform-flux coefficient set is used
    uf = vf = wf = 0;
}

fnpf_ddweno_f_nug::~fnpf_ddweno_f_nug()
{
}

double fnpf_ddweno_f_nug::dswenox(slice& f, double uw)
{
    if(uw>0.0)
    {
        isqmin(f);
        return weno_min_x();
    }
    else if(uw<0.0)
    {
        isqmax(f);
        return weno_max_x();
    }
    else
        return 0.0;
}

double fnpf_ddweno_f_nug::dswenoy(slice& f, double uw)
{
    if(uw>0.0)
    {
        jsqmin(f);
        return weno_min_y();
    }
    else if(uw<0.0)
    {
        jsqmax(f);
        return weno_max_y();
    }
    else
        return 0.0;
}

void fnpf_ddweno_f_nug::dsdiffx(slice &f, slice &dq)
{
    // faces i-3 .. i+2 of the cells 0 .. knox-1
    for(int ii=-3; ii<p->knox+2; ++ii)
    for(int jj=0; jj<p->knoy; ++jj)
    dq(ii,jj) = (f(ii+1,jj)-f(ii,jj))/p->DXP[ii+marge];
}

void fnpf_ddweno_f_nug::dsdiffy(slice &f, slice &dq)
{
    // faces j-3 .. j+2 of the cells 0 .. knoy-1
    for(int ii=0; ii<p->knox; ++ii)
    for(int jj=-3; jj<p->knoy+2; ++jj)
    dq(ii,jj) = (f(ii,jj+1)-f(ii,jj))/p->DYP[jj+marge];
}

void fnpf_ddweno_f_nug::isqmin(slice& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->wet[Im2J]>0 && p->wet[Im3J]>0)
    if(p->wet[Im1J]>0 && p->wet[Im2J]>0)
    if(p->wet[Im1J]>0 && p->wet[IJ]>0)
    if(p->wet[Ip1J]>0 && p->wet[IJ]>0)
    if(p->wet[Ip2J]>0 && p->wet[Ip1J]>0)
    {
        if(p->wet[Im2J]>0 && p->wet[Im3J]>0)
        q1 = (f(i-2,j)-f(i-3,j))/p->DXP[IM3];

        if(p->wet[Im1J]>0 && p->wet[Im2J]>0)
        q2 = (f(i-1,j)-f(i-2,j))/p->DXP[IM2];

        if(p->wet[Im1J]>0 && p->wet[IJ]>0)
        q3 = (f(i,j)-f(i-1,j))/p->DXP[IM1];

        if(p->wet[Ip1J]>0 && p->wet[IJ]>0)
        q4 = (f(i+1,j)-f(i,j))/p->DXP[IP];

        if(p->wet[Ip2J]>0 && p->wet[Ip1J]>0)
        q5 = (f(i+2,j)-f(i+1,j))/p->DXP[IP1];
    }
}

void fnpf_ddweno_f_nug::jsqmin(slice& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->wet[IJm2]>0 && p->wet[IJm3]>0)
    if(p->wet[IJm1]>0 && p->wet[IJm2]>0)
    if(p->wet[IJ]>0 && p->wet[IJm1]>0)
    if(p->wet[IJp1]>0 && p->wet[IJ]>0)
    if(p->wet[IJp2]>0 && p->wet[IJp1]>0)
    {
        if(p->wet[IJm2]>0 && p->wet[IJm3]>0)
        q1 = (f(i,j-2)-f(i,j-3))/p->DYP[JM3];

        if(p->wet[IJm1]>0 && p->wet[IJm2]>0)
        q2 = (f(i,j-1)-f(i,j-2))/p->DYP[JM2];

        if(p->wet[IJ]>0 && p->wet[IJm1]>0)
        q3 = (f(i,j)-f(i,j-1))/p->DYP[JM1];

        if(p->wet[IJp1]>0 && p->wet[IJ]>0)
        q4 = (f(i,j+1)-f(i,j))/p->DYP[JP];

        if(p->wet[IJp2]>0 && p->wet[IJp1]>0)
        q5 = (f(i,j+2)-f(i,j+1))/p->DYP[JP1];
    }
}

void fnpf_ddweno_f_nug::isqmax(slice& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->wet[Im1J]>0 && p->wet[Im2J]>0)
    if(p->wet[IJ]>0 && p->wet[Im1J]>0)
    if(p->wet[Ip1J]>0 && p->wet[IJ]>0)
    if(p->wet[Ip2J]>0 && p->wet[Ip1J]>0)
    if(p->wet[Ip3J]>0 && p->wet[Ip2J]>0)
    {
        if(p->wet[Im1J]>0 && p->wet[Im2J]>0)
        q1 = (f(i-1,j)-f(i-2,j))/p->DXP[IM2];

        if(p->wet[IJ]>0 && p->wet[Im1J]>0)
        q2 = (f(i,j)-f(i-1,j))/p->DXP[IM1];

        if(p->wet[Ip1J]>0 && p->wet[IJ]>0)
        q3 = (f(i+1,j)-f(i,j))/p->DXP[IP];

        if(p->wet[Ip2J]>0 && p->wet[Ip1J]>0)
        q4 = (f(i+2,j)-f(i+1,j))/p->DXP[IP1];

        if(p->wet[Ip3J]>0 && p->wet[Ip2J]>0)
        q5 = (f(i+3,j)-f(i+2,j))/p->DXP[IP2];
    }
}

void fnpf_ddweno_f_nug::jsqmax(slice& f)
{
    q1=q2=q3=q4=q5=0.0;

    if(p->wet[IJm1]>0 && p->wet[IJm2]>0)
    if(p->wet[IJ]>0 && p->wet[IJm1]>0)
    if(p->wet[IJp1]>0 && p->wet[IJ]>0)
    if(p->wet[IJp2]>0 && p->wet[IJp1]>0)
    if(p->wet[IJp3]>0 && p->wet[IJp2]>0)
    {
        if(p->wet[IJm1]>0 && p->wet[IJm2]>0)
        q1 = (f(i,j-1)-f(i,j-2))/p->DYP[JM2];

        if(p->wet[IJ]>0 && p->wet[IJm1]>0)
        q2 = (f(i,j)-f(i,j-1))/p->DYP[JM1];

        if(p->wet[IJp1]>0 && p->wet[IJ]>0)
        q3 = (f(i,j+1)-f(i,j))/p->DYP[JP];

        if(p->wet[IJp2]>0 && p->wet[IJp1]>0)
        q4 = (f(i,j+2)-f(i,j+1))/p->DYP[JP1];

        if(p->wet[IJp3]>0 && p->wet[IJp2]>0)
        q5 = (f(i,j+3)-f(i,j+2))/p->DYP[JP2];
    }
}
