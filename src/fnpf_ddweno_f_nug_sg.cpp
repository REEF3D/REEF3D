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
}

fnpf_ddweno_f_nug::~fnpf_ddweno_f_nug()
{
}

double fnpf_ddweno_f_nug::dswenox(slice& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf = 0;

    if(uw==0.0)
    {
        return 0.0;
    }
    else if(uw>0.0)
    {
        isqmin(f);
        is_min_x();
        weight_min_x();

        return w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
             + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
             + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
    }
    else if(uw<0.0)
    {
        isqmax(f);
        is_max_x();
        weight_max_x();

        return w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
             + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
             + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
    }
}

double fnpf_ddweno_f_nug::dswenoy(slice& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf = 0;

    if(uw==0.0)
    {
        return 0.0;
    }
    else if(uw>0.0)
    {
        jsqmin(f);
        is_min_y();
        weight_min_y();

        return w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
             + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
             + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
    }
    else if(uw<0.0)
    {
        jsqmax(f);
        is_max_y();
        weight_max_y();

        return w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
             + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
             + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
    }
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
        q1 = (f(i-2,j)-f(i-3,j))/DX[IM3];

        if(p->wet[Im1J]>0 && p->wet[Im2J]>0)
        q2 = (f(i-1,j)-f(i-2,j))/DX[IM2];

        if(p->wet[Im1J]>0 && p->wet[IJ]>0)
        q3 = (f(i,j)-f(i-1,j))/DX[IM1];

        if(p->wet[Ip1J]>0 && p->wet[IJ]>0)
        q4 = (f(i+1,j)-f(i,j))/DX[IP];

        if(p->wet[Ip2J]>0 && p->wet[Ip1J]>0)
        q5 = (f(i+2,j)-f(i+1,j))/DX[IP1];
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
        q1 = (f(i,j-2)-f(i,j-3))/DY[JM3];

        if(p->wet[IJm1]>0 && p->wet[IJm2]>0)
        q2 = (f(i,j-1)-f(i,j-2))/DY[JM2];

        if(p->wet[IJ]>0 && p->wet[IJm1]>0)
        q3 = (f(i,j)-f(i,j-1))/DY[JM1];

        if(p->wet[IJp1]>0 && p->wet[IJ]>0)
        q4 = (f(i,j+1)-f(i,j))/DY[JP];

        if(p->wet[IJp2]>0 && p->wet[IJp1]>0)
        q5 = (f(i,j+2)-f(i,j+1))/DY[JP1];
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
        q1 = (f(i-1,j)-f(i-2,j))/DX[IM2];

        if(p->wet[IJ]>0 && p->wet[Im1J]>0)
        q2 = (f(i,j)-f(i-1,j))/DX[IM1];

        if(p->wet[Ip1J]>0 && p->wet[IJ]>0)
        q3 = (f(i+1,j)-f(i,j))/DX[IP];

        if(p->wet[Ip2J]>0 && p->wet[Ip1J]>0)
        q4 = (f(i+2,j)-f(i+1,j))/DX[IP1];

        if(p->wet[Ip3J]>0 && p->wet[Ip2J]>0)
        q5 = (f(i+3,j)-f(i+2,j))/DX[IP2];
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
        q1 = (f(i,j-1)-f(i,j-2))/DY[JM2];

        if(p->wet[IJ]>0 && p->wet[IJm1]>0)
        q2 = (f(i,j)-f(i,j-1))/DY[JM1];

        if(p->wet[IJp1]>0 && p->wet[IJ]>0)
        q3 = (f(i,j+1)-f(i,j))/DY[JP];

        if(p->wet[IJp2]>0 && p->wet[IJp1]>0)
        q4 = (f(i,j+2)-f(i,j+1))/DY[JP1];

        if(p->wet[IJp3]>0 && p->wet[IJp2]>0)
        q5 = (f(i,j+3)-f(i,j+2))/DY[JP2];
    }
}
