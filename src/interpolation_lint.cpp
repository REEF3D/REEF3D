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

#include "fdm.h"
#include "field.h"
#include "interpolation.h"
#include "lexer.h"

double interpolation::lint1(field &b, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag1[IJK] > TOPO_FLAG)
        v1 = b(i, j, k);
    if(p->flag1[IJp1K] > TOPO_FLAG)
        v2 = b(i, j + 1, k);
    if(p->flag1[Ip1JK] > TOPO_FLAG)
        v3 = b(i + 1, j, k);
    if(p->flag1[Ip1Jp1K] > TOPO_FLAG)
        v4 = b(i + 1, j + 1, k);
    if(p->flag1[IJKp1] > TOPO_FLAG)
        v5 = b(i, j, k + 1);
    if(p->flag1[IJp1Kp1] > TOPO_FLAG)
        v6 = b(i, j + 1, k + 1);
    if(p->flag1[Ip1JKp1] > TOPO_FLAG)
        v7 = b(i + 1, j, k + 1);
    if(p->flag1[Ip1Jp1Kp1] > TOPO_FLAG)
        v8 = b(i + 1, j + 1, k + 1);

    // cout<<" v1: "<<v1<<" v3: "<<v3<<" ";

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

double interpolation::lint1c(field &b, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag1[IJK] > OBJ_FLAG)
        v1 = b(i, j, k);
    if(p->flag1[IJp1K] > OBJ_FLAG)
        v2 = b(i, j + 1, k);
    if(p->flag1[Ip1JK] > OBJ_FLAG)
        v3 = b(i + 1, j, k);
    if(p->flag1[Ip1Jp1K] > OBJ_FLAG)
        v4 = b(i + 1, j + 1, k);
    if(p->flag1[IJKp1] > OBJ_FLAG)
        v5 = b(i, j, k + 1);
    if(p->flag1[IJp1Kp1] > OBJ_FLAG)
        v6 = b(i, j + 1, k + 1);
    if(p->flag1[Ip1JKp1] > OBJ_FLAG)
        v7 = b(i + 1, j, k + 1);
    if(p->flag1[Ip1Jp1Kp1] > OBJ_FLAG)
        v8 = b(i + 1, j + 1, k + 1);

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

double interpolation::lint2(field &b, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag2[IJK] > TOPO_FLAG)
        v1 = b(i, j, k);
    if(p->flag2[IJp1K] > TOPO_FLAG)
        v2 = b(i, j + 1, k);
    if(p->flag2[Ip1JK] > TOPO_FLAG)
        v3 = b(i + 1, j, k);
    if(p->flag2[Ip1Jp1K] > TOPO_FLAG)
        v4 = b(i + 1, j + 1, k);
    if(p->flag2[IJKp1] > TOPO_FLAG)
        v5 = b(i, j, k + 1);
    if(p->flag2[IJp1Kp1] > TOPO_FLAG)
        v6 = b(i, j + 1, k + 1);
    if(p->flag2[Ip1JKp1] > TOPO_FLAG)
        v7 = b(i + 1, j, k + 1);
    if(p->flag2[Ip1Jp1Kp1] > TOPO_FLAG)
        v8 = b(i + 1, j + 1, k + 1);

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

double interpolation::lint2c(field &b, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag2[IJK] > OBJ_FLAG)
        v1 = b(i, j, k);
    if(p->flag2[IJp1K] > OBJ_FLAG)
        v2 = b(i, j + 1, k);
    if(p->flag2[Ip1JK] > OBJ_FLAG)
        v3 = b(i + 1, j, k);
    if(p->flag2[Ip1Jp1K] > OBJ_FLAG)
        v4 = b(i + 1, j + 1, k);
    if(p->flag2[IJKp1] > OBJ_FLAG)
        v5 = b(i, j, k + 1);
    if(p->flag2[IJp1Kp1] > OBJ_FLAG)
        v6 = b(i, j + 1, k + 1);
    if(p->flag2[Ip1JKp1] > OBJ_FLAG)
        v7 = b(i + 1, j, k + 1);
    if(p->flag2[Ip1Jp1Kp1] > OBJ_FLAG)
        v8 = b(i + 1, j + 1, k + 1);

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

double interpolation::lint3(field &b, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag3[IJK] > TOPO_FLAG)
        v1 = b(i, j, k);
    if(p->flag3[IJp1K] > TOPO_FLAG)
        v2 = b(i, j + 1, k);
    if(p->flag3[Ip1JK] > TOPO_FLAG)
        v3 = b(i + 1, j, k);
    if(p->flag3[Ip1Jp1K] > TOPO_FLAG)
        v4 = b(i + 1, j + 1, k);
    if(p->flag3[IJKp1] > TOPO_FLAG)
        v5 = b(i, j, k + 1);
    if(p->flag3[IJp1Kp1] > TOPO_FLAG)
        v6 = b(i, j + 1, k + 1);
    if(p->flag3[Ip1JKp1] > TOPO_FLAG)
        v7 = b(i + 1, j, k + 1);
    if(p->flag3[Ip1Jp1Kp1] > TOPO_FLAG)
        v8 = b(i + 1, j + 1, k + 1);

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

double interpolation::lint3c(field &b, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag3[IJK] > OBJ_FLAG)
        v1 = b(i, j, k);
    if(p->flag3[IJp1K] > OBJ_FLAG)
        v2 = b(i, j + 1, k);
    if(p->flag3[Ip1JK] > OBJ_FLAG)
        v3 = b(i + 1, j, k);
    if(p->flag3[Ip1Jp1K] > OBJ_FLAG)
        v4 = b(i + 1, j + 1, k);
    if(p->flag3[IJKp1] > OBJ_FLAG)
        v5 = b(i, j, k + 1);
    if(p->flag3[IJp1Kp1] > OBJ_FLAG)
        v6 = b(i, j + 1, k + 1);
    if(p->flag3[Ip1JKp1] > OBJ_FLAG)
        v7 = b(i + 1, j, k + 1);
    if(p->flag3[Ip1Jp1Kp1] > OBJ_FLAG)
        v8 = b(i + 1, j + 1, k + 1);

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

double interpolation::lint4(field &f, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag4[IJK] > TOPO_FLAG)
        v1 = f(i, j, k);
    if(p->flag4[IJp1K] > TOPO_FLAG)
        v2 = f(i, j + 1, k);
    if(p->flag4[Ip1JK] > TOPO_FLAG)
        v3 = f(i + 1, j, k);
    if(p->flag4[Ip1Jp1K] > TOPO_FLAG)
        v4 = f(i + 1, j + 1, k);
    if(p->flag4[IJKp1] > TOPO_FLAG)
        v5 = f(i, j, k + 1);
    if(p->flag4[IJp1Kp1] > TOPO_FLAG)
        v6 = f(i, j + 1, k + 1);
    if(p->flag4[Ip1JKp1] > TOPO_FLAG)
        v7 = f(i + 1, j, k + 1);
    if(p->flag4[Ip1Jp1Kp1] > TOPO_FLAG)
        v8 = f(i + 1, j + 1, k + 1);

    // cout<<" v1: "<<v1<<" v3: "<<v3<<" ";

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

double interpolation::lint_a(field &f, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag2[IJK] > OBJ_FLAG)
        v1 = b(i, j, k);
    if(p->flag2[IJp1K] > OBJ_FLAG)
        v2 = b(i, j + 1, k);
    if(p->flag2[Ip1JK] > OBJ_FLAG)
        v3 = b(i + 1, j, k);
    if(p->flag2[Ip1Jp1K] > OBJ_FLAG)
        v4 = b(i + 1, j + 1, k);
    if(p->flag2[IJKp1] > OBJ_FLAG)
        v5 = b(i, j, k + 1);
    if(p->flag2[IJp1Kp1] > OBJ_FLAG)
        v6 = b(i, j + 1, k + 1);
    if(p->flag2[Ip1JKp1] > OBJ_FLAG)
        v7 = b(i + 1, j, k + 1);
    if(p->flag2[Ip1Jp1Kp1] > OBJ_FLAG)
        v8 = b(i + 1, j + 1, k + 1);

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

double interpolation::lint3(field &b, int &i, int &j, int &k, double wa, double wb, double wc)
{
    v1 = v2 = v3 = v4 = v5 = v6 = v7 = v8 = 0.0;

    if(p->flag3[IJK] > TOPO_FLAG)
        v1 = b(i, j, k);
    if(p->flag3[IJp1K] > TOPO_FLAG)
        v2 = b(i, j + 1, k);
    if(p->flag3[Ip1JK] > TOPO_FLAG)
        v3 = b(i + 1, j, k);
    if(p->flag3[Ip1Jp1K] > TOPO_FLAG)
        v4 = b(i + 1, j + 1, k);
    if(p->flag3[IJKp1] > TOPO_FLAG)
        v5 = b(i, j, k + 1);
    if(p->flag3[IJp1Kp1] > TOPO_FLAG)
        v6 = b(i, j + 1, k + 1);
    if(p->flag3[Ip1JKp1] > TOPO_FLAG)
        v7 = b(i + 1, j, k + 1);
    if(p->flag3[Ip1Jp1Kp1] > TOPO_FLAG)
        v8 = b(i + 1, j + 1, k + 1);

    x1 = wa * v1 + (1.0 - wa) * v3;
    x2 = wa * v2 + (1.0 - wa) * v4;

    x3 = wa * v5 + (1.0 - wa) * v7;
    x4 = wa * v6 + (1.0 - wa) * v8;

    y1 = wb * x1 + (1.0 - wb) * x2;
    y2 = wb * x3 + (1.0 - wb) * x4;

    value = wc * y1 + (1.0 - wc) * y2;

    return value;
}

