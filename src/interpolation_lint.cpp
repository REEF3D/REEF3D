/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"interpolation.h"
#include"field.h"
#include"lexer.h"
#include"fdm.h"

double interpolation::lint1(field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

pip=4;
    if(p->flag1[IJK]>TOPO)
    v1=b(i,j,k);
    if(p->flag1[IJp1K]>TOPO)
    v2=b(i,j+1,k);
    if(p->flag1[Ip1JK]>TOPO)
    v3=b(i+1,j,k);
    if(p->flag1[Ip1Jp1K]>TOPO)
    v4=b(i+1,j+1,k);
    if(p->flag1[IJKp1]>TOPO)
    v5=b(i,j,k+1);
    if(p->flag1[IJp1Kp1]>TOPO)
    v6=b(i,j+1,k+1);
    if(p->flag1[Ip1JKp1]>TOPO)
    v7=b(i+1,j,k+1);
    if(p->flag1[Ip1Jp1Kp1]>TOPO)
    v8=b(i+1,j+1,k+1);
pip=0;

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5 + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;

    value = wc*y1 +(1.0-wc)*y2;

    return value;
}

double interpolation::lint2(field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

pip=4;
    if(p->flag2[IJK]>TOPO)
    v1=b(i,j,k);
    if(p->flag2[IJp1K]>TOPO)
    v2=b(i,j+1,k);
    if(p->flag2[Ip1JK]>TOPO)
    v3=b(i+1,j,k);
    if(p->flag2[Ip1Jp1K]>TOPO)
    v4=b(i+1,j+1,k);
    if(p->flag2[IJKp1]>TOPO)
    v5=b(i,j,k+1);
    if(p->flag2[IJp1Kp1]>TOPO)
    v6=b(i,j+1,k+1);
    if(p->flag2[Ip1JKp1]>TOPO)
    v7=b(i+1,j,k+1);
    if(p->flag2[Ip1Jp1Kp1]>TOPO)
    v8=b(i+1,j+1,k+1);
pip=0;


    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5 + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;

    value = wc*y1 +(1.0-wc)*y2;

 return value;
}

double interpolation::lint3(field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

pip=4;
    if(p->flag3[IJK]>TOPO)
    v1=b(i,j,k);
    if(p->flag3[IJp1K]>TOPO)
    v2=b(i,j+1,k);
    if(p->flag3[Ip1JK]>TOPO)
    v3=b(i+1,j,k);
    if(p->flag3[Ip1Jp1K]>TOPO)
    v4=b(i+1,j+1,k);
    if(p->flag3[IJKp1]>TOPO)
    v5=b(i,j,k+1);
    if(p->flag3[IJp1Kp1]>TOPO)
    v6=b(i,j+1,k+1);
    if(p->flag3[Ip1JKp1]>TOPO)
    v7=b(i+1,j,k+1);
    if(p->flag3[Ip1Jp1Kp1]>TOPO)
    v8=b(i+1,j+1,k+1);
pip=0;

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5 + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;

    value = wc*y1 +(1.0-wc)*y2;

 return value;
}

double interpolation::lint4(field& f, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(p->flag4[IJK]>TOPO)
    v1=f(i,j,k);
    if(p->flag4[IJp1K]>TOPO)
    v2=f(i,j+1,k);
    if(p->flag4[Ip1JK]>TOPO)
    v3=f(i+1,j,k);
    if(p->flag4[Ip1Jp1K]>TOPO)
    v4=f(i+1,j+1,k);
    if(p->flag4[IJKp1]>TOPO)
    v5=f(i,j,k+1);
    if(p->flag4[IJp1Kp1]>TOPO)
    v6=f(i,j+1,k+1);
    if(p->flag4[Ip1JKp1]>TOPO)
    v7=f(i+1,j,k+1);
    if(p->flag4[Ip1Jp1Kp1]>TOPO)
    v8=f(i+1,j+1,k+1);
    pip=0;


    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5 + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;


    value = wc*y1 +(1.0-wc)*y2;

    pip=0;
 return value;

}

double interpolation::lint4phi(fdm *a, field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{	
    double epphi=1.6*p->DXM;
	double epphi2=0.6*p->DXM;
    v1=v2=v3=v4=v5=v6=v7=v8= p->phimean-p->pos_z();

	pip=4;
    if(a->topo(i,j,k)>-epphi && a->fb(i,j,k)>-epphi2)
    v1=b(i,j,k);
    if(a->topo(i,j+1,k)>-epphi && a->fb(i,j+1,k)>-epphi2)
    v2=b(i,j+1,k);
    if(a->topo(i+1,j,k)>-epphi && a->fb(i+1,j,k)>-epphi2)
    v3=b(i+1,j,k);
    if(a->topo(i+1,j+1,k)>-epphi && a->fb(i+1,j+1,k)>-epphi2)
    v4=b(i+1,j+1,k);
    if(a->topo(i,j,k+1)>-epphi && a->fb(i,j,k+1)>-epphi2)
    v5=b(i,j,k+1);
    if(a->topo(i,j+1,k+1)>-epphi && a->fb(i,j+1,k+1)>-epphi2)
    v6=b(i,j+1,k+1);
    if(a->topo(i+1,j,k+1)>-epphi && a->fb(i+1,j,k+1)>-epphi2)
    v7=b(i+1,j,k+1);
    if(a->topo(i+1,j+1,k+1)>-epphi && a->fb(i+1,j+1,k+1)>-epphi2)
    v8=b(i+1,j+1,k+1);
    pip=0;

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5 + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;

    value = wc*y1 +(1.0-wc)*y2;


 return value;

}

double interpolation::lint_a(field& f, int& i,int& j, int& k, double wa, double wb, double wc)
{
pip=4;

    x1 = wa*f(i,j,k)   + (1.0-wa)*f(i+1,j,k);
    x2 = wa*f(i,j+1,k) + (1.0-wa)*f(i+1,j+1,k);

    x3 = wa*f(i,j,k+1)   + (1.0-wa)*f(i+1,j,k+1);
    x4 = wa*f(i,j+1,k+1) + (1.0-wa)*f(i+1,j+1,k+1);

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;
pip=0;

    value = wc*y1 +(1.0-wc)*y2;

 return value;

}

double interpolation::lint4b(field& f, int& i,int& j, int& k, double wa, double wb, double wc)
{
     v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(p->flag4[IJK]>TOPO)
    v1=f(i,j,k);
    if(p->flag4[IJp1K]>TOPO)
    v2=f(i,j+1,k);
    if(p->flag4[Ip1JK]>TOPO)
    v3=f(i+1,j,k);
    if(p->flag4[Ip1Jp1K]>TOPO)
    v4=f(i+1,j+1,k);
    
    if(p->flag4[IJK]<=TOPO)
    v1=f(i,j,k+1);
    if(p->flag4[IJp1K]<=TOPO)
    v2=f(i,j+1,k+1);
    if(p->flag4[Ip1JK]<=TOPO)
    v3=f(i+1,j,k+1);
    if(p->flag4[Ip1Jp1K]<=TOPO)
    v4=f(i+1,j+1,k+1);
    
    v5=f(i,j,k+1);
    v6=f(i,j+1,k+1);
    v7=f(i+1,j,k+1);
    v8=f(i+1,j+1,k+1);
    
    
    
    pip=0;


    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5 + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;


    value = wc*y1 +(1.0-wc)*y2;

    pip=0;
 return value;

}


double interpolation::lint4kin(field& f, int& i,int& j, int& k, double wa, double wb, double wc)
{
    pip=4;
    v1=f(i,j,k);
    v2=f(i,j+1,k);
    v3=f(i+1,j,k);
    v4=f(i+1,j+1,k);
    v5=f(i,j,k+1);
    v6=f(i,j+1,k+1);
    v7=f(i+1,j,k+1);
    v8=f(i+1,j+1,k+1);

    
    if(p->flagsf4[IJK]<0)
    v1=f(i,j,k+1);
    if(p->flagsf4[IJp1K]<0)
    v2=f(i,j+1,k+1);
    if(p->flagsf4[Ip1JK]<0)
    v3=f(i+1,j,k+1);
    if(p->flagsf4[IJp1K]<0)
    v4=f(i+1,j+1,k+1);
    if(p->flagsf4[IJKp1]<0)
    v5=f(i,j,k+2);
    if(p->flagsf4[IJp1Kp1]<0)
    v6=f(i,j+1,k+2);
    if(p->flagsf4[Ip1JKp1]<0)
    v7=f(i+1,j,k+2);
    if(p->flagsf4[Ip1Jp1Kp1]<0)
    v8=f(i+1,j+1,k+2);
    pip=0;


    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5 + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;


    value = wc*y1 +(1.0-wc)*y2;

    pip=0;
 return value;

}
