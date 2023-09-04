/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

double interpolation::lint1_2D(field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=0.0;
    
    jj=j;
    j=0;

pip=4;
    if(p->flag1[IJK]>TOPO)
    v1=b(i,j,k);
    if(p->flag1[IJKp1]>TOPO)
    v2=b(i,j,k+1);
    if(p->flag1[Ip1JK]>TOPO)
    v3=b(i+1,j,k);
    if(p->flag1[Ip1JKp1]>TOPO)
    v4=b(i+1,j,k+1);
pip=0;
    j=jj;

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    value = wc*x1 +(1.0-wc)*x2;

    return value;
}

double interpolation::lint2_2D(field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=0.0;
    
    jj=j;
    j=0;

pip=4;
    if(p->flag2[IJK]>TOPO)
    v1=b(i,j,k);
    if(p->flag2[IJKp1]>TOPO)
    v2=b(i,j,k+1);
    if(p->flag2[Ip1JK]>TOPO)
    v3=b(i+1,j,k);
    if(p->flag2[Ip1JKp1]>TOPO)
    v4=b(i+1,j,k+1);
pip=0;
    j=jj;

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    value = wc*x1 +(1.0-wc)*x2;

 return value;
}

double interpolation::lint3_2D(field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=0.0;
    
    jj=j;
    j=0;

pip=4;
    if(p->flag3[IJK]>TOPO)
    v1=b(i,j,k);
    if(p->flag3[IJKp1]>TOPO)
    v2=b(i,j,k+1);
    if(p->flag3[Ip1JK]>TOPO)
    v3=b(i+1,j,k);
    if(p->flag3[Ip1JKp1]>TOPO)
    v4=b(i+1,j,k+1);
pip=0;
    j=jj;

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    value = wc*x1 +(1.0-wc)*x2;

    return value;
}

double interpolation::lint4_2D(field& f, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=0.0;
    
    jj=j;
    j=0;

    pip=4;
    if(p->flag4[IJK]>TOPO)
    v1=f(i,j,k);
    if(p->flag4[IJKp1]>TOPO)
    v2=f(i,j,k+1);
    if(p->flag4[Ip1JK]>TOPO)
    v3=f(i+1,j,k);
    if(p->flag4[Ip1JKp1]>TOPO)
    v4=f(i+1,j,k+1);
    pip=0;
    j=jj;

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;
    
    value = wc*x1 +(1.0-wc)*x2;

    return value;
}

double interpolation::lint4phi_2D(fdm *a, field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{	
    double epphi=1.6*p->DXM;
	double epphi2=0.6*p->DXM;
    v1=v2=v3=v4=p->phimean-p->pos_z();
    
    jj=j;
    j=0;

	pip=4;
    if(a->topo(i,j,k)>-epphi && a->fb(i,j,k)>-epphi2)
    v1=b(i,j,k);
    if(a->topo(i,j,k+1)>-epphi && a->fb(i,j,k+1)>-epphi2)
    v2=b(i,j,k+1);
    if(a->topo(i+1,j,k)>-epphi && a->fb(i+1,j,k)>-epphi2)
    v3=b(i+1,j,k);
    if(a->topo(i+1,j,k+1)>-epphi && a->fb(i+1,j,k+1)>-epphi2)
    v4=b(i+1,j,k+1);
    pip=0;
    
    j=jj;

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;
    
    value = wc*x1 +(1.0-wc)*x2;

    return value;
}

double interpolation::lint_a_2D(field& f, int& i,int& j, int& k, double wa, double wb, double wc)
{
    jj=j;
    j=0;
        
    pip=4;

    x1 = wa*f(i,j,k)   + (1.0-wa)*f(i+1,j,k);
    x2 = wa*f(i,j,k+1) + (1.0-wa)*f(i+1,j,k+1);

    pip=0;
    j=jj;

    value = wc*x1 +(1.0-wc)*x2;

 return value;

}
