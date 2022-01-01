/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
--------------------------------------------------------------------*/

#include"ccipol.h"
#include"lexer.h"
#include"fdm.h"

ccipol::ccipol(lexer* p):dx(p->DXM)
{
}

ccipol::~ccipol()
{
}

double ccipol::ccipol1(fdm* a,field& f, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-1.0);
    j=int((yp)/dx-0.5);
    k=int((zp)/dx-0.5);

    wa=((double(i) + 2.0)-xp/dx);
    wb=((double(j) + 1.5)-yp/dx);
    wc=((double(k) + 1.5)-zp/dx);

    value = lint(a,f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double ccipol::ccipol2(fdm* a,field& f, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-0.5);
    j=int((yp)/dx-1.0);
    k=int((zp)/dx-0.5);

    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 2.0)-yp/dx);
    wc=((double(k) + 1.5)-zp/dx);

    value = lint(a,f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double ccipol::ccipol3(fdm* a,field& f, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-0.5);
    j=int((yp)/dx-0.5);
    k=int((zp)/dx-1.0);

    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 1.5)-yp/dx);
    wc=((double(k) + 2.0)-zp/dx);

    value = lint(a,f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double ccipol::ccipol4(fdm* a,field& f, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-0.5);
    j=int((yp)/dx-0.5);
    k=int((zp)/dx-0.5);

    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 1.5)-yp/dx);
    wc=((double(k) + 1.5)-zp/dx);

    value =  lint(a,f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double ccipol::ccipol4_a(fdm* a,field& f, double& xp, double& yp, double& zp)
{
    ii=i;
    jj=j;
    kk=k;

    i=int((xp)/dx-0.5);
    j=int((yp)/dx-0.5);
    k=int((zp)/dx-0.5);

    wa=((double(i) + 1.5)-xp/dx);
    wb=((double(j) + 1.5)-yp/dx);
    wc=((double(k) + 1.5)-zp/dx);

    value =  lint_a(f,i,j,k,wa,wb,wc);

    i=ii;
    j=jj;
    k=kk;

    return value;
}

double ccipol::lint(fdm *a,field& b, int& i,int& j, int& k, double wa, double wb, double wc)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(a->topo(i,j,k)>0.0)
    v1=b(i,j,k);
    if(a->topo(i,j+1,k)>0.0)
    v2=b(i,j+1,k);
    if(a->topo(i+1,j,k)>0.0)
    v3=b(i+1,j,k);
    if(a->topo(i+1,j+1,k)>0.0)
    v4=b(i+1,j+1,k);
    if(a->topo(i,j,k+1)>0.0)
    v5=b(i,j,k+1);
    if(a->topo(i,j+1,k+1)>0.0)
    v6=b(i,j+1,k+1);
    if(a->topo(i+1,j,k+1)>0.0)
    v7=b(i+1,j,k+1);
    if(a->topo(i+1,j+1,k+1)>0.0)
    v8=b(i+1,j+1,k+1);
    pip=0;


pip=4;

    x1 = wa*v1   + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5   + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;

    value = wc*y1 +(1.0-wc)*y2;

pip=0;
 return value;

}

double ccipol::lint_a(field& f, int& i,int& j, int& k, double wa, double wb, double wc)
{
pip=4;

    x1 = wa*f(i,j,k)   + (1.0-wa)*f(i+1,j,k);
    x2 = wa*f(i,j+1,k) + (1.0-wa)*f(i+1,j+1,k);

    x3 = wa*f(i,j,k+1)   + (1.0-wa)*f(i+1,j,k+1);
    x4 = wa*f(i,j+1,k+1) + (1.0-wa)*f(i+1,j+1,k+1);

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;

    value = wc*y1 +(1.0-wc)*y2;

pip=0;
 return value;

}

double ccipol::ipol1(fdm* a,lexer* p, field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(a->topo(i,j,k)>0.0)
    v1=b(i,j,k);
    if(a->topo(i,j+1,k)>0.0)
    v2=b(i,j+1,k);
    if(a->topo(i,j,k+1)>0.0)
    v3=b(i,j,k+1);
    if(a->topo(i,j+1,k+1)>0.0)
    v4=b(i,j+1,k+1);
    pip=0;

    val= 0.25*(v1+v2+v3+v4);

    if(p->flag5[Ip1JK]==-4)
    {
    pip=4;
    if(a->topo(i+1,j,k)>0.0)
    v5=b(i+1,j,k);
    if(a->topo(i+1,j+1,k)>0.0)
    v6=b(i+1,j+1,k);
    if(a->topo(i+1,j,k+1)>0.0)
    v7=b(i+1,j,k+1);
    if(a->topo(i+1,j+1,k+1)>0.0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    val= 0.5*(val + 0.25*(v5+v6+v7+v8));
    }

    if(p->flag5[IJK]==-1)
    {
    pip=4;
    if(a->topo(i+1,j,k)>0.0)
    v5=b(i+1,j,k);
    if(a->topo(i+1,j+1,k)>0.0)
    v6=b(i+1,j+1,k);
    if(a->topo(i+1,j,k+1)>0.0)
    v7=b(i+1,j,k+1);
    if(a->topo(i+1,j+1,k+1)>0.0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    val= 0.5*(val + 0.25*(v5+v6+v7+v8));
    }


    return val;
}

double ccipol::ipol2(fdm* a,lexer* p, field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(a->topo(i,j,k)>0.0)
    v1=b(i,j,k);
    if(a->topo(i+1,j,k)>0.0)
    v2=b(i+1,j,k);
    if(a->topo(i,j,k+1)>0.0)
    v3=b(i,j,k+1);
    if(a->topo(i+1,j,k+1)>0.0)
    v4=b(i+1,j,k+1);
    pip=0;

    val= 0.25*(v1+v2+v3+v4);

    if(p->flag5[IJp1K]==-2)
    {
    pip=4;
    if(a->topo(i,j+1,k)>0.0)
    v5=b(i,j+1,k);
    if(a->topo(i+1,j+1,k)>0.0)
    v6=b(i+1,j+1,k);
    if(a->topo(i,j+1,k+1)>0.0)
    v7=b(i,j+1,k+1);
    if(a->topo(i+1,j+1,k+1)>0.0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    val= 0.5*(val + 0.25*(v5+v6+v7+v8));
    }

    if( p->flag5[IJK]==-3)
    {
    pip=4;
    if(a->topo(i,j+1,k)>0.0)
    v5=b(i,j+1,k);
    if(a->topo(i+1,j+1,k)>0.0)
    v6=b(i+1,j+1,k);
    if(a->topo(i,j+1,k+1)>0.0)
    v7=b(i,j+1,k+1);
    if(a->topo(i+1,j+1,k+1)>0.0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    val= 0.5*(val + 0.25*(v5+v6+v7+v8));
    }

    return val;
}

double ccipol::ipol3(fdm* a,lexer* p, field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(a->topo(i,j,k)>0.0)
    v1=b(i,j,k);
    if(a->topo(i,j+1,k)>0.0)
    v2=b(i,j+1,k);
    if(a->topo(i+1,j,k)>0.0)
    v3=b(i+1,j,k);
    if(a->topo(i+1,j+1,k)>0.0)
    v4=b(i+1,j+1,k);
    pip=0;

    val= 0.25*(v1+v2+v3+v4);

    if(p->flag5[IJKp1]==-6)
    {
     pip=4;
    if(a->topo(i,j,k+1)>0.0)
    v5=b(i,j,k+1);
    if(a->topo(i,j+1,k+1)>0.0)
    v6=b(i,j+1,k+1);
    if(a->topo(i+1,j,k+1)>0.0)
    v7=b(i+1,j,k+1);
    if(a->topo(i+1,j+1,k+1)>0.0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    val= 0.5*(val + 0.25*(v5+v6+v7+v8));
    }

     if(p->flag5[IJK]==-5)
    {
     pip=4;
    if(a->topo(i,j,k+1)>0.0)
    v5=b(i,j,k+1);
    if(a->topo(i,j+1,k+1)>0.0)
    v6=b(i,j+1,k+1);
    if(a->topo(i+1,j,k+1)>0.0)
    v7=b(i+1,j,k+1);
    if(a->topo(i+1,j+1,k+1)>0.0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    val= 0.5*(val + 0.25*(v5+v6+v7+v8));
    }

    return val;
}


double ccipol::ipol4(fdm* a,lexer* p, field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(a->topo(i,j,k)>0.0 && a->fb(i,j,k)>0.0)
    v1=b(i,j,k);
    if(a->topo(i,j+1,k)>0.0 && a->fb(i,j+1,k)>0.0)
    v2=b(i,j+1,k);
    if(a->topo(i+1,j,k)>0.0 && a->fb(i+1,j,k)>0.0)
    v3=b(i+1,j,k);
    if(a->topo(i+1,j+1,k)>0.0 && a->fb(i+1,j+1,k)>0.0)
    v4=b(i+1,j+1,k);
    if(a->topo(i,j,k+1)>0.0 && a->fb(i,j,k+1)>0.0)
    v5=b(i,j,k+1);
    if(a->topo(i,j+1,k+1)>0.0 && a->fb(i,j+1,k+1)>0.0)
    v6=b(i,j+1,k+1);
    if(a->topo(i+1,j,k+1)>0.0 && a->fb(i+1,j,k+1)>0.0)
    v7=b(i+1,j,k+1);
    if(a->topo(i+1,j+1,k+1)>0.0 && a->fb(i+1,j+1,k+1)>0.0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    val=0.125*(v1+v2+v3+v4+v5+v6+v7+v8);
/*
    pip=4;
    val=0.125*(b(i,j,k)+b(i,j+1,k)+b(i+1,j,k)+b(i+1,j+1,k) +b(i,j,k+1)+b(i,j+1,k+1)+b(i+1,j,k+1)+b(i+1,j+1,k+1));
pip=0;
*/
    return val;
}

double ccipol::ipol4_a(fdm* a,lexer* p, field& b)
{

pip=4;
    val=0.125*(b(i,j,k)+b(i,j+1,k)+b(i+1,j,k)+b(i+1,j+1,k) +b(i,j,k+1)+b(i,j+1,k+1)+b(i+1,j,k+1)+b(i+1,j+1,k+1));
pip=0;

    return val;
}


