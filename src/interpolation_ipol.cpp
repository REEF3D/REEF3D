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
#include"fdm.h"
#include"field.h"
#include"lexer.h"


double interpolation::ipol1(field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(p->flag[IJK]>0)
    v1=b(i,j,k);
    if(p->flag[IJp1K]>0)
    v2=b(i,j+1,k);
    if(p->flag[IJKp1]>0)
    v3=b(i,j,k+1);
    if(p->flag[IJp1Kp1]>0)
    v4=b(i,j+1,k+1);
    pip=0;

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[Ip1JK]==-4)
    {
    pip=4;
    if(p->flag[Ip1JK]>0)
    v5=b(i+1,j,k);
    if(p->flag[Ip1Jp1K]>0)
    v6=b(i+1,j+1,k);
    if(p->flag[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    if(p->flag5[IJK]==-1)
    {
    pip=4;
    if(p->flag[Ip1JK]>0)
    v5=b(i+1,j,k);
    if(p->flag[Ip1Jp1K]>0)
    v6=b(i+1,j+1,k);
    if(p->flag[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    return value;
}

double interpolation::ipol2( field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(p->flag[IJK]>0)
    v1=b(i,j,k);
    if(p->flag[Ip1JK]>0)
    v2=b(i+1,j,k);
    if(p->flag[IJKp1]>0)
    v3=b(i,j,k+1);
    if(p->flag[Ip1JKp1]>0)
    v4=b(i+1,j,k+1);
    pip=0;

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[IJp1K]==-2)
    {
    pip=4;
    if(p->flag[IJp1K]>0)
    v5=b(i,j+1,k);
    if(p->flag[Ip1Jp1K]>0)
    v6=b(i+1,j+1,k);
    if(p->flag[IJp1Kp1]>0)
    v7=b(i,j+1,k+1);
    if(p->flag[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    if( p->flag5[IJK]==-3)
    {
    pip=4;
    if(p->flag[IJp1K]>0)
    v5=b(i,j+1,k);
    if(p->flag[Ip1Jp1K]>0)
    v6=b(i+1,j+1,k);
    if(p->flag[IJp1Kp1]>0)
    v7=b(i,j+1,k+1);
    if(p->flag[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    return value;
}

double interpolation::ipol3( field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(p->flag[IJK]>0)
    v1=b(i,j,k);
    if(p->flag[IJp1K]>0)
    v2=b(i,j+1,k);
    if(p->flag[Ip1JK]>0)
    v3=b(i+1,j,k);
    if(p->flag[Ip1Jp1K]>0)
    v4=b(i+1,j+1,k);
    pip=0;

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[IJKp1]==-6)
    {
     pip=4;
    if(p->flag[IJKp1]>0)
    v5=b(i,j,k+1);
    if(p->flag[IJp1Kp1]>0)
    v6=b(i,j+1,k+1);
    if(p->flag[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

     if(p->flag5[IJK]==-5)
    {
     pip=4;
    if(p->flag[IJKp1]>0)
    v5=b(i,j,k+1);
    if(p->flag[IJp1Kp1]>0)
    v6=b(i,j+1,k+1);
    if(p->flag[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }
    
    if(p->flag5[IJKp1]==3)
    {
     pip=4;
    v5=b(i,j,k+1);
    v6=b(i,j+1,k+1);
    v7=b(i+1,j,k+1);
    v8=b(i+1,j+1,k+1);
    pip=0;

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    return value;
}

double interpolation::ipol4( field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;
    
    
    
    if(p->j_dir==0)
    {
    jj=j;
    j=0;
    pip=4;
    if(p->flag4[IJK]>0)
    v1=b(i,j,k);
    if(p->flag4[Ip1JK]>0)
    v2=b(i+1,j,k);
    if(p->flag4[IJKp1]>0)
    v3=b(i,j,k+1);
    if(p->flag4[Ip1JKp1]>0)
    v4=b(i+1,j,k+1);
    pip=0;
    j=jj;
    
    value=0.25*(v1+v2+v3+v4);
    }
    
    if(p->j_dir==1)
    {
    pip=4;
    if(p->flag4[IJK]>0)
    v1=b(i,j,k);
    if(p->flag4[IJp1K]>0)
    v2=b(i,j+1,k);
    if(p->flag4[Ip1JK]>0)
    v3=b(i+1,j,k);
    if(p->flag4[Ip1Jp1K]>0)
    v4=b(i+1,j+1,k);
    if(p->flag4[IJKp1]>0)
    v5=b(i,j,k+1);
    if(p->flag4[IJp1Kp1]>0)
    v6=b(i,j+1,k+1);
    if(p->flag4[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag4[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);
    pip=0;

    value=0.125*(v1+v2+v3+v4+v5+v6+v7+v8);
    }

    return value;
}

double interpolation::ipol4press( field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;
    
    pip=4;
    v1=b(i,j,k);
    v2=b(i,j+1,k);
    v3=b(i+1,j,k);
    v4=b(i+1,j+1,k);
    v5=b(i,j,k+1);
    v6=b(i,j+1,k+1);
    v7=b(i+1,j,k+1);
    v8=b(i+1,j+1,k+1);
    pip=0;

    value=0.125*(v1+v2+v3+v4+v5+v6+v7+v8);

    return value;
}

double interpolation::ipol4ro(fdm *a, field& b)
{
	double phival,H,roval;
	double epsi=1.6*p->DXM;
	
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    pip=4;
    if(p->flag[IJK]>0)
    v1=a->phi(i,j,k);
    if(p->flag[IJp1K]>0)
    v2=a->phi(i,j+1,k);
    if(p->flag[Ip1JK]>0)
    v3=a->phi(i+1,j,k);
    if(p->flag[Ip1Jp1K]>0)
    v4=a->phi(i+1,j+1,k);
    if(p->flag[IJKp1]>0)
    v5=a->phi(i,j,k+1);
    if(p->flag[IJp1Kp1]>0)
    v6=a->phi(i,j+1,k+1);
    if(p->flag[Ip1JKp1]>0)
    v7=a->phi(i+1,j,k+1);
    if(p->flag[Ip1Jp1Kp1]>0)
    v8=a->phi(i+1,j+1,k+1);
    pip=0;

    phival=0.125*(v1+v2+v3+v4+v5+v6+v7+v8);
	
	if(phival>epsi)
	H=1.0;

	if(phival<-epsi)
	H=0.0;

	if(fabs(phival)<=epsi)
	H=0.5*(1.0 + phival/epsi + (1.0/PI)*sin((PI*phival)/epsi));
		
	roval = p->W1*H + p->W3*(1.0-H);
	
    return roval;
}


double interpolation::ipol4phi(fdm *a, field& b)
{
    double epphi=1.6*p->DXM;
	double epphi2=0.6*p->DXM;
    v1=v2=v3=v4=v5=v6=v7=v8 = p->phimean-p->pos_z()-0.5*p->DXM;

	pip=4;
    if(a->topo(i,j,k)>-epphi)// && a->fb(i,j,k)>-epphi2)
    v1=b(i,j,k);
    if(a->topo(i,j+1,k)>-epphi)// && a->fb(i,j+1,k)>-epphi2)
    v2=b(i,j+1,k);
    if(a->topo(i+1,j,k)>-epphi)// && a->fb(i+1,j,k)>-epphi2)
    v3=b(i+1,j,k);
    if(a->topo(i+1,j+1,k)>-epphi)// && a->fb(i+1,j+1,k)>-epphi2)
    v4=b(i+1,j+1,k);
    if(a->topo(i,j,k+1)>-epphi)// && a->fb(i,j,k+1)>-epphi2)
    v5=b(i,j,k+1);
    if(a->topo(i,j+1,k+1)>-epphi)// && a->fb(i,j+1,k+1)>-epphi2)
    v6=b(i,j+1,k+1);
    if(a->topo(i+1,j,k+1)>-epphi)// && a->fb(i+1,j,k+1)>-epphi2)
    v7=b(i+1,j,k+1);
    if(a->topo(i+1,j+1,k+1)>-epphi)// && a->fb(i+1,j+1,k+1)>-epphi2)
    v8=b(i+1,j+1,k+1);
    pip=0;
	
	 value=0.125*(v1+v2+v3+v4+v5+v6+v7+v8);
	 
	 return value;
	
}

double interpolation::ipol4_a( field& b)
{

pip=4;
    value=0.125*(b(i,j,k)+b(i,j+1,k)+b(i+1,j,k)+b(i+1,j+1,k) +
                 b(i,j,k+1)+b(i,j+1,k+1)+b(i+1,j,k+1)+b(i+1,j+1,k+1));
pip=0;

    return value;
}

double interpolation::ipol4topo(fdm *a, field& b)
{
    double epphi=2.6*p->DXM;

    v1=v2=v3=v4=v5=v6=v7=v8 = p->S57-p->pos_z()-0.5*p->DXM;

	pip=4;
    if(a->solid(i,j,k)>-epphi)
    v1=b(i,j,k);
    if(a->solid(i,j+1,k)>-epphi)
    v2=b(i,j+1,k);
    if(a->solid(i+1,j,k)>-epphi)
    v3=b(i+1,j,k);
    if(a->solid(i+1,j+1,k)>-epphi)
    v4=b(i+1,j+1,k);
    if(a->solid(i,j,k+1)>-epphi)
    v5=b(i,j,k+1);
    if(a->solid(i,j+1,k+1)>-epphi)
    v6=b(i,j+1,k+1);
    if(a->solid(i+1,j,k+1)>-epphi)
    v7=b(i+1,j,k+1);
    if(a->solid(i+1,j+1,k+1)>-epphi)
    v8=b(i+1,j+1,k+1);
    pip=0;
	
    value=0.125*(v1+v2+v3+v4+v5+v6+v7+v8);
	 
    return value;
	
}


