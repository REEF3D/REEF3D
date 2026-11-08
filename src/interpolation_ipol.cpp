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

#include"interpolation.h"
#include"fdm.h"
#include"field.h"
#include"lexer.h"
#include"heaviside.h"


double interpolation::ipol1(field& b)
{
    
    if(p->j_dir==0)
    {
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    if(j==-1)
    {
    if(p->flag1[IJp1K]>0)
    v1=b(i,j+1,k);
    
    if(p->flag1[IJp1Kp1]>0)
    v2=b(i,j+1,k+1);
    }
    
    if(j==0)
    {
    if(p->flag1[IJK]>0)
    v1=b(i,j,k);
    
    if(p->flag1[IJKp1]>0)
    v2=b(i,j,k+1);
    }
    
    value= 0.5*(v1+v2);
    
    }
    
    if(p->j_dir==1)
    {
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    if(p->flag1[IJK]>0)
    v1=b(i,j,k);
    if(p->flag1[IJp1K]>0)
    v2=b(i,j+1,k);
    if(p->flag1[IJKp1]>0)
    v3=b(i,j,k+1);
    if(p->flag1[IJp1Kp1]>0)
    v4=b(i,j+1,k+1);
    

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[Ip1JK]==-4)
    {
    if(p->flag1[Ip1JK]>0)
    v5=b(i+1,j,k);
    if(p->flag1[Ip1Jp1K]>0)
    v6=b(i+1,j+1,k);
    if(p->flag1[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag1[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    if(p->flag5[IJK]==-1)
    {
    if(p->flag1[Ip1JK]>0)
    v5=b(i+1,j,k);
    if(p->flag1[Ip1Jp1K]>0)
    v6=b(i+1,j+1,k);
    if(p->flag1[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag1[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);


    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }
    
    }

    return value;
}

double interpolation::ipol2( field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    if(p->flag2[IJK]>0)
    v1=b(i,j,k);
    if(p->flag2[Ip1JK]>0)
    v2=b(i+1,j,k);
    if(p->flag2[IJKp1]>0)
    v3=b(i,j,k+1);
    if(p->flag2[Ip1JKp1]>0)
    v4=b(i+1,j,k+1);

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[IJp1K]==-2)
    {
    if(p->flag2[IJp1K]>0)
    v5=b(i,j+1,k);
    if(p->flag2[Ip1Jp1K]>0)
    v6=b(i+1,j+1,k);
    if(p->flag2[IJp1Kp1]>0)
    v7=b(i,j+1,k+1);
    if(p->flag2[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    if( p->flag5[IJK]==-3)
    {
    if(p->flag2[IJp1K]>0)
    v5=b(i,j+1,k);
    if(p->flag2[Ip1Jp1K]>0)
    v6=b(i+1,j+1,k);
    if(p->flag2[IJp1Kp1]>0)
    v7=b(i,j+1,k+1);
    if(p->flag2[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    return value;
}

double interpolation::ipol3( field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;
    
    if(p->j_dir==0)
    { 
    if(j==-1)
    {
    if(p->flag3[IJp1K]>0)
    v1=b(i,j+1,k);
    if(p->flag3[Ip1Jp1K]>0)
    v2=b(i+1,j+1,k);
    }
    
    if(j==0)
    {
    if(p->flag3[IJK]>0)
    v1=b(i,j,k);
    if(p->flag3[Ip1JK]>0)
    v2=b(i+1,j,k);
    }

    value= 0.5*(v1+v2);
    }

    if(p->j_dir==1)
    { 
    if(p->flag3[IJK]>0)
    v1=b(i,j,k);
    if(p->flag3[IJp1K]>0)
    v2=b(i,j+1,k);
    if(p->flag3[Ip1JK]>0)
    v3=b(i+1,j,k);
    if(p->flag3[Ip1Jp1K]>0)
    v4=b(i+1,j+1,k);

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[IJKp1]==-6)
    {
    if(p->flag3[IJKp1]>0)
    v5=b(i,j,k+1);
    if(p->flag3[IJp1Kp1]>0)
    v6=b(i,j+1,k+1);
    if(p->flag3[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag3[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

     if(p->flag5[IJK]==-5)
    {
    if(p->flag3[IJKp1]>0)
    v5=b(i,j,k+1);
    if(p->flag3[IJp1Kp1]>0)
    v6=b(i,j+1,k+1);
    if(p->flag3[Ip1JKp1]>0)
    v7=b(i+1,j,k+1);
    if(p->flag3[Ip1Jp1Kp1]>0)
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }
    
    if(p->flag5[IJKp1]==3)
    {
    v5=b(i,j,k+1);
    v6=b(i,j+1,k+1);
    v7=b(i+1,j,k+1);
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }
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
    if(p->flag4[IJK]>0)
    v1=b(i,j,k);
    if(p->flag4[Ip1JK]>0)
    v2=b(i+1,j,k);
    if(p->flag4[IJKp1]>0)
    v3=b(i,j,k+1);
    if(p->flag4[Ip1JKp1]>0)
    v4=b(i+1,j,k+1);
    j=jj;
    
    value=0.25*(v1+v2+v3+v4);
    }
    
    if(p->j_dir==1)
    {
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

    value=0.125*(v1+v2+v3+v4+v5+v6+v7+v8);
    }

    return value;
}

double interpolation::ipol4press( field& b)
{
    double denom = 0.0;
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;


    if(p->flag4[IJK]>0)
    {
    v1=b(i,j,k);
    denom+=1.0;
    }
    
    if(p->flag4[Ip1JK]>0)
    {
    v2=b(i+1,j,k);
    denom+=1.0;
    }
    
    if(p->flag4[Ip1Jp1K]>0)
    {
    v3=b(i+1,j+1,k);
    denom+=1.0;
    }
    
    if(p->flag4[IJp1K]>0)
    {
    v4=b(i,j+1,k);
    denom+=1.0;
    }
    
    if(p->flag4[IJKp1]>0)
    {
    v5=b(i,j,k+1);
    denom+=1.0;
    }
    
    if(p->flag4[Ip1JKp1]>0)
    {
    v6=b(i+1,j,k+1);
    denom+=1.0;
    }
    
    if(p->flag4[Ip1Jp1Kp1]>0)
    {
    v7=b(i+1,j+1,k+1);
    denom+=1.0;
    }
    
    if(p->flag4[IJp1Kp1]>0)
    {
    v8=b(i,j+1,k+1);
    denom+=1.0;
    }
    
    denom = denom>0.0?denom:1.0e20;
    
	 value=(1.0/denom)*(v1+v2+v3+v4+v5+v6+v7+v8);
	 
	 return value;
}


double interpolation::ipol4phi(fdm *a, field& b)
{
    double denom = 0.0;
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    if(p->flag4[IJK]>0)
    {
    v1=b(i,j,k);
    denom+=1.0;
    }
    
    if(p->flag4[Ip1JK]>0)
    {
    v2=b(i+1,j,k);
    denom+=1.0;
    }
    
    if(p->flag4[Ip1Jp1K]>0)
    {
    v3=b(i+1,j+1,k);
    denom+=1.0;
    }
    
    if(p->flag4[IJp1K]>0)
    {
    v4=b(i,j+1,k);
    denom+=1.0;
    }
    
    if(p->flag4[IJKp1]>0)
    {
    v5=b(i,j,k+1);
    denom+=1.0;
    }
    
    if(p->flag4[Ip1JKp1]>0)
    {
    v6=b(i+1,j,k+1);
    denom+=1.0;
    }
    
    if(p->flag4[Ip1Jp1Kp1]>0)
    {
    v7=b(i+1,j+1,k+1);
    denom+=1.0;
    }
    
    if(p->flag4[IJp1Kp1]>0)
    {
    v8=b(i,j+1,k+1);
    denom+=1.0;
    }

    
	 value=(1.0/denom)*(v1+v2+v3+v4+v5+v6+v7+v8);
     
     value=0.125*(b(i,j,k)+b(i,j+1,k)+b(i+1,j,k)+b(i+1,j+1,k) +
                 b(i,j,k+1)+b(i,j+1,k+1)+b(i+1,j,k+1)+b(i+1,j+1,k+1));
	 
	 return value;
}

double interpolation::ipol4_a( field& b)
{
    if(p->j_dir==0)
    {
    jj=j;
    j=0;
    value=0.25*(b(i,j,k)+b(i+1,j,k)+b(i,j,k+1)+b(i+1,j,k+1));
    j=jj;
    }
                 
                 
    if(p->j_dir==1)
    value=0.125*(b(i,j,k)+b(i,j+1,k)+b(i+1,j,k)+b(i+1,j+1,k) +
                 b(i,j,k+1)+b(i,j+1,k+1)+b(i+1,j,k+1)+b(i+1,j+1,k+1));
                 
    return value;
}

double interpolation::ipol4_a_slice(slice& b)
{
    if(p->j_dir==0)
    {
    jj=j;
    j=0;
    value=0.5*(b(i,j)+b(i+1,j));
    j=jj;
    }
                 
                 
    if(p->j_dir==1)
    value=0.25*(b(i,j)+b(i,j+1)+b(i+1,j)+b(i+1,j+1));
                 
    return value;
}

double interpolation::ipol4_a_sliceint(sliceint& b)
{
    if(p->j_dir==0)
    {
    jj=j;
    j=0;
    value=0.5*double(b(i,j)+b(i+1,j));
    j=jj;
    }
                 
                 
    if(p->j_dir==1)
    value=0.25*double(b(i,j)+b(i,j+1)+b(i+1,j)+b(i+1,j+1));
                 
    return value;
}

// as ipol1/2/3, for the printers (VTU, LAGOON store): a face on an open boundary (inflow,
// outflow, patch: the ghost cell next to it has IO 1-3) is not a fluid face (flag < 0) but holds
// the boundary velocity; ipol1/2/3 leave it out, so the velocity on the in- and outflow
// boundary showed 0. The flow solvers and probes keep ipol1/2/3.
namespace
{
inline bool ipol_open(lexer *p, int a, int b)
{
    return (p->flag4[a]>0 && p->flag4[b]<0 && p->IO[b]>=1 && p->IO[b]<=3)
        || (p->flag4[b]>0 && p->flag4[a]<0 && p->IO[a]>=1 && p->IO[a]<=3);
}
}

#define F1OPEN(n) (p->flag1[n]>0 || ipol_open(p,(n),(n)+p->jmax*p->kmax))
#define F2OPEN(n) (p->flag2[n]>0 || ipol_open(p,(n),(n)+p->kmax))
#define F3OPEN(n) (p->flag3[n]>0 || ipol_open(p,(n),(n)+1))

double interpolation::ipol1_print(field& b)
{
    
    if(p->j_dir==0)
    {
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    if(j==-1)
    {
    if(F1OPEN(IJp1K))
    v1=b(i,j+1,k);
    
    if(F1OPEN(IJp1Kp1))
    v2=b(i,j+1,k+1);
    }
    
    if(j==0)
    {
    if(F1OPEN(IJK))
    v1=b(i,j,k);
    
    if(F1OPEN(IJKp1))
    v2=b(i,j,k+1);
    }
    
    value= 0.5*(v1+v2);
    
    }
    
    if(p->j_dir==1)
    {
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    if(F1OPEN(IJK))
    v1=b(i,j,k);
    if(F1OPEN(IJp1K))
    v2=b(i,j+1,k);
    if(F1OPEN(IJKp1))
    v3=b(i,j,k+1);
    if(F1OPEN(IJp1Kp1))
    v4=b(i,j+1,k+1);
    

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[Ip1JK]==-4 && !F1OPEN(IJK))
    {
    if(F1OPEN(Ip1JK))
    v5=b(i+1,j,k);
    if(F1OPEN(Ip1Jp1K))
    v6=b(i+1,j+1,k);
    if(F1OPEN(Ip1JKp1))
    v7=b(i+1,j,k+1);
    if(F1OPEN(Ip1Jp1Kp1))
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    if(p->flag5[IJK]==-1 && !F1OPEN(IJK))
    {
    if(F1OPEN(Ip1JK))
    v5=b(i+1,j,k);
    if(F1OPEN(Ip1Jp1K))
    v6=b(i+1,j+1,k);
    if(F1OPEN(Ip1JKp1))
    v7=b(i+1,j,k+1);
    if(F1OPEN(Ip1Jp1Kp1))
    v8=b(i+1,j+1,k+1);


    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }
    
    }

    return value;
}

double interpolation::ipol2_print(field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    if(F2OPEN(IJK))
    v1=b(i,j,k);
    if(F2OPEN(Ip1JK))
    v2=b(i+1,j,k);
    if(F2OPEN(IJKp1))
    v3=b(i,j,k+1);
    if(F2OPEN(Ip1JKp1))
    v4=b(i+1,j,k+1);

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[IJp1K]==-2 && !F2OPEN(IJK))
    {
    if(F2OPEN(IJp1K))
    v5=b(i,j+1,k);
    if(F2OPEN(Ip1Jp1K))
    v6=b(i+1,j+1,k);
    if(F2OPEN(IJp1Kp1))
    v7=b(i,j+1,k+1);
    if(F2OPEN(Ip1Jp1Kp1))
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    if(p->flag5[IJK]==-3 && !F2OPEN(IJK))
    {
    if(F2OPEN(IJp1K))
    v5=b(i,j+1,k);
    if(F2OPEN(Ip1Jp1K))
    v6=b(i+1,j+1,k);
    if(F2OPEN(IJp1Kp1))
    v7=b(i,j+1,k+1);
    if(F2OPEN(Ip1Jp1Kp1))
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

    return value;
}

double interpolation::ipol3_print(field& b)
{
    v1=v2=v3=v4=v5=v6=v7=v8=0.0;
    
    if(p->j_dir==0)
    { 
    if(j==-1)
    {
    if(F3OPEN(IJp1K))
    v1=b(i,j+1,k);
    if(F3OPEN(Ip1Jp1K))
    v2=b(i+1,j+1,k);
    }
    
    if(j==0)
    {
    if(F3OPEN(IJK))
    v1=b(i,j,k);
    if(F3OPEN(Ip1JK))
    v2=b(i+1,j,k);
    }

    value= 0.5*(v1+v2);
    }

    if(p->j_dir==1)
    { 
    if(F3OPEN(IJK))
    v1=b(i,j,k);
    if(F3OPEN(IJp1K))
    v2=b(i,j+1,k);
    if(F3OPEN(Ip1JK))
    v3=b(i+1,j,k);
    if(F3OPEN(Ip1Jp1K))
    v4=b(i+1,j+1,k);

    value= 0.25*(v1+v2+v3+v4);

    if(p->flag5[IJKp1]==-6 && !F3OPEN(IJK))
    {
    if(F3OPEN(IJKp1))
    v5=b(i,j,k+1);
    if(F3OPEN(IJp1Kp1))
    v6=b(i,j+1,k+1);
    if(F3OPEN(Ip1JKp1))
    v7=b(i+1,j,k+1);
    if(F3OPEN(Ip1Jp1Kp1))
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }

     if(p->flag5[IJK]==-5 && !F3OPEN(IJK))
    {
    if(F3OPEN(IJKp1))
    v5=b(i,j,k+1);
    if(F3OPEN(IJp1Kp1))
    v6=b(i,j+1,k+1);
    if(F3OPEN(Ip1JKp1))
    v7=b(i+1,j,k+1);
    if(F3OPEN(Ip1Jp1Kp1))
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }
    
    if(p->flag5[IJKp1]==3 && !F3OPEN(IJK))
    {
    v5=b(i,j,k+1);
    v6=b(i,j+1,k+1);
    v7=b(i+1,j,k+1);
    v8=b(i+1,j+1,k+1);

    value= 0.5*(value + 0.25*(v5+v6+v7+v8));
    }
    }

    return value;
}

#undef F1OPEN
#undef F2OPEN
#undef F3OPEN
