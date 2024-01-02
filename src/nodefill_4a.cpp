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

#include"nodefill.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"norm_vec.h"

void nodefill::nodefill4a(lexer *p, fdm *a, ghostcell *pgc, field &f, field &eta) 
{
    pip=5;
	TPLOOP
	eta(i,j,k)=0.0;

	double val,factor,denom;
	int q;

	ALOOP
	{
		val = f(i,j,k);
		
		eta(i,j,k) += 0.125*val;
		eta(i,j-1,k) += 0.125*val;
		eta(i,j,k-1) += 0.125*val;
		eta(i,j-1,k-1) += 0.125*val;
		eta(i-1,j,k) += 0.125*val;
		eta(i-1,j-1,k) += 0.125*val;
		eta(i-1,j,k-1) += 0.125*val;
		eta(i-1,j-1,k-1) += 0.125*val;
	}

	
	GC4ALOOP
	{
		i=p->gcb4a[n][0];
		j=p->gcb4a[n][1];
		k=p->gcb4a[n][2];
		
		if(p->gcb4a[n][3]==1)
		{
		val = f(i-1,j,k);
		
		eta(i-1,j,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j-1,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j-1,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		}
			
		if(p->gcb4a[n][3]==4)
		{
		val = f(i+1,j,k);
		
		eta(i,j,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i,j-1,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i,j,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i,j-1,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		}
		
		
		if(p->gcb4a[n][3]==3)
		{
		val = f(i,j-1,k);
		
		eta(i,j-1,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j-1,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i,j-1,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j-1,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		}
			
		if(p->gcb4a[n][3]==2)
		{
		val = f(i,j+1,k);
		
		eta(i,j,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i,j,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		}
		
		
		if(p->gcb4a[n][3]==5)
		{
		val = f(i,j,k-1);
		
		eta(i,j,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i,j-1,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j-1,k-1) += 0.125*val * (1.0/double(p->gcside4[n]));
		}
			
		if(p->gcb4a[n][3]==6)
		{
		val = f(i,j,k+1);
		
		eta(i,j,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i,j-1,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		eta(i-1,j-1,k) += 0.125*val * (1.0/double(p->gcside4[n]));
		}
	}
	
	// Para
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
	
		val = f(i-1,j,k);
		
		eta(i-1,j,k) += 0.125*val;
		eta(i-1,j-1,k) += 0.125*val;
		eta(i-1,j,k-1) += 0.125*val;
		eta(i-1,j-1,k-1) += 0.125*val;
    }
	
	for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];
	
		val = f(i+1,j,k);
		
		eta(i,j,k) += 0.125*val;
		eta(i,j-1,k) += 0.125*val;
		eta(i,j,k-1) += 0.125*val;
		eta(i,j-1,k-1) += 0.125*val;
    }
	
	for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];
	
		val = f(i,j-1,k);
		
		eta(i,j-1,k) += 0.125*val;
		eta(i-1,j-1,k) += 0.125*val;
		eta(i,j-1,k-1) += 0.125*val;
		eta(i-1,j-1,k-1) += 0.125*val;
    }
	
	for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];
	
		val = f(i,j+1,k);
		
		eta(i,j,k) += 0.125*val;
		eta(i-1,j,k) += 0.125*val;
		eta(i,j,k-1) += 0.125*val;
		eta(i-1,j,k-1) += 0.125*val;
    }
	
	for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
	
		val = f(i,j,k-1);
		
		eta(i,j,k-1) += 0.125*val;
		eta(i-1,j,k-1) += 0.125*val;
		eta(i,j-1,k-1) += 0.125*val;
		eta(i-1,j-1,k-1) += 0.125*val;
    }
	
	for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];
	
		val = f(i,j,k+1);
		
		eta(i,j,k) += 0.125*val;
		eta(i-1,j,k) += 0.125*val;
		eta(i,j-1,k) += 0.125*val;
		eta(i-1,j-1,k) += 0.125*val;
    }
	
	
// Paraco
	int sidesum=0;
	int aa,bb,cc;
	
	// 1
	for(n=0;n<p->gcparaco1_count;++n)
    {
    i=p->gcparaco1[n][0];
    j=p->gcparaco1[n][1];
    k=p->gcparaco1[n][2];
	
	aa=bb=cc=0;
	
	if(p->gcparaco1[n][3]==3 || p->gcparaco1[n][4]==3)
	bb=-1;
	
	if(p->gcparaco1[n][3]==2 || p->gcparaco1[n][4]==2)
	bb=1;
	
	if(p->gcparaco1[n][3]==5 || p->gcparaco1[n][4]==5)
	cc=-1;
	
	if(p->gcparaco1[n][3]==6 || p->gcparaco1[n][4]==6)
	cc=1;
	
	sidesum= fabs(aa) + fabs(bb) + fabs(cc);
		
		val = f(i-1,j,k);
		factor = 1.0/double(p->gcparaco1[n][5]);
        
        denom = 1.0e20;
		
		if(double(p->gcparaco1[n][5])==1)
		denom = 1.0;
		
		if(double(p->gcparaco1[n][5])==2)
		denom = 2.0;
		
		if(double(p->gcparaco1[n][5])==3)
		denom = 3.0;
		
		if(sidesum==1)
		{
			if(bb==-1)
			{
			eta(i-1,j,k-1) += 0.125*factor*val;
			eta(i-1,j,k)   += 0.125*factor*val;
			}
			
			if(bb==1)
			{
			eta(i-1,j-1,k-1) += 0.125*factor*val;
			eta(i-1,j-1,k)   += 0.125*factor*val;
			}
			
			if(cc==-1)
			{
			eta(i-1,j,k) 	+= 0.125*factor*val;
			eta(i-1,j-1,k)   += 0.125*factor*val;
			}
			
			if(cc==1)
			{
			eta(i-1,j,k-1) += 0.125*factor*val;
			eta(i-1,j-1,k-1) += 0.125*factor*val;
			}
		}
		
		if(sidesum==2)
		{
			if(bb==-1 && cc==-1)
			eta(i-1,j,k) += 0.125*val/denom;
			
			if(bb==-1 && cc==1)
			eta(i-1,j,k-1) += 0.125*val/denom;
			
			if(bb==1 && cc==-1)
			eta(i-1,j-1,k) += 0.125*val/denom;
			
			if(bb==1 && cc==1)
			eta(i-1,j-1,k-1) += 0.125*val/denom;
		}
    }
	
	// 4
	for(n=0;n<p->gcparaco4_count;++n)
    {
    i=p->gcparaco4[n][0];
    j=p->gcparaco4[n][1];
    k=p->gcparaco4[n][2];
	
	aa=bb=cc=0;
	
	if(p->gcparaco4[n][3]==3 || p->gcparaco4[n][4]==3)
	bb=-1;
	
	if(p->gcparaco4[n][3]==2 || p->gcparaco4[n][4]==2)
	bb=1;
	
	if(p->gcparaco4[n][3]==5 || p->gcparaco4[n][4]==5)
	cc=-1;
	
	if(p->gcparaco4[n][3]==6 || p->gcparaco4[n][4]==6)
	cc=1;
	
	sidesum= fabs(aa) + fabs(bb) + fabs(cc);
	
		val = f(i+1,j,k);
		factor = 1.0/double(p->gcparaco4[n][5]);
		
        denom = 1.0e20;
        
		if(double(p->gcparaco4[n][5])==1)
		denom = 1.0;
		
		if(double(p->gcparaco4[n][5])==2)
		denom = 2.0;
		
		if(double(p->gcparaco4[n][5])==3)
		denom = 3.0;
		
		
		if(sidesum==1)
		{
			if(bb==-1)
			{
			eta(i,j,k-1) += 0.125*factor*val;
			eta(i,j,k)   += 0.125*factor*val;
			}
			
			if(bb==1)
			{
			eta(i,j-1,k-1) += 0.125*factor*val;
			eta(i,j-1,k)   += 0.125*factor*val;
			}
			
			if(cc==-1)
			{
			eta(i,j,k) 		+= 0.125*factor*val;
			eta(i,j-1,k)    += 0.125*factor*val;
			}
			
			if(cc==1)
			{
			eta(i,j,k-1) 	+= 0.125*factor*val;
			eta(i,j-1,k-1)  += 0.125*factor*val;
			}
		}
		
		if(sidesum==2)
		{
			if(bb==-1 && cc==-1)
			eta(i,j,k) += 0.125*val/denom;
			
			if(bb==-1 && cc==1)
			eta(i,j,k-1) += 0.125*val/denom;
			
			if(bb==1 && cc==-1)
			eta(i,j-1,k) += 0.125*val/denom;
			
			if(bb==1 && cc==1)
			eta(i,j-1,k-1) += 0.125*val/denom;
		}
    }
	
	
	// 3
	for(n=0;n<p->gcparaco3_count;++n)
    {
    i=p->gcparaco3[n][0];
    j=p->gcparaco3[n][1];
    k=p->gcparaco3[n][2];
	
	aa=bb=cc=0;
	
	if(p->gcparaco3[n][3]==1 || p->gcparaco3[n][4]==1)
	aa=-1;
	
	if(p->gcparaco3[n][3]==4 || p->gcparaco3[n][4]==4)
	aa=1;
	
	if(p->gcparaco3[n][3]==5 || p->gcparaco3[n][4]==5)
	cc=-1;
	
	if(p->gcparaco3[n][3]==6 || p->gcparaco3[n][4]==6)
	cc=1;
	
	sidesum= fabs(aa) + fabs(bb) + fabs(cc);
	
		val = f(i,j-1,k);
		factor = 1.0/double(p->gcparaco3[n][5]);
		
        denom = 1.0e20;
        
		if(double(p->gcparaco3[n][5])==1)
		denom = 1.0;
		
		if(double(p->gcparaco3[n][5])==2)
		denom = 2.0;
		
		if(double(p->gcparaco3[n][5])==3)
		denom = 3.0;
		
		if(sidesum==1)
		{
			if(aa==-1)
			{
			eta(i,j-1,k-1) += 0.125*factor*val;
			eta(i,j-1,k)   += 0.125*factor*val;
			}
			
			if(aa==1)
			{
			eta(i-1,j-1,k-1) += 0.125*factor*val;
			eta(i-1,j-1,k)   += 0.125*factor*val;
			}
			
			if(cc==-1)
			{
			eta(i,j-1,k) 	+= 0.125*factor*val;
			eta(i-1,j-1,k)   += 0.125*factor*val;
			}
			
			if(cc==1)
			{
			eta(i,j-1,k-1) += 0.125*factor*val;
			eta(i-1,j-1,k-1) += 0.125*factor*val;
			}
		}
		
		if(sidesum==2)
		{
			if(aa==-1 && cc==-1)
			eta(i,j-1,k) += 0.125*val/denom;
			
			if(aa==-1 && cc==1)
			eta(i,j-1,k-1) += 0.125*val/denom;
			
			if(aa==1 && cc==-1)
			eta(i-1,j-1,k) += 0.125*val/denom;
			
			if(aa==1 && cc==1)
			eta(i-1,j-1,k-1) += 0.125*val/denom;
		}
    }
	
	// 2
	for(n=0;n<p->gcparaco2_count;++n)
    {
    i=p->gcparaco2[n][0];
    j=p->gcparaco2[n][1];
    k=p->gcparaco2[n][2];
	
	aa=bb=cc=0;
	
	if(p->gcparaco2[n][3]==1 || p->gcparaco2[n][4]==1)
	aa=-1;
	
	if(p->gcparaco2[n][3]==4 || p->gcparaco2[n][4]==4)
	aa=1;
	
	if(p->gcparaco2[n][3]==5 || p->gcparaco2[n][4]==5)
	cc=-1;
	
	if(p->gcparaco2[n][3]==6 || p->gcparaco2[n][4]==6)
	cc=1;
	
	sidesum= fabs(aa) + fabs(bb) + fabs(cc);
	
		val = f(i,j+1,k);
		factor = 1.0/double(p->gcparaco2[n][5]);
		
        denom = 1.0e20;
        
		if(double(p->gcparaco2[n][5])==1)
		denom = 1.0;
		
		if(double(p->gcparaco2[n][5])==2)
		denom = 2.0;
		
		if(double(p->gcparaco2[n][5])==3)
		denom = 3.0;
		
		
		if(sidesum==1)
		{
			if(aa==-1)
			{
			eta(i,j,k-1) += 0.125*factor*val;
			eta(i,j,k)   += 0.125*factor*val;
			}
			
			if(aa==1)
			{
			eta(i-1,j,k-1) += 0.125*factor*val;
			eta(i-1,j,k)   += 0.125*factor*val;
			}
			
			if(cc==-1)
			{
			eta(i,j,k) 		+= 0.125*factor*val;
			eta(i-1,j,k)    += 0.125*factor*val;
			}
			
			if(cc==1)
			{
			eta(i,j,k-1) += 0.125*factor*val;
			eta(i-1,j,k-1) += 0.125*factor*val;
			}
		}
		
		if(sidesum==2)
		{
			if(aa==-1 && cc==-1)
			eta(i,j,k) += 0.125*val/denom;
			
			if(aa==-1 && cc==1)
			eta(i,j,k-1) += 0.125*val/denom;
			
			if(aa==1 && cc==-1)
			eta(i-1,j,k) += 0.125*val/denom;
			
			if(aa==1 && cc==1)
			eta(i-1,j,k-1) += 0.125*val/denom;
		}
    }
	
	
	// 5
	for(n=0;n<p->gcparaco5_count;++n)
    {
    i=p->gcparaco5[n][0];
    j=p->gcparaco5[n][1];
    k=p->gcparaco5[n][2];
	
	aa=bb=cc=0;
	
	if(p->gcparaco5[n][3]==1 || p->gcparaco5[n][4]==1)
	aa=-1;
	
	if(p->gcparaco5[n][3]==4 || p->gcparaco5[n][4]==4)
	aa=1;
	
	if(p->gcparaco5[n][3]==3 || p->gcparaco5[n][4]==3)
	bb=-1;
	
	if(p->gcparaco5[n][3]==2 || p->gcparaco5[n][4]==2)
	bb=1;
	
	sidesum= fabs(aa) + fabs(bb) + fabs(cc);
	
		val = f(i,j,k-1);
		factor = 1.0/double(p->gcparaco5[n][5]);
		
        denom = 1.0e20;
        
		if(double(p->gcparaco5[n][5])==1)
		denom = 1.0;
		
		if(double(p->gcparaco5[n][5])==2)
		denom = 2.0;
		
		if(double(p->gcparaco5[n][5])==3)
		denom = 3.0;
		
		if(sidesum==1)
		{
			if(aa==-1)
			{
			eta(i,j-1,k-1) 	+= 0.125*factor*val;
			eta(i,j,k-1)   += 0.125*factor*val;
			}
			
			if(aa==1)
			{
			eta(i-1,j-1,k-1) += 0.125*factor*val;
			eta(i-1,j,k-1)   += 0.125*factor*val;
			}
			
			if(bb==-1)
			{
			eta(i,j,k-1) 	+= 0.125*factor*val;
			eta(i-1,j,k-1)   += 0.125*factor*val;
			}
			
			if(bb==1)
			{
			eta(i,j-1,k-1) += 0.125*factor*val;
			eta(i-1,j-1,k-1) += 0.125*factor*val;
			}
			
		}
		
		if(sidesum==2)
		{
			if(aa==-1 && bb==-1)
			eta(i,j,k-1) += 0.125*val/denom;
			
			if(aa==-1 && bb==1)
			eta(i,j-1,k-1) += 0.125*val/denom;
			
			if(aa==1 && bb==-1)
			eta(i-1,j,k-1) += 0.125*val/denom;
			
			if(aa==1 && bb==1)
			eta(i-1,j-1,k-1) += 0.125*val/denom;
		}
    }
	
	// 6
	for(n=0;n<p->gcparaco6_count;++n)
    {
    i=p->gcparaco6[n][0];
    j=p->gcparaco6[n][1];
    k=p->gcparaco6[n][2];
	
	aa=bb=cc=0;
	
	if(p->gcparaco6[n][3]==1 || p->gcparaco6[n][4]==1)
	aa=-1;
	
	if(p->gcparaco6[n][3]==4 || p->gcparaco6[n][4]==4)
	aa=1;
	
	if(p->gcparaco6[n][3]==3 || p->gcparaco6[n][4]==3)
	bb=-1;
	
	if(p->gcparaco6[n][3]==2 || p->gcparaco6[n][4]==2)
	bb=1;
	
	sidesum= fabs(aa) + fabs(bb) + fabs(cc);
		
		val = f(i,j,k+1);
		
		
		factor = 1.0/double(p->gcparaco6[n][5]);
		
        denom = 1.0e20;
        
		if(double(p->gcparaco6[n][5])==1)
		denom = 1.0;
		
		if(double(p->gcparaco6[n][5])==2)
		denom = 2.0;
		
		if(double(p->gcparaco6[n][5])==3)
		denom = 3.0;
		
		if(sidesum==1)
		{
			if(aa==-1)
			{
			eta(i,j-1,k) += 0.125*factor*val;
			eta(i,j,k)   += 0.125*factor*val;
			}
			
			if(aa==1)
			{
			eta(i-1,j-1,k) += 0.125*factor*val;
			eta(i-1,j,k)   += 0.125*factor*val;
			}
			
			if(bb==-1)
			{
			eta(i,j,k) 	   += 0.125*factor*val;
			eta(i-1,j,k)   += 0.125*factor*val;
			}
			
			if(bb==1)
			{
			eta(i,j-1,k) += 0.125*factor*val;
			eta(i-1,j-1,k) += 0.125*factor*val;
			}
		}
		
		if(sidesum==2)
		{
			if(aa==-1 && bb==-1)
			eta(i,j,k) += 0.125*val/denom;
			
			if(aa==-1 && bb==1)
			eta(i,j-1,k) += 0.125*val/denom;
			
			if(aa==1 && bb==-1)
			eta(i-1,j,k) += 0.125*val/denom;
			
			if(aa==1 && bb==1)
			eta(i-1,j-1,k) += 0.125*val/denom;
		}
    }
	
	
	//----------------------------------------------------------------
	// DGC
    int di,dj,dk;
    int dii,djj,dkk;

    for(n=0;n<p->dgc4_count;++n)
    {
        i=p->dgc4[n][0];
        j=p->dgc4[n][1];
        k=p->dgc4[n][2];
        
        
        di=p->dgc4[n][3];
        dj=p->dgc4[n][4];
        dk=p->dgc4[n][5];
        
        
        if(di==0)
        {
        djj = dj>0?0:dj;  
        dkk = dk>0?0:dk;  
        
        eta(i-1,j+djj,k+dkk) += 0.125*f(i,j+dj,k+dk);
        eta(i,j+djj,k+dkk)   += 0.125*f(i,j+dj,k+dk);
        }
        
        if(dj==0)
        {
        dii = di>0?0:di;
        dkk = dk>0?0:dk;    
        
        eta(i+dii,j-1,k+dkk) += 0.125*f(i+di,j,k+dk);  //-------
        eta(i+dii,j,k+dkk)   += 0.125*f(i+di,j,k+dk);
        }
        
        if(dk==0)
        {
        dii = di>0?0:di; 
        djj = dj>0?0:dj;   
        
        eta(i+dii,j+djj,k-1) += 0.125*f(i+di,j+dj,k);
        eta(i+dii,j+djj,k)   += 0.125*f(i+di,j+dj,k);
        }
        
        
        if(di!=0 && dj!=0 && dk!=0)
        {
        dii = di>0?0:di; 
        djj = dj>0?0:dj;  
        dkk = dk>0?0:dk; 

        eta(i+dii,j+djj,k+dkk) += 0.125*f(i+di,j+dj,k+dk);

        }        
    }
pip=0;
}
