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

#include"mgccart.h"
#include"cart4.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fieldint5.h"

mgccart::mgccart(lexer *p)
{
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;
}

mgccart::~mgccart()
{
}

void mgccart::startmgc(lexer* p)
{
}

void mgccart::makemgc(lexer* p, int *mgc)
{

}

void mgccart::mgcsetup(lexer* p, int *mgc)
{
	for(i=0;i<imax*jmax*kmax;++i)
	mgc[i]=0;

	LOOP
	mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin]=1;
}

void mgccart::fillmgc(lexer* p, int gcb_count, int **gcb, int *mgc)
{
	int q,n;
	p->gcextra4=10;
//--------------------------
//WALL1
	QGC4LOOP
	{
	    i=p->gcb4[q][0];
		j=p->gcb4[q][1];
		k=p->gcb4[q][2];
		
		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
        mgc[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
		mgc[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
        mgc[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		mgc[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}

//--------------------------
//WALL2
	QGC4LOOP
	{
        i=p->gcb4[q][0];
		j=p->gcb4[q][1];
		k=p->gcb4[q][2];

		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
		if( mgc[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& mgc[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			mgc[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
        if( mgc[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]>1
		&& mgc[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]<10)
        {
			mgc[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
		if(mgc[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]>1
		&& mgc[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]<10)
        {
			mgc[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		if( mgc[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]>1
		&& mgc[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]<10)
        {
			mgc[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		if( mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]>1
		&& mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]<10)
        {
			mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]=p->gcextra4;
			++p->gcextra4;
        }

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		if(mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]>1
		&& mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]<10)
        {
			mgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]=p->gcextra4;
			++p->gcextra4;
        }
	}
}

void mgccart::gcdirfill(lexer* p, int gcb_count, int **gcb, int *mgc)
{

}

void mgccart::gcsidefill(lexer *p)
{
	
	fieldint5 side(p);
	int n;
	
	p->Iresize(p->gcside4,p->gcside4_size, p->gcb4_count); 
	p->gcside4_size=p->gcb4_count;
	
	GC4LOOP
	p->gcside4[n]=0;
	
	
	BLOOP
	side(i,j,k)=0;
	
	
	
	GC4LOOP
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
		if(p->gcb4[n][3]==1)
		side(i-1,j,k) +=1;
		
		if(p->gcb4[n][3]==4)
		side(i+1,j,k) +=1;
	}
	
	BLOOP
	if(side(i,j,k)>0)
	side(i,j,k)=10;
	
	GC4LOOP
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
		if(p->gcb4[n][3]==3)
		side(i,j-1,k) +=1;
		
		if(p->gcb4[n][3]==2)
		side(i,j+1,k) +=1;
	}
	
	BLOOP
	{
	if(side(i,j,k)>0 && side(i,j,k)<=2)
	side(i,j,k)=10;
	
	if(side(i,j,k)>10)
	side(i,j,k)=20;
	}
	
	GC4LOOP
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
		if(p->gcb4[n][3]==5)
		side(i,j,k-1) +=1;
		
		if(p->gcb4[n][3]==6)
		side(i,j,k+1) +=1;
	}
	
	BLOOP
	{
	if(side(i,j,k)>0 && side(i,j,k)<=2)
	side(i,j,k)=10;
	
	if(side(i,j,k)>10&& side(i,j,k)<=12)
	side(i,j,k)=20;
	
	if(side(i,j,k)>20)
	side(i,j,k)=30;	
	}
	
	BLOOP
	{
	if(side(i,j,k)==10)
	side(i,j,k)=1;
	
	if(side(i,j,k)==20)
	side(i,j,k)=2;
	
	if(side(i,j,k)==30)
	side(i,j,k)=3;	
	}
	
	
	
	GC4LOOP
	{
		i=p->gcb4[n][0];
		j=p->gcb4[n][1];
		k=p->gcb4[n][2];
		
		if(p->gcb4[n][3]==1)
		p->gcside4[n] = side(i-1,j,k);

		if(p->gcb4[n][3]==4)
		p->gcside4[n] = side(i+1,j,k);
		
		
		if(p->gcb4[n][3]==3)
		p->gcside4[n] = side(i,j-1,k);

		if(p->gcb4[n][3]==2)
		p->gcside4[n] = side(i,j+1,k);
		
		
		if(p->gcb4[n][3]==5)
		p->gcside4[n] = side(i,j,k-1);

		if(p->gcb4[n][3]==6)
		p->gcside4[n] = side(i,j,k+1);
	}
	
}
