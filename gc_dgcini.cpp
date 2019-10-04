/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::dgcini1(lexer* p)
{/*
    int count=0;
    int *num;
    fieldint5 dcount(p);
	
    ULOOP
    dcount(i,j,k)=0;
    
    count=0;
    GCB1
    {
    
    i=p->gcb1[n][0];
    j=p->gcb1[n][1];
    k=p->gcb1[n][2];
    
    ++dcount(i,j,k);
    }

    ULOOP
    if(dcount(i,j,k)==2)
    {
        if(p->flag1[UIm1JK]<0
        && p->flag1[UIp1JK]<0)
        dcount(i,j,k)=0;

        if(p->flag1[UIJm1K]<0
        && p->flag1[UIJp1K]<0)
        dcount(i,j,k)=0;

        if(p->flag1[UIJKm1]<0
        && p->flag1[UIJKp1]<0)
        dcount(i,j,k)=0;
    }

    count=1;
    ULOOP
    if(dcount(i,j,k)>1)
    ++count;
	

	p->Iresize(p->dgc1,p->dgc1_count,count,12,12);
	p->dgc1_count=count;
	
    p->Iarray(num,p->dgc1_count);

        count=0;
        ULOOP
        if(dcount(i,j,k)>1)
        {
        p->dgc1[count][0]=i;
        p->dgc1[count][1]=j;
        p->dgc1[count][2]=k;
        p->dgc1[count][3]=dcount(i,j,k);
        num[count]=dcount(i,j,k);
        ++count;
        }

        for(n=0;n<p->dgc1_count;++n)
        {
        i=p->dgc1[n][0];
        j=p->dgc1[n][1];
        k=p->dgc1[n][2];


            QGC1LOOP
            if(p->gcb1[q][0]==i && p->gcb1[q][1]==j && p->gcb1[q][2]==k && num[n]>0)
            {
            p->dgc1[n][4+p->dgc1[n][3]-num[n]] = p->gcb1[q][3];
            --num[n];
			
			if(num[n]<0)
			cout<<p->mpirank<<" NUM dgc 1 -!"<<endl;
            }
        }
		
        p->del_Iarray(num,p->dgc1_count);	*/
}

void ghostcell::dgcini2(lexer* p)
{/*
    int count=0;
    int *num;
    fieldint5 dcount(p);
	
    VLOOP
    dcount(i,j,k)=0;

    count=0;
    GC2LOOP
    {
    i=p->gcb2[n][0];
    j=p->gcb2[n][1];
    k=p->gcb2[n][2];
    ++dcount(i,j,k);
    }

    VLOOP
    if(dcount(i,j,k)==2)
    {
        if(p->flag2[VIm1JK]<0
        && p->flag2[VIp1JK]<0)
        dcount(i,j,k)=0;

        if(p->flag2[VIJm1K]<0
        && p->flag2[VIJp1K]<0)
        dcount(i,j,k)=0;

        if(p->flag2[VIJKm1]<0
        && p->flag2[VIJKp1]<0)
        dcount(i,j,k)=0;
    }

    count=0;
    VLOOP
    if(dcount(i,j,k)>1)
    ++count;
	

    p->Iresize(p->dgc2,p->dgc2_count,count,12,12);
	p->dgc2_count=count;
    p->Iarray(num,p->dgc2_count);

        count=0;
        VLOOP
        if(dcount(i,j,k)>1)
        {
        p->dgc2[count][0]=i;
        p->dgc2[count][1]=j;
        p->dgc2[count][2]=k;
        p->dgc2[count][3]=dcount(i,j,k);
        num[count]=dcount(i,j,k);
        ++count;
        }

        for(n=0;n<p->dgc2_count;++n)
        {
        i=p->dgc2[n][0];
        j=p->dgc2[n][1];
        k=p->dgc2[n][2];

            QGC2LOOP
            if(p->gcb2[q][0]==i && p->gcb2[q][1]==j && p->gcb2[q][2]==k && num[n]>0)
            {
            p->dgc2[n][4+p->dgc2[n][3]-num[n]] = p->gcb2[q][3];
            --num[n];
                
                if(num[n]<0)
                {
                cout<<p->mpirank<<" NUM dgc 2 -!  "<<num[n]<<endl;
                cout<<p->mpirank<<"i: "<<p->gcb2[q][0]<<" "<<i<<"  j: "<<p->gcb2[q][1]<<" "<<j<<"  k: "<<p->gcb2[q][2]<<" "<<k<<"   dir:"<<p->gcb2[q][3]<<" dcount: "<<p->dgc2[n][3]<<endl;
                }
            
            }
        }
		
        p->del_Iarray(num,p->dgc2_count);		*/
}

void ghostcell::dgcini3(lexer* p)
{/*
    int count=0;
    int *num;
    fieldint5 dcount(p);
	
    WLOOP
    dcount(i,j,k)=0;

    count=0;
    GC3LOOP
    {
    i=p->gcb3[n][0];
    j=p->gcb3[n][1];
    k=p->gcb3[n][2];
    ++dcount(i,j,k);
    }

    WLOOP
    if(dcount(i,j,k)==2)
    {
        if(p->flag3[WIm1JK]<0
        && p->flag3[WIp1JK]<0)
        dcount(i,j,k)=0;

        if(p->flag3[WIJm1K]<0
        && p->flag3[WIJp1K]<0)
        dcount(i,j,k)=0;

        if(p->flag3[WIJKm1]<0
        && p->flag3[WIJKp1]<0)
        dcount(i,j,k)=0;
    }

    count=0;
    WLOOP
    if(dcount(i,j,k)>1)
    ++count;

    p->Iresize(p->dgc3,p->dgc3_count,count,12,12);
	p->dgc3_count=count;
    p->Iarray(num,p->dgc3_count);

        count=0;
        WLOOP
        if(dcount(i,j,k)>1)
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=dcount(i,j,k);
        num[count]=dcount(i,j,k);
        ++count;
        }

        for(n=0;n<p->dgc3_count;++n)
        {
        i=p->dgc3[n][0];
        j=p->dgc3[n][1];
        k=p->dgc3[n][2];

            QGC3LOOP
            if(p->gcb3[q][0]==i && p->gcb3[q][1]==j && p->gcb3[q][2]==k && num[n]>0)
            {
            p->dgc3[n][4+p->dgc3[n][3]-num[n]] = p->gcb3[q][3];
            --num[n];
			
			if(num[n]<0)
			cout<<p->mpirank<<" NUM dgc 3 -!"<<endl;
            }
        }
		
        p->del_Iarray(num,p->dgc3_count);	*/	
}

void ghostcell::dgcini4(lexer* p)
{/*
    int count=0;
    int *num;
    fieldint5 dcount(p);
	
    LOOP
    dcount(i,j,k)=0;

    count=0;
    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
    ++dcount(i,j,k);
    }

    LOOP
    if(dcount(i,j,k)==2)
    {
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0
        && p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]<0)
        dcount(i,j,k)=0;

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]<0
        && p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]<0)
        dcount(i,j,k)=0;

        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]<0
        && p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]<0)
        dcount(i,j,k)=0;
    }

    count=0;
    LOOP
    if(dcount(i,j,k)>1)
    ++count;

    p->Iresize(p->dgc4,p->dgc4_count,count,12,12);
	p->dgc4_count=count;
    p->Iarray(num,p->dgc4_count);

        count=0;
        LOOP
        if(dcount(i,j,k)>1)
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=dcount(i,j,k);
        num[count]=dcount(i,j,k);
        ++count;
        }

        for(n=0;n<p->dgc4_count;++n)
        {
        i=p->dgc4[n][0];
        j=p->dgc4[n][1];
        k=p->dgc4[n][2];
		
		

            QGC4LOOP
            if(p->gcb4[q][0]==i && p->gcb4[q][1]==j && p->gcb4[q][2]==k && num[n]>0)
            {
            p->dgc4[n][4+p->dgc4[n][3]-num[n]] = p->gcb4[q][3];
            --num[n];
			
			if(num[n]<0)
			cout<<p->mpirank<<" NUM dgc 4 -!"<<endl;
            }
        }
		
        p->del_Iarray(num,p->dgc4_count);*/
}
