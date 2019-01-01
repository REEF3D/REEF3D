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
#include"sliceint5.h"

void ghostcell::dgcslini1(lexer* p)
{
    int count=0;
    int *num;
    sliceint5 dcount(p);
	
    SLICELOOP1
    dcount(i,j)=0;

    count=0;
    GCSL1LOOP
    {
    i=p->gcbsl1[n][0];
    j=p->gcbsl1[n][1];
    ++dcount(i,j);
    }

    SLICELOOP1
    if(dcount(i,j)==2)
    {
        if(p->flagslice1[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice1[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0)
        dcount(i,j)=0;

        if(p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice1[(i-p->imin)*p->jmax + (j-p->jmin+1)]<0)
        dcount(i,j)=0;
    }

    count=1;
    SLICELOOP1
    if(dcount(i,j)>1)
    ++count;
	

	p->Iresize(p->dgcsl1,p->dgcsl1_count,count,12,12);
	p->dgcsl1_count=count;
	
    p->Iarray(num,p->dgcsl1_count);

        count=0;
        SLICELOOP1
        if(dcount(i,j)>1)
        {
        p->dgcsl1[count][0]=i;
        p->dgcsl1[count][1]=j;
        p->dgcsl1[count][3]=dcount(i,j);
        num[count]=dcount(i,j);
        ++count;
        }

        for(n=0;n<p->dgcsl1_count;++n)
        {
        i=p->dgcsl1[n][0];
        j=p->dgcsl1[n][1];


            QGCSL1LOOP
            if(p->gcbsl1[q][0]==i && p->gcbsl1[q][1]==j && num[n]>0)
            {
            p->dgcsl1[n][4+p->dgcsl1[n][3]-num[n]] = p->gcbsl1[q][3];
            --num[n];
			
			if(num[n]<0)
			cout<<p->mpirank<<" NUM dgc 1 -!"<<endl;
            }
        }
		
        p->del_Iarray(num,p->dgcsl1_count);	

}

void ghostcell::dgcslini2(lexer* p)
{
    int count=0;
    int *num;
    sliceint5 dcount(p);
	
    SLICELOOP2
    dcount(i,j)=0;

    count=0;
    GCSL2LOOP
    {
    i=p->gcbsl2[n][0];
    j=p->gcbsl2[n][1];
    ++dcount(i,j);
    }

    SLICELOOP2
    if(dcount(i,j)==2)
    {
        if(p->flagslice2[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice2[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0)
        dcount(i,j)=0;

        if(p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice2[(i-p->imin)*p->jmax + (j-p->jmin+1)]<0)
        dcount(i,j)=0;
    }

    count=0;
    SLICELOOP2
    if(dcount(i,j)>1)
    ++count;
	

    p->Iresize(p->dgcsl2,p->dgcsl2_count,count,12,12);
	p->dgcsl2_count=count;
    p->Iarray(num,p->dgcsl2_count);

        count=0;
        SLICELOOP2
        if(dcount(i,j)>1)
        {
        p->dgcsl2[count][0]=i;
        p->dgcsl2[count][1]=j;
        p->dgcsl2[count][3]=dcount(i,j);
        num[count]=dcount(i,j);
        ++count;
        }

        for(n=0;n<p->dgcsl2_count;++n)
        {
        i=p->dgcsl2[n][0];
        j=p->dgcsl2[n][1];

            QGCSL2LOOP
            if(p->gcbsl2[q][0]==i && p->gcbsl2[q][1]==j && num[n]>0)
            {
            p->dgcsl2[n][4+p->dgcsl2[n][3]-num[n]] = p->gcbsl2[q][3];
            --num[n];
                
                if(num[n]<0)
                {
                cout<<p->mpirank<<" NUM dgc 2 -!  "<<num[n]<<endl;
                cout<<p->mpirank<<"i: "<<p->gcbsl2[q][0]<<" "<<i<<"  j: "<<p->gcbsl2[q][1]<<" "<<j<<"  k: "<<p->gcbsl2[q][2]<<" "<<k<<"   dir:"<<p->gcbsl2[q][3]<<" dcount: "<<p->dgcsl2[n][3]<<endl;
                }
            
            }
        }
		
        p->del_Iarray(num,p->dgcsl2_count);		
}

void ghostcell::dgcslini3(lexer* p)
{
    int count=0;
    int *num;
    sliceint5 dcount(p);
	
    SLICELOOP3
    dcount(i,j)=0;

    count=0;
    GCSL3LOOP
    {
    i=p->gcbsl3[n][0];
    j=p->gcbsl3[n][1];
    ++dcount(i,j);
    }

    SLICELOOP3
    if(dcount(i,j)==2)
    {
        if(p->flagslice3[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice3[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0)
        dcount(i,j)=0;

        if(p->flagslice3[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice3[(i-p->imin)*p->jmax + (j-p->jmin+1)]<0)
        dcount(i,j)=0;
    }

    count=0;
    SLICELOOP3
    if(dcount(i,j)>1)
    ++count;

    p->Iresize(p->dgcsl3,p->dgcsl3_count,count,12,12);
	p->dgcsl3_count=count;
    p->Iarray(num,p->dgcsl3_count);

        count=0;
        SLICELOOP3
        if(dcount(i,j)>1)
        {
        p->dgcsl3[count][0]=i;
        p->dgcsl3[count][1]=j;
        p->dgcsl3[count][3]=dcount(i,j);
        num[count]=dcount(i,j);
        ++count;
        }

        for(n=0;n<p->dgcsl3_count;++n)
        {
        i=p->dgcsl3[n][0];
        j=p->dgcsl3[n][1];
        k=p->dgcsl3[n][2];

            QGCSL3LOOP
            if(p->gcbsl3[q][0]==i && p->gcbsl3[q][1]==j && num[n]>0)
            {
            p->dgcsl3[n][4+p->dgcsl3[n][3]-num[n]] = p->gcbsl3[q][3];
            --num[n];
			
			if(num[n]<0)
			cout<<p->mpirank<<" NUM dgc 3 -!"<<endl;
            }
        }
		
        p->del_Iarray(num,p->dgcsl3_count);
		
}

void ghostcell::dgcslini4(lexer* p)
{
    int count=0;
    int *num;
    sliceint5 dcount(p);
	
    SLICELOOP4
    dcount(i,j)=0;

    count=0;
    GCSL4LOOP
    {
    i=p->gcbsl4[n][0];
    j=p->gcbsl4[n][1];
    ++dcount(i,j);
    }

    SLICELOOP4
    if(dcount(i,j)==2)
    {
        if(p->flagslice4[(i-p->imin-1)*p->jmax + (j-p->jmin)]<0
        && p->flagslice4[(i-p->imin+1)*p->jmax + (j-p->jmin)]<0)
        dcount(i,j)=0;

        if(p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin-1)]<0
        && p->flagslice4[(i-p->imin)*p->jmax + (j-p->jmin+1)]<0)
        dcount(i,j)=0;
    }

    count=0;
    SLICELOOP4
    if(dcount(i,j)>1)
    ++count;

    p->Iresize(p->dgcsl4,p->dgcsl4_count,count,12,12);
	p->dgcsl4_count=count;
    p->Iarray(num,p->dgcsl4_count);

        count=0;
        SLICELOOP4
        if(dcount(i,j)>1)
        {
        p->dgcsl4[count][0]=i;
        p->dgcsl4[count][1]=j;
        p->dgcsl4[count][3]=dcount(i,j);
        num[count]=dcount(i,j);
        ++count;
        }

        for(n=0;n<p->dgcsl4_count;++n)
        {
        i=p->dgcsl4[n][0];
        j=p->dgcsl4[n][1];
        

            QGCSL4LOOP
            if(p->gcbsl4[q][0]==i && p->gcbsl4[q][1]==j && num[n]>0)
            {
            p->dgcsl4[n][4+p->dgcsl4[n][3]-num[n]] = p->gcbsl4[q][3];
            --num[n];
			
			if(num[n]<0)
			cout<<p->mpirank<<" NUM dgc 4 -!"<<endl;
            }
        }
		
        p->del_Iarray(num,p->dgcsl4_count);
}
