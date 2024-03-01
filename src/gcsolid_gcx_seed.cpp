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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::gcxsd_seed(lexer *p, fdm *a)
{
	int count[6];
	
	for(q=0;q<6;++q)
	count[q]=0;
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
		
		if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID || p->flag4[IJK]==SOLID)
        ++count[0];
    }

    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];

    
		if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]==SOLID || p->flag4[IJK]==SOLID)
		++count[1];
    }

    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];

		if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==SOLID || p->flag4[IJK]==SOLID)
		++count[2];
    }

    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

		if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID || p->flag4[IJK]==SOLID)
		++count[3];

    }

    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
		
		if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==SOLID || p->flag4[IJK]==SOLID)
		++count[4];
    }

    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

		if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]==SOLID || p->flag4[IJK]==SOLID)
		++count[5];
    }
    
    p->Iresize(gcxsd,6,6,gcxsd_count,count,6,6); 
	
	for(q=0;q<6;++q)
	gcxsd_count[q]=count[q];
	
	//for(q=0;q<6;++q)
	//cout<<p->mpirank<<" GXCFB_COUNT_"<<q+1<<"  "<<count[q]<<endl;
    
 //--   
    for(q=0;q<6;++q)
	count[q]=0;
	
	for(n=0;n<p->gcpara1_count;++n)
    {
    i=p->gcpara1[n][0];
    j=p->gcpara1[n][1];
    k=p->gcpara1[n][2];
		
		if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID || p->flag4[IJK]==SOLID)
        {
			 gcxsd[0][count[0]][0]=i;
            gcxsd[0][count[0]][1]=j;
            gcxsd[0][count[0]][2]=k;
			++count[0];
		}
    }

    for(n=0;n<p->gcpara2_count;++n)
    {
    i=p->gcpara2[n][0];
    j=p->gcpara2[n][1];
    k=p->gcpara2[n][2];

    
		if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]==SOLID || p->flag4[IJK]==SOLID)
		{
			 gcxsd[1][count[1]][0]=i;
            gcxsd[1][count[1]][1]=j;
            gcxsd[1][count[1]][2]=k;
			++count[1];
		}
    }

    for(n=0;n<p->gcpara3_count;++n)
    {
    i=p->gcpara3[n][0];
    j=p->gcpara3[n][1];
    k=p->gcpara3[n][2];

		if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]==SOLID || p->flag4[IJK]==SOLID)
		{
			 gcxsd[2][count[2]][0]=i;
            gcxsd[2][count[2]][1]=j;
            gcxsd[2][count[2]][2]=k;
			++count[2];
		}
    }

    for(n=0;n<p->gcpara4_count;++n)
    {
    i=p->gcpara4[n][0];
    j=p->gcpara4[n][1];
    k=p->gcpara4[n][2];

		if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID || p->flag4[IJK]==SOLID)
		{
			 gcxsd[3][count[3]][0]=i;
            gcxsd[3][count[3]][1]=j;
            gcxsd[3][count[3]][2]=k;
			++count[3];
		}
    }

    for(n=0;n<p->gcpara5_count;++n)
    {
    i=p->gcpara5[n][0];
    j=p->gcpara5[n][1];
    k=p->gcpara5[n][2];
		
		if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]==SOLID || p->flag4[IJK]==SOLID)
		{
			 gcxsd[4][count[4]][0]=i;
            gcxsd[4][count[4]][1]=j;
            gcxsd[4][count[4]][2]=k;
			++count[4];
		}
    }

    for(n=0;n<p->gcpara6_count;++n)
    {
    i=p->gcpara6[n][0];
    j=p->gcpara6[n][1];
    k=p->gcpara6[n][2];

		if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]==SOLID || p->flag4[IJK]==SOLID)
		{
			 gcxsd[5][count[5]][0]=i;
            gcxsd[5][count[5]][1]=j;
            gcxsd[5][count[5]][2]=k;
			++count[5];
		}
    }
	
}
