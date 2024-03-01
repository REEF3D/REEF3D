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

void ghostcell::gcbsd_seed(lexer *p, fdm *a)
{
    int count[6];
	
	for(q=0;q<6;++q)
	count[q]=0;
	
	for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count[3];
        
        /*if(p->flag4[IJK]>0)
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID)
        ++count[0];*/
    }
    
    for(q=0;q<p->gcpara2_count;++q)
    {
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]>0)
        ++count[2];
    }
    
    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]>0)
        ++count[1];
    }
    
    for(q=0;q<p->gcpara4_count;++q)
    {
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count[0];  
        
        /*if(p->flag4[IJK]>0)
        if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID)
        ++count[3];  */
    }
    
    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]>0)
        ++count[5];
    }
    
    for(q=0;q<p->gcpara6_count;++q)
    {
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]>0)
        ++count[4];
    }
	
	
	p->Iresize(gcbsd,6,6,gcbsd_count,count,6,6); 
	
	for(q=0;q<6;++q)
	gcbsd_count[q]=count[q];
	
	//for(q=0;q<6;++q)
	//cout<<p->mpirank<<" GCBSD_COUNT_"<<q+1<<"  "<<count[q]<<endl;
	
	
	for(q=0;q<6;++q)
	count[q]=0;

    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        {
		gcbsd[3][count[3]][0]=i-1;
		gcbsd[3][count[3]][1]=j;
		gcbsd[3][count[3]][2]=k;
		gcbsd[3][count[3]][3]=21;
		++count[3];
        }
        /*
        if(p->flag4[IJK]>0)
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID)
        {
        gcbsd[0][count[0]][0]=i;
		gcbsd[0][count[0]][1]=j;
		gcbsd[0][count[0]][2]=k;
		gcbsd[0][count[0]][3]=21;
		++count[0];
        }*/
    }
    
    for(q=0;q<p->gcpara2_count;++q)
    {
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]>0)
        {
        gcbsd[2][count[2]][0]=j;
		gcbsd[2][count[2]][1]=j+1;
		gcbsd[2][count[2]][2]=k;
		gcbsd[2][count[2]][3]=21;
		++count[2];
        }
    }
    
    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]>0)
        {
        gcbsd[1][count[1]][0]=i;
		gcbsd[1][count[1]][1]=j-1;
		gcbsd[1][count[1]][2]=k;
		gcbsd[1][count[1]][3]=21;
		++count[1];
        }
    }
    
    for(q=0;q<p->gcpara4_count;++q)
    {
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        {
        gcbsd[0][count[0]][0]=i+1;
		gcbsd[0][count[0]][1]=j;
		gcbsd[0][count[0]][2]=k;
		gcbsd[0][count[0]][3]=21;
		++count[0];
        }
        
        /*
        if(p->flag4[IJK]>0)
        if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]==SOLID)
        {
		gcbsd[3][count[3]][0]=i;
		gcbsd[3][count[3]][1]=j;
		gcbsd[3][count[3]][2]=k;
		gcbsd[3][count[3]][3]=21;
		++count[3];
        }*/
    }
    
    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]>0)
        {
        gcbsd[5][count[5]][0]=i;
		gcbsd[5][count[5]][1]=j;
		gcbsd[5][count[5]][2]=k-1;
		gcbsd[5][count[5]][3]=21;
		++count[5];
        }
    }
    
    for(q=0;q<p->gcpara6_count;++q)
    {
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->flag4[IJK]==SOLID)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]>0)
        {
        gcbsd[4][count[4]][0]=i;
		gcbsd[4][count[4]][1]=j;
		gcbsd[4][count[4]][2]=k+1;
		gcbsd[4][count[4]][3]=21;
		++count[4];
        }
    }
    
    
}