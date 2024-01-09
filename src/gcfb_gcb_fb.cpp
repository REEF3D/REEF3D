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

void ghostcell::gcfb_update_extra_gcb(lexer *p, fdm *a, field &f)
{
    gcfb_b_paraseed(p,a);
    //gcfb_x_paraseed(p,a);
    
    //gcparax_generic(a, p, f, gcxfb_count, gcxfb);
    gcb_generic_fbpress(p, f, gcbfb_count, gcbfb);
}

void ghostcell::gcfb_b_paraseed(lexer *p, fdm *a)
{
	int count[6];
	
	for(q=0;q<6;++q)
	count[q]=0;
	
	for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count[3];
    }
    
    for(q=0;q<p->gcpara2_count;++q)
    {
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]>0)
        ++count[2];
    }
    
    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]>0)
        ++count[1];
    }
    
    for(q=0;q<p->gcpara4_count;++q)
    {
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        ++count[0];   
    }
    
    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]>0)
        ++count[5];
    }
    
    for(q=0;q<p->gcpara6_count;++q)
    {
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]>0)
        ++count[4];
    }
	
	
	p->Iresize(gcbfb,6,6,gcbfb_count,count,6,6); 
	
	for(q=0;q<6;++q)
	gcbfb_count[q]=count[q];
	
	//for(q=0;q<6;++q)
	//cout<<p->mpirank<<" GCBFB_COUNT_"<<q+1<<"  "<<count[q]<<endl;
	
	
	for(q=0;q<6;++q)
	count[q]=0;

    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin-1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        {
		gcbfb[3][count[3]][0]=i-1;
		gcbfb[3][count[3]][1]=j;
		gcbfb[3][count[3]][2]=k;
		gcbfb[3][count[3]][3]=41;
		++count[3];
        }
    }
    
    for(q=0;q<p->gcpara2_count;++q)
    {
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin+1)*p->kmax + k-p->kmin]>0)
        {
        gcbfb[2][count[2]][0]=j;
		gcbfb[2][count[2]][1]=j+1;
		gcbfb[2][count[2]][2]=k;
		gcbfb[2][count[2]][3]=41;
		++count[2];
        }
    }
    
    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin-1)*p->kmax + k-p->kmin]>0)
        {
        gcbfb[1][count[1]][0]=i;
		gcbfb[1][count[1]][1]=j-1;
		gcbfb[1][count[1]][2]=k;
		gcbfb[1][count[1]][3]=41;
		++count[1];
        }
    }
    
    for(q=0;q<p->gcpara4_count;++q)
    {
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin+1)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin]>0)
        {
        gcbfb[0][count[0]][0]=i+1;
		gcbfb[0][count[0]][1]=j;
		gcbfb[0][count[0]][2]=k;
		gcbfb[0][count[0]][3]=41;
		++count[0];
        }
    }
    
    for(q=0;q<p->gcpara5_count;++q)
    {
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin-1]>0)
        {
        gcbfb[5][count[5]][0]=i;
		gcbfb[5][count[5]][1]=j;
		gcbfb[5][count[5]][2]=k-1;
		gcbfb[5][count[5]][3]=41;
		++count[5];
        }
    }
    
    for(q=0;q<p->gcpara6_count;++q)
    {
    i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->flag4[IJK]==FLT)
        if(p->flag4[(i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin+1]>0)
        {
        gcbfb[4][count[4]][0]=i;
		gcbfb[4][count[4]][1]=j;
		gcbfb[4][count[4]][2]=k+1;
		gcbfb[4][count[4]][3]=41;
		++count[4];
        }
    }
}

