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

#include"mgc3.h"
#include"lexer.h"

void mgc3::make_dgc(lexer* p)
{
    p->Iarray(hgc,imax*jmax*kmax);
    
    for(i=0;i<imax*jmax*kmax;++i)
    hgc[i]=0;
}

void mgc3::fill_dgc(lexer* p)
{
    int q, count;
    
    QGC3LOOP
	{
        i=p->gcb3[q][0];
        j=p->gcb3[q][1];
        k=p->gcb3[q][2];
		
		if(p->gcb3[q][3]==1)
		for(n=0;n<p->margin;++n)
        hgc[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb3[q][3]==4)
		for(n=0;n<p->margin;++n)
		hgc[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb3[q][3]==3)
		for(n=0;n<p->margin;++n)
        hgc[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;

		if(p->gcb3[q][3]==2)
		for(n=0;n<p->margin;++n)
		hgc[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;

		if(p->gcb3[q][3]==5)
		for(n=0;n<p->margin;++n)
		hgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;

		if(p->gcb3[q][3]==6)
		for(n=0;n<p->margin;++n)
		hgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}
    /*
    for(q=0;q<p->gcpara1_count;++q)
    {
    i=p->gcpara1[q][0];
    j=p->gcpara1[q][1];
    k=p->gcpara1[q][2];
        
        if(p->gcpara1[q][5]==1)
        for(n=0;n<p->margin;++n)
        hgc[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;
    }

    for(q=0;q<p->gcpara3_count;++q)
    {
    i=p->gcpara3[q][0];
    j=p->gcpara3[q][1];
    k=p->gcpara3[q][2];
        
        if(p->gcpara3[q][5]==1)
        for(n=0;n<p->margin;++n)
        hgc[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;
    }

	for(q=0;q<p->gcpara5_count;++q)
	{
    i=p->gcpara5[q][0];
    j=p->gcpara5[q][1];
    k=p->gcpara5[q][2];
        
        if(p->gcpara5[q][5]==1)
        for(n=0;n<p->margin;++n)
        hgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;
	}

	for(q=0;q<p->gcpara4_count;++q)
	{
    i=p->gcpara4[q][0];
    j=p->gcpara4[q][1];
    k=p->gcpara4[q][2];
        
        if(p->gcpara4[q][5]==1)
        for(n=0;n<p->margin;++n)
        hgc[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;
	}

	for(q=0;q<p->gcpara2_count;++q)
	{
    i=p->gcpara2[q][0];
    j=p->gcpara2[q][1];
    k=p->gcpara2[q][2];
        
        if(p->gcpara2[q][5]==1)
        for(n=0;n<p->margin;++n)
        hgc[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;
	}

	for(q=0;q<p->gcpara6_count;++q)
	{
	i=p->gcpara6[q][0];
    j=p->gcpara6[q][1];
    k=p->gcpara6[q][2];
        
        if(p->gcpara6[q][5]==1)
        for(n=0;n<p->margin;++n)
        hgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}
    */
    count=0;
    WLOOP
    {
        // i-j
        if(p->flag3[Im1Jm1K]<0 && p->flag3[Im1JK]<0 && p->flag3[IJm1K]<0 && hgc[Im1Jm1K]==0)
        if((i>=0 || p->nb1<0) && (j>=0 || p->nb3<0))
        ++count;
        
        if(p->flag3[Ip1Jm1K]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJm1K]<0 && hgc[Ip1Jm1K]==0)
        if((i<p->knox || p->nb4<0) && (j>=0 || p->nb3<0))
        ++count;
        
        if(p->flag3[Ip1Jp1K]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJp1K]<0 && hgc[Ip1Jp1K]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0))   
        ++count;
        
        if(p->flag3[Im1Jp1K]<0 && p->flag3[Im1JK]<0  && p->flag3[IJp1K]<0 && hgc[Im1Jp1K]==0)
        if((i>=0 || p->nb1<0) && (j<p->knoy || p->nb2<0)) 
        ++count;
        
        // i-k
        if(p->flag3[Im1JKm1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJKm1]<0 && hgc[Im1JKm1]==0)
        if((i>=0 || p->nb1<0) && (k>=0 || p->nb5<0))
        ++count;
        
        if(p->flag3[Ip1JKm1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJKm1]<0 && hgc[Ip1JKm1]==0)
        if((i<p->knox || p->nb4<0) && (k>=0 || p->nb5<0))
        ++count;
    
        if(p->flag3[Im1JKp1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJKp1]<0 && hgc[Im1JKp1]==0)
        if((i>=0 || p->nb1<0) && (k<p->knoz || p->nb6<0))
        ++count;
        
        if(p->flag3[Ip1JKp1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJKp1]<0 && hgc[Ip1JKp1]==0)
        if((i<p->knox || p->nb4<0) && (k<p->knoz || p->nb6<0))
        ++count;
        
        // j-k
        if(p->flag3[IJp1Kp1]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKp1]<0 && hgc[IJp1Kp1]==0)
        if((j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))   
        ++count;
        
        if(p->flag3[IJm1Kp1]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKp1]<0 && hgc[IJm1Kp1]==0)
        if((j>=0 || p->nb3<0) && (k<p->knoz || p->nb6<0)) 
        ++count;
        
        if(p->flag3[IJp1Km1]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKm1]<0 && hgc[IJp1Km1]==0)
        if((j<p->knoy || p->nb2<0) && (k>=0 || p->nb5<0))   
        ++count;
        
        if(p->flag3[IJm1Km1]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKm1]<0 && hgc[IJm1Km1]==0)
        if((j>=0 || p->nb3<0) && (k>=0 || p->nb5<0)) 
        ++count;
        
        
        // diag
        if(p->flag3[Im1Jm1Km1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKm1]<0 && hgc[Im1Jm1Km1]==0)
        if((i>=0 || p->nb1<0) && (j>=0 || p->nb3<0) && (k>=0 || p->nb5<0)) 
        ++count;
        
        if(p->flag3[Ip1Jm1Km1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKm1]<0 && hgc[Ip1Jm1Km1]==0)
        if((i<p->knox || p->nb4<0) && (j>=0 || p->nb3<0) && (k>=0 || p->nb5<0))  
        ++count;
        
        if(p->flag3[Ip1Jp1Km1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKm1]<0 && hgc[Ip1Jp1Km1]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0) && (k>=0 || p->nb5<0))  
        ++count;
        
        if(p->flag3[Im1Jp1Km1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKm1]<0 && hgc[Im1Jp1Km1]==0)
        if((i>=0 || p->nb1<0) && (j<p->knoy || p->nb2<0) && (k>=0 || p->nb5<0))  
        ++count;
        
        
        if(p->flag3[Im1Jm1Kp1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKp1]<0 && hgc[Im1Jm1Kp1]==0)
        if((i>=0 || p->nb1<0) && (j>=0 || p->nb3<0) && (k<p->knoz || p->nb6<0))   
        ++count;
        
        if(p->flag3[Ip1Jm1Kp1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKp1]<0 && hgc[Ip1Jm1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (j>=0 || p->nb3<0) && (k<p->knoz || p->nb6<0))  
        ++count;
        
        if(p->flag3[Ip1Jp1Kp1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKp1]<0 && hgc[Ip1Jp1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))
        ++count;
        
        if(p->flag3[Im1Jp1Kp1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKp1]<0 && hgc[Im1Jp1Kp1]==0)
        if((i>=0 || p->nb1<0) && (j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))   
        ++count;
    }
    
    p->Iresize(p->dgc3,p->dgc3_count,count,8,8);
    
    p->dgc3_count = count;
    
    //cout<<p->mpirank<<"  DGC3_count: "<<count<<endl;
    
    count=0;
    WLOOP
    {
        //i-j
        if(p->flag3[Im1Jm1K]<0 && p->flag3[Im1JK]<0 && p->flag3[IJm1K]<0 && hgc[Im1Jm1K]==0)
        if((i>=0 || p->nb1<0) && (j>=0 || p->nb3<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=-1;
        p->dgc3[count][4]=-1;
        p->dgc3[count][5]=0;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Ip1Jm1K]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJm1K]<0 && hgc[Ip1Jm1K]==0)
        if((i<p->knox || p->nb4<0) && (j>=0 || p->nb3<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=1;
        p->dgc3[count][4]=-1;
        p->dgc3[count][5]=0;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Ip1Jp1K]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJp1K]<0 && hgc[Ip1Jp1K]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=1;
        p->dgc3[count][4]=1;
        p->dgc3[count][5]=0;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Im1Jp1K]<0 && p->flag3[Im1JK]<0 && p->flag3[IJp1K]<0 && hgc[Im1Jp1K]==0)
        if((i>=0 || p->nb1<0) && (j<p->knoy || p->nb2<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=-1;
        p->dgc3[count][4]=1;
        p->dgc3[count][5]=0;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        // i-k
        if(p->flag3[Im1JKm1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJKm1]<0 && hgc[Im1JKm1]==0)
        if((i>=0 || p->nb1<0) && (k>=0 || p->nb5<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=-1;
        p->dgc3[count][4]=0;
        p->dgc3[count][5]=-1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Ip1JKm1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJKm1]<0 && hgc[Ip1JKm1]==0)
        if((i<p->knox || p->nb4<0) && (k>=0 || p->nb5<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=1;
        p->dgc3[count][4]=0;
        p->dgc3[count][5]=-1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
    
        if(p->flag3[Im1JKp1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJKp1]<0 && hgc[Im1JKp1]==0)
        if((i>=0 || p->nb1<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=-1;
        p->dgc3[count][4]=0;
        p->dgc3[count][5]=1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Ip1JKp1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJKp1]<0 && hgc[Ip1JKp1]==0)
        if((i<p->knox || p->nb4<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=1;
        p->dgc3[count][4]=0;
        p->dgc3[count][5]=1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        // j-k
        if(p->flag3[IJp1Kp1]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKp1]<0 && hgc[IJp1Kp1]==0)
        if((j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))   
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=0;
        p->dgc3[count][4]=1;
        p->dgc3[count][5]=1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[IJm1Kp1]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKp1]<0 && hgc[IJm1Kp1]==0)
        if((j>=0 || p->nb3<0) && (k<p->knoz || p->nb6<0)) 
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=0;
        p->dgc3[count][4]=-1;
        p->dgc3[count][5]=1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[IJp1Km1]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKm1]<0 && hgc[IJp1Km1]==0)
        if((j<p->knoy || p->nb2<0) && (k>=0 || p->nb5<0))   
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=0;
        p->dgc3[count][4]=1;
        p->dgc3[count][5]=-1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[IJm1Km1]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKm1]<0 && hgc[IJm1Km1]==0)
        if((j>=0 || p->nb3<0) && (k>=0 || p->nb5<0)) 
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=0;
        p->dgc3[count][4]=-1;
        p->dgc3[count][5]=-1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        
        // diag
        if(p->flag3[Im1Jm1Km1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKm1]<0 && hgc[Im1Jm1Km1]==0)
        if((i>=0 || p->nb1<0) && (j>=0 || p->nb3<0) && (k>=0 || p->nb5<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=-1;
        p->dgc3[count][4]=-1;
        p->dgc3[count][5]=-1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Ip1Jm1Km1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKm1]<0 && hgc[Ip1Jm1Km1]==0)
        if((i<p->knox || p->nb4<0) && (j>=0 || p->nb3<0) && (k>=0 || p->nb5<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=1;
        p->dgc3[count][4]=-1;
        p->dgc3[count][5]=-1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Ip1Jp1Km1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKm1]<0 && hgc[Ip1Jp1Km1]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0) && (k>=0 || p->nb5<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=1;
        p->dgc3[count][4]=1;
        p->dgc3[count][5]=-1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Im1Jp1Km1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKm1]<0 && hgc[Im1Jp1Km1]==0)
        if((i>=0 || p->nb1<0) && (j<p->knoy || p->nb2<0) && (k>=0 || p->nb5<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=-1;
        p->dgc3[count][4]=1;
        p->dgc3[count][5]=-1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        
        if(p->flag3[Im1Jm1Kp1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKp1]<0 && hgc[Im1Jm1Kp1]==0)
        if((i>=0 || p->nb1<0) && (j>=0 || p->nb3<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=-1;
        p->dgc3[count][4]=-1;
        p->dgc3[count][5]=1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Ip1Jm1Kp1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJm1K]<0 && p->flag3[IJKp1]<0 && hgc[Ip1Jm1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (j>=0 || p->nb3<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=1;
        p->dgc3[count][4]=-1;
        p->dgc3[count][5]=1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Ip1Jp1Kp1]<0 && p->flag3[Ip1JK]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKp1]<0 && hgc[Ip1Jp1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=1;
        p->dgc3[count][4]=1;
        p->dgc3[count][5]=1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
        
        if(p->flag3[Im1Jp1Kp1]<0 && p->flag3[Im1JK]<0 && p->flag3[IJp1K]<0 && p->flag3[IJKp1]<0 && hgc[Im1Jp1Kp1]==0)
        if((i>=0 || p->nb1<0) && (j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc3[count][0]=i;
        p->dgc3[count][1]=j;
        p->dgc3[count][2]=k;
        p->dgc3[count][3]=-1;
        p->dgc3[count][4]=1;
        p->dgc3[count][5]=1;
        p->dgc3[count][6]=1;
            
        ++count;
        }
    }
    
    
    p->del_Iarray(hgc,imax*jmax*kmax);
}










