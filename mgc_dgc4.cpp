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

#include"mgc4.h"
#include"lexer.h"

void mgc4::make_dgc(lexer* p)
{
    p->Iarray(hgc,imax*jmax*kmax);
    
    for(i=0;i<imax*jmax*kmax;++i)
    hgc[i]=0;
}

void mgc4::fill_dgc(lexer* p)
{
    int q, count;
    
    QGC4LOOP
	{
        i=p->gcb4[q][0];
        j=p->gcb4[q][1];
        k=p->gcb4[q][2];
		
		if(p->gcb4[q][3]==1)
		for(n=0;n<p->margin;++n)
        hgc[(i-imin-n-1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==4)
		for(n=0;n<p->margin;++n)
		hgc[(i-imin+n+1)*jmax*kmax + (j-jmin)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==3)
		for(n=0;n<p->margin;++n)
        hgc[(i-imin)*jmax*kmax + (j-jmin-n-1)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==2)
		for(n=0;n<p->margin;++n)
		hgc[(i-imin)*jmax*kmax + (j-jmin+n+1)*kmax + k-kmin]+=1;

		if(p->gcb4[q][3]==5)
		for(n=0;n<p->margin;++n)
		hgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin-n-1]+=1;

		if(p->gcb4[q][3]==6)
		for(n=0;n<p->margin;++n)
		hgc[(i-imin)*jmax*kmax + (j-jmin)*kmax + k-kmin+n+1]+=1;
	}
    
    count=0;
    LOOP
    {
        
        if(p->flag4[Im1Jm1K]<0 && hgc[Im1Jm1K]==0)
        if((i>=p->knox || p->nb1<0) && (j>=p->knoy || p->nb3<0))
        ++count;
        
        if(p->flag4[Ip1Jm1K]<0 && hgc[Ip1Jm1K]==0)
        if((i<p->knox || p->nb4<0) && (j>=p->knoy || p->nb3<0))
        ++count;
        
        if(p->flag4[Ip1Jp1K]<0 && hgc[Ip1Jp1K]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0))   
        ++count;
        
        if(p->flag4[Im1Jp1K]<0 && hgc[Im1Jp1K]==0)
        if((i>=p->knox || p->nb1<0) && (j<p->knoy || p->nb2<0)) 
        ++count;
        
        // k+-
        
        if(p->flag4[Im1JKm1]<0 && hgc[Im1JKm1]==0)
        if((i>=p->knox || p->nb1<0) && (k>=p->knoz || p->nb5<0))
        ++count;
        
        if(p->flag4[Ip1JKm1]<0 && hgc[Ip1JKm1]==0)
        if((i<p->knox || p->nb4<0) && (k>=p->knoz || p->nb5<0))
        ++count;
        
        if(p->flag4[IJp1Kp1]<0 && hgc[IJp1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (k<p->knoz || p->nb6<0))   
        ++count;
        
        if(p->flag4[IJp1Kp1]<0 && hgc[IJp1Kp1]==0)
        if((i>=p->knox || p->nb1<0) && (k<p->knoz || p->nb6<0)) 
        ++count;
        
        
        if(p->flag4[Im1Jm1Km1]<0 && hgc[Im1Jm1Km1]==0)
        if((i>=p->knox || p->nb1<0) && (j>=p->knoy || p->nb3<0) && (k>=p->knoz || p->nb5<0)) 
        ++count;
        
        if(p->flag4[Ip1Jm1Km1]<0 && hgc[Ip1Jm1Km1]==0)
        if((i<p->knox || p->nb4<0) && (j>=p->knoy || p->nb3<0) && (k>=p->knoz || p->nb5<0))  
        ++count;
        
        if(p->flag4[Ip1Jp1Km1]<0 && hgc[Ip1Jp1Km1]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0) && (k>=p->knoz || p->nb5<0))  
        ++count;
        
        if(p->flag4[Im1Jp1Km1]<0 && hgc[Im1Jp1Km1]==0)
        if((i>=p->knox || p->nb1<0) && (j<p->knoy || p->nb2<0) && (k>=p->knoz || p->nb5<0))  
        ++count;
        
        
        if(p->flag4[Im1Jm1Kp1]<0 && hgc[Im1Jm1Kp1]==0)
        if((i>=p->knox || p->nb1<0) && (j>=p->knoy || p->nb3<0) && (k<p->knoz || p->nb6<0))   
        ++count;
        
        if(p->flag4[Ip1Jm1Kp1]<0 && hgc[Ip1Jm1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (j>=p->knoy || p->nb3<0) && (k<p->knoz || p->nb6<0))  
        ++count;
        
        if(p->flag4[Ip1Jp1Kp1]<0 && hgc[Ip1Jp1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))
        ++count;
        
        if(p->flag4[Im1Jp1Kp1]<0 && hgc[Im1Jp1Kp1]==0)
        if((i>=p->knox || p->nb1<0) && (j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))   
        ++count;
    }
    
    p->Iresize(p->dgc4,p->dgc4_count,count,8,8);
    
    p->dgc4_count = count;
    
    count=0;
    LOOP
    {
        
        if(p->flag4[Im1Jm1K]<0 && hgc[Im1Jm1K]==0)
        if((i>=p->knox || p->nb1<0) && (j>=p->knoy || p->nb3<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=-1;
        p->dgc4[count][4]=-1;
        p->dgc4[count][5]=0;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Ip1Jm1K]<0 && hgc[Ip1Jm1K]==0)
        if((i<p->knox || p->nb4<0) && (j>=p->knoy || p->nb3<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=1;
        p->dgc4[count][4]=-1;
        p->dgc4[count][5]=0;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Ip1Jp1K]<0 && hgc[Ip1Jp1K]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=1;
        p->dgc4[count][4]=1;
        p->dgc4[count][5]=0;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Im1Jp1K]<0 && hgc[Im1Jp1K]==0)
        if((i>=p->knox || p->nb1<0) && (j<p->knoy || p->nb2<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=-1;
        p->dgc4[count][4]=1;
        p->dgc4[count][5]=0;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        //---
        if(p->flag4[Im1JKm1]<0 && hgc[Im1JKm1]==0)
        if((i>=p->knox || p->nb1<0) && (k>=p->knoz || p->nb5<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=-1;
        p->dgc4[count][4]=-1;
        p->dgc4[count][5]=0;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Ip1Jm1K]<0 && hgc[Ip1Jm1K]==0)
        if((i<p->knox || p->nb4<0) && (j>=p->knoy || p->nb3<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=1;
        p->dgc4[count][4]=-1;
        p->dgc4[count][5]=0;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Ip1Jp1K]<0 && hgc[Ip1Jp1K]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=1;
        p->dgc4[count][4]=1;
        p->dgc4[count][5]=0;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Im1Jp1K]<0 && hgc[Im1Jp1K]==0)
        if((i>=p->knox || p->nb1<0) && (j<p->knoy || p->nb2<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=-1;
        p->dgc4[count][4]=1;
        p->dgc4[count][5]=0;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        
        if(p->flag4[Im1Jm1Km1]<0 && hgc[Im1Jm1Km1]==0)
        if((i>=p->knox || p->nb1<0) && (j>=p->knoy || p->nb3<0) && (k>=p->knoz || p->nb5<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=-1;
        p->dgc4[count][4]=-1;
        p->dgc4[count][5]=-1;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Ip1Jm1Km1]<0 && hgc[Ip1Jm1Km1]==0)
        if((i<p->knox || p->nb4<0) && (j>=p->knoy || p->nb3<0) && (k>=p->knoz || p->nb5<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=1;
        p->dgc4[count][4]=-1;
        p->dgc4[count][5]=-1;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Ip1Jp1Km1]<0 && hgc[Ip1Jp1Km1]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0) && (k>=p->knoz || p->nb5<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=1;
        p->dgc4[count][4]=1;
        p->dgc4[count][5]=-1;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Im1Jp1Km1]<0 && hgc[Im1Jp1Km1]==0)
        if((i>=p->knox || p->nb1<0) && (j<p->knoy || p->nb2<0) && (k>=p->knoz || p->nb5<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=-1;
        p->dgc4[count][4]=1;
        p->dgc4[count][5]=-1;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        
        if(p->flag4[Im1Jm1Kp1]<0 && hgc[Im1Jm1Kp1]==0)
        if((i>=p->knox || p->nb1<0) && (j>=p->knoy || p->nb3<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=-1;
        p->dgc4[count][4]=-1;
        p->dgc4[count][5]=1;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Ip1Jm1Kp1]<0 && hgc[Ip1Jm1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (j>=p->knoy || p->nb3<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=1;
        p->dgc4[count][4]=-1;
        p->dgc4[count][5]=1;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Ip1Jp1Kp1]<0 && hgc[Ip1Jp1Kp1]==0)
        if((i<p->knox || p->nb4<0) && (j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=1;
        p->dgc4[count][4]=1;
        p->dgc4[count][5]=1;
        p->dgc4[count][6]=1;
            
        ++count;
        }
        
        if(p->flag4[Im1Jp1Kp1]<0 && hgc[Im1Jp1Kp1]==0)
        if((i>=p->knox || p->nb1<0) && (j<p->knoy || p->nb2<0) && (k<p->knoz || p->nb6<0))
        {
        p->dgc4[count][0]=i;
        p->dgc4[count][1]=j;
        p->dgc4[count][2]=k;
        p->dgc4[count][3]=-1;
        p->dgc4[count][4]=1;
        p->dgc4[count][5]=1;
        p->dgc4[count][6]=1;
            
        ++count;
        }
    }
    
    cout<<p->mpirank<<"  DGC_count: "<<count<<endl;
    
    
    p->del_Iarray(hgc,imax*jmax*kmax);
}










