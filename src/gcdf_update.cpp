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

#include"ghostcell.h"
void gcbl_reset_all(lexer*);
#include"lexer.h"
#include"fdm.h"
#include"fieldint4.h"

void ghostcell::gcdf_update(lexer *p, fdm *a)
{
    gcbl_reset_all(p);

    double psi;
    
    // -----------------------------------------------------------
    // FLAG
    if(p->G505==1)
    {
    BASELOOP
    {
    if (p->j_dir==0)
    psi = -p->X41*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]);
	
    if (p->j_dir==1)
    psi = -p->X41*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]); 
    
    
    if( (a->solid(i,j,k)>=psi || p->solidread==0) && (a->topo(i,j,k)>=psi || p->toporead==0))
    p->flag4[IJK]=10;
    
    if( (a->solid(i,j,k)<psi && p->solidread==1) || (a->topo(i,j,k)<psi && p->toporead==1))
    p->flag4[IJK]=-10;
    }
    
    flagx(p,p->flag4);
    
    BASELOOP
    {
    p->flag1[IJK]=p->flag4[IJK];
    p->flag2[IJK]=p->flag4[IJK];
    p->flag3[IJK]=p->flag4[IJK];
    
    if(p->flag4[IJK]>0 && p->flag4[Ip1JK]<0)
    p->flag1[IJK]=-10;
    
    if(p->flag4[IJK]>0 && p->flag4[IJp1K]<0)
    p->flag2[IJK]=-10;
    
    if(p->flag4[IJK]>0 && p->flag4[IJKp1]<0)
    p->flag3[IJK]=-10;
    }
    
    flagx(p,p->flag1);
    flagx(p,p->flag2);
    flagx(p,p->flag3);
    }
    
    /*
    count=0;
    LOOP
    ++count;
    
    count=globalisum(count);

    if(p->mpirank==0)
    cout<<"number of active cells: "<<count<<endl;*/
    
    // -----------------------------------------------------------
    // FLAGSF
    BASELOOP
    {
    if((a->fb(i,j,k)>0.0 || p->X10==0) && (a->solid(i,j,k)>0.0 || p->solidread==0) && (a->topo(i,j,k)>0.0 || p->toporead==0))
    p->flagsf4[IJK]=1;
    
    if((a->fb(i,j,k)<0.0 && p->X10==1) || (a->solid(i,j,k)<0.0 && p->solidread==1) || (a->topo(i,j,k)<0.0 && p->toporead==1))
    p->flagsf4[IJK]=-1;
    }
    
    flagx(p,p->flagsf4);
    
    BASELOOP
    {
    p->flagsf1[IJK]=p->flagsf4[IJK];
    p->flagsf2[IJK]=p->flagsf4[IJK];
    p->flagsf3[IJK]=p->flagsf4[IJK];
    
    if(p->flagsf4[IJK]>0 && p->flagsf4[Ip1JK]<0)
    p->flagsf1[IJK]=-1;
    
    if(p->flagsf4[IJK]>0 && p->flagsf4[IJp1K]<0)
    p->flagsf2[IJK]=-1;
    
    if(p->flagsf4[IJK]>0 && p->flagsf4[IJKp1]<0)
    p->flagsf3[IJK]=-1;
    }
    
    flagx(p,p->flagsf1);
    flagx(p,p->flagsf2);
    flagx(p,p->flagsf3);
    
    
    // -----------------------------------------------------------
    // count gcdf entries
    count=0;
    
    
    // gcdf count
    BASELOOP
    if(p->flagsf4[IJK]>0)
    {
     
        if(p->flagsf4[Im1JK]<0)
        ++count;
        
        if(p->flagsf4[Ip1JK]<0)
        ++count;
        
        if(p->flagsf4[IJm1K]<0)
        ++count;
        
        if(p->flagsf4[IJp1K]<0)
        ++count;

        if(p->flagsf4[IJKm1]<0)
        ++count;
        
        if(p->flagsf4[IJKp1]<0)
        ++count;        
    }
    
    if(p->gcdf4_count!=count)
    {
    p->Iresize(p->gcdf4,p->gcdf4_count,count,6,6);
    
    p->gcdf4_count=count;
    }
    
    //cout<<p->mpirank<<" p->gcdf4_count: "<<p->gcdf4_count<<endl;
    
    // assign gcdf entries
    count=0;
    
    BASELOOP
    if(p->flagsf4[IJK]>0)
    {
        if(p->flagsf4[Im1JK]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=1;
        p->gcdf4[count][4]=48;
        ++count;
        }
        
        if(p->flagsf4[Ip1JK]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=4;
        p->gcdf4[count][4]=48;
        ++count;
        }
        
        if(p->flagsf4[IJm1K]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=3;
        p->gcdf4[count][4]=48;
        ++count;
        }
        
        if(p->flagsf4[IJp1K]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=2;
        p->gcdf4[count][4]=48;
        ++count;
        }

        if(p->flagsf4[IJKm1]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=5;
        p->gcdf4[count][4]=48;
        ++count;
        }
        
        if(p->flagsf4[IJKp1]<0)
        {
        p->gcdf4[count][0]=i;
        p->gcdf4[count][1]=j;
        p->gcdf4[count][2]=k;
        p->gcdf4[count][3]=6;
        p->gcdf4[count][4]=48;
        ++count;
        }       
    }
    
    fieldint4 cval(p);
    
    count=0;

    BASELOOP
	{
    cval(i,j,k)=count;
    
    ++count;
	}
    
    GC4LOOP
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];
	p->gcb4[n][5]=cval(i,j,k);
	}
    

    GCDF4LOOP
    {
        i=p->gcdf4[n][0];
        j=p->gcdf4[n][1];
        k=p->gcdf4[n][2];
        p->gcdf4[n][5]=cval(i,j,k);
    }

    gcdf_update_impl(p, p->flagsf1, p->gcdf1, p->gcdf1_count);
    gcdf_update_impl(p, p->flagsf2, p->gcdf2, p->gcdf2_count);
    gcdf_update_impl(p, p->flagsf3, p->gcdf3, p->gcdf3_count);
}

void ghostcell::gcdf_update_impl(lexer *p, const int *flagsf, int **&gcdf, int &gcdf_count)
{
    // -----------------------
    // flagsf

    count = 0;

    BASELOOP
    if(flagsf[IJK]>0)
    {
        if(flagsf[Im1JK]<0)
            ++count;

        if(flagsf[Ip1JK]<0)
            ++count;

        if(flagsf[IJm1K]<0)
            ++count;

        if(flagsf[IJp1K]<0)
            ++count;

        if(flagsf[IJKm1]<0)
            ++count;

        if(flagsf[IJKp1]<0)
            ++count;
    }

    if(gcdf_count!=count)
    {
        p->Iresize(gcdf,gcdf_count,count,6,6);

        gcdf_count=count;
    }
    
    //cout<<p->mpirank<<" p->gcdf3_count: "<<p->gcdf3_count<<endl;
    
    // assign gcdf entries
    count=0;
    
    BASELOOP
    if(flagsf[IJK]>0)
    {
        if(flagsf[Im1JK]<0)
        {
            gcdf[count][0]=i;
            gcdf[count][1]=j;
            gcdf[count][2]=k;
            gcdf[count][3]=1;
            gcdf[count][4]=48;
            ++count;
        }

        if(flagsf[Ip1JK]<0)
        {
            gcdf[count][0]=i;
            gcdf[count][1]=j;
            gcdf[count][2]=k;
            gcdf[count][3]=4;
            gcdf[count][4]=48;
            ++count;
        }

        if(flagsf[IJm1K]<0)
        {
            gcdf[count][0]=i;
            gcdf[count][1]=j;
            gcdf[count][2]=k;
            gcdf[count][3]=3;
            gcdf[count][4]=48;
            ++count;
        }

        if(flagsf[IJp1K]<0)
        {
            gcdf[count][0]=i;
            gcdf[count][1]=j;
            gcdf[count][2]=k;
            gcdf[count][3]=2;
            gcdf[count][4]=48;
            ++count;
        }

        if(flagsf[IJKm1]<0)
        {
            gcdf[count][0]=i;
            gcdf[count][1]=j;
            gcdf[count][2]=k;
            gcdf[count][3]=5;
            gcdf[count][4]=48;
            ++count;
        }

        if(flagsf[IJKp1]<0)
        {
            gcdf[count][0]=i;
            gcdf[count][1]=j;
            gcdf[count][2]=k;
            gcdf[count][3]=6;
            gcdf[count][4]=48;
            ++count;
        }
    }


}