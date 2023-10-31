/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include"patchBC_interface.h"

int iowave::iozonecheck(lexer *p, fdm*a)
{	
	int check=1;
	
	dg = distgen(p);
	db = distbeach(p);
	
	if(p->B98==2)
	if(dg<dist1 || db<dist2)
	check=0;

	return check;		
}

void iowave::gcio_update(lexer *p, fdm *a, ghostcell *pgc)
{
	
    int count1,count2;
	
	count1=0;
    count2=0;
    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        ++count1;

        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        ++count2;
    }
	
	//cout<<p->mpirank<<"  gcin_count: "<<p->gcin_count<<" count1: "<<count1<<"  gcout_count: "<<p->gcout_count<<" count2: "<<count2<<endl;
	
	p->Iresize(p->gcin,p->gcin_count, count1, 6, 6); 
	p->Iresize(p->gcout,p->gcout_count, count2, 6, 6); 

    count1=0;
    count2=0;
    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
        p->gcin[count1][0]=p->gcb4[n][0];
        p->gcin[count1][1]=p->gcb4[n][1];
        p->gcin[count1][2]=p->gcb4[n][2];
        p->gcin[count1][3]=p->gcb4[n][3];
        p->gcin[count1][5]=p->gcb4[n][5];
        ++count1;
        }

        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        {
        p->gcout[count2][0]=p->gcb4[n][0];
        p->gcout[count2][1]=p->gcb4[n][1];
        p->gcout[count2][2]=p->gcb4[n][2];
        p->gcout[count2][3]=p->gcb4[n][3];
        p->gcout[count2][5]=p->gcb4[n][5];
        ++count2;
        }
    }

    p->gcin_count=count1;
    p->gcout_count=count2;
    
     // 4a ---------------
	
	count1=0;
    count2=0;
    GC4ALOOP
    {
        if(p->gcb4a[n][4]==1)
        {
        p->gcin4a[count1][0]=p->gcb4a[n][0];
        p->gcin4a[count1][1]=p->gcb4a[n][1];
        p->gcin4a[count1][2]=p->gcb4a[n][2];
        p->gcin4a[count1][3]=p->gcb4a[n][3];
        p->gcin4a[count1][5]=p->gcb4a[n][5];
        ++count1;
        }

        if(p->gcb4a[n][4]==2)
        {
        p->gcout4a[count2][0]=p->gcb4a[n][0];
        p->gcout4a[count2][1]=p->gcb4a[n][1];
        p->gcout4a[count2][2]=p->gcb4a[n][2];
        p->gcout4a[count2][3]=p->gcb4a[n][3];
        p->gcout4a[count2][5]=p->gcb4a[n][5];
        ++count2;
        }
    }

    p->gcin4a_count=count1;
    p->gcout4a_count=count2;

    if(p->I10==1)
    velini(p,a,pgc);
	
	if(p->B98>=3)
	gen_ini(p,a,pgc);
	
	if(p->B99==3||p->B99==4||p->B99==5)
	awa_ini(p,a,pgc);
    
    
    // BC update
    MALOOP
    p->BC[IJK] = 0;
    
    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
        i = p->gcb4[n][0];
        j = p->gcb4[n][1];
        k = p->gcb4[n][2];
        
        // inflow
        if(p->gcb4[n][3]==1 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[Im1JK] = 1;
        
        if(p->gcb4[n][3]==4 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[Ip1JK] = 1;
        
        if(p->gcb4[n][3]==3 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[IJm1K] = 1;
        
        if(p->gcb4[n][3]==2 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[IJp1K] = 1;
        
        if(p->gcb4[n][3]==5 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[IJKm1] = 1;
        
        if(p->gcb4[n][3]==6 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[IJKp1] = 1;
        }

        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        {
        // outflow
        if(p->gcb4[n][3]==1 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[Im1JK] = 2;
        
        if(p->gcb4[n][3]==4 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[Ip1JK] = 2;
        
        if(p->gcb4[n][3]==3 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[IJm1K] = 2;
        
        if(p->gcb4[n][3]==2 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[IJp1K] = 2;
        
        if(p->gcb4[n][3]==5 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[IJKm1] = 2;
        
        if(p->gcb4[n][3]==6 && (p->gcb4[n][4]==1 || p->gcb4[n][4]==6))
        p->BC[IJKp1] = 2;
        }
    }
       
    for(int qq=0;qq<pBC->obj_count;++qq)
    for(n=0;n<pBC->patch[qq]->gcb_count;++n)
    {
    
    if(pBC->patch[qq]->gcb[n][3]==1)
    p->BC[Im1JK] = 1;
    
    if(pBC->patch[qq]->gcb[n][3]==4)
    p->BC[Ip1JK] = 1;
    
    if(pBC->patch[qq]->gcb[n][3]==3)
    p->BC[IJm1K] = 1;
    
    if(pBC->patch[qq]->gcb[n][3]==2)
    p->BC[IJp1K] = 1;
    
    if(pBC->patch[qq]->gcb[n][3]==5)
    p->BC[IJKm1] = 1;
    
    if(pBC->patch[qq]->gcb[n][3]==6)
    p->BC[IJKp1] = 1;
    }
    
}

void iowave::iogcb_update(lexer *p, fdm *a, ghostcell *pgc)
{
    int count1,count2;
	
	count1=0;
    count2=0;
    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        ++count1;

        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        ++count2;
    }
	
	
	p->Iresize(p->gcin,p->gcin_count, count1, 6, 6); 
	p->Iresize(p->gcout,p->gcout_count, count2, 6, 6); 

    count1=0;
    count2=0;
    GC4LOOP
    {
        if(p->gcb4[n][4]==1 || p->gcb4[n][4]==6)
        {
        p->gcin[count1][0]=p->gcb4[n][0];
        p->gcin[count1][1]=p->gcb4[n][1];
        p->gcin[count1][2]=p->gcb4[n][2];
        p->gcin[count1][3]=p->gcb4[n][3];
        p->gcin[count1][5]=p->gcb4[n][5];
        ++count1;
        }

        if(p->gcb4[n][4]==2 || p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        {
        p->gcout[count2][0]=p->gcb4[n][0];
        p->gcout[count2][1]=p->gcb4[n][1];
        p->gcout[count2][2]=p->gcb4[n][2];
        p->gcout[count2][3]=p->gcb4[n][3];
        p->gcout[count1][5]=p->gcb4[n][5];
        ++count2;
        }
    }

    p->gcin_count=count1;
    p->gcout_count=count2;
	
	if(p->I10==1)
    velini(p,a,pgc);
	
	if(p->B98==4)
	gen_ini(p,a,pgc);
	
	if(p->B99==3||p->B99==4||p->B99==5)
	awa_ini(p,a,pgc);
}

void iowave::awa_ini(lexer *p, fdm *a, ghostcell *pgc)
{
    int count1,count2,count3,count4;
	int flag,q;
	
	count=0;
    for(n=0;n<p->gcb4_count;++n)
    {
        if(p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        ++count;
    }
	
	p->Iresize(gcawa1, gcawa1_count,count, 4, 4); 
	p->Iresize(gcawa2, gcawa2_count,count, 4, 4); 
	p->Iresize(gcawa3, gcawa3_count,count, 4, 4); 
	p->Iresize(gcawa4, gcawa4_count,count, 4, 4); 
	gcawa1_count=count;
	gcawa2_count=count;
	gcawa3_count=count;
	gcawa4_count=count;
		
	// 1
    count1=0;
    for(n=0;n<p->gcb1_count;++n)
    {
        if(p->gcb1[n][4]==7 || p->gcb1[n][4]==8)
        {
			flag=1;
			for(q=0;q<count1;++q)
			if(gcawa1[q][0]==p->gcb1[n][0] && gcawa1[q][1]==p->gcb1[n][1] && gcawa1[q][2]==p->gcb1[n][3])
			flag=0;
			
				if(flag==1)
				{
				gcawa1[count1][0]=p->gcb1[n][0];
				gcawa1[count1][1]=p->gcb1[n][1];
				gcawa1[count1][2]=p->gcb1[n][3];
				++count1;
				}
        }
    }
	
	// 2
    count2=0;
    for(n=0;n<p->gcb2_count;++n)
    {
        if(p->gcb2[n][4]==7 || p->gcb2[n][4]==8)
        {
			flag=1;
			for(q=0;q<count2;++q)
			if(gcawa2[q][0]==p->gcb2[n][0] && gcawa2[q][1]==p->gcb2[n][1] && gcawa2[q][2]==p->gcb2[n][3])
			flag=0;
			
				if(flag==1)
				{
				gcawa2[count2][0]=p->gcb2[n][0];
				gcawa2[count2][1]=p->gcb2[n][1];
				gcawa2[count2][2]=p->gcb2[n][3];
				++count2;
				}
        }
    }
	
	// 3
    count3=0;
    for(n=0;n<p->gcb3_count;++n)
    {
        if(p->gcb3[n][4]==7 || p->gcb3[n][4]==8)
        {
			flag=1;
			for(q=0;q<count3;++q)
			if(gcawa3[q][0]==p->gcb3[n][0] && gcawa3[q][1]==p->gcb3[n][1] && gcawa3[q][2]==p->gcb3[n][3])
			flag=0;
			
				if(flag==1)
				{
				gcawa3[count3][0]=p->gcb3[n][0];
				gcawa3[count3][1]=p->gcb3[n][1];
				gcawa3[count3][2]=p->gcb3[n][3];
				++count3;
				}
        }
    }
	
	// 4
    count4=0;
    for(n=0;n<p->gcb4_count;++n)
    {
        if(p->gcb4[n][4]==7 || p->gcb4[n][4]==8)
        {
            
			flag=1;
			for(q=0;q<count4;++q)
			if(gcawa4[q][0]==p->gcb4[n][0] && gcawa4[q][1]==p->gcb4[n][1] && gcawa4[q][2]==p->gcb4[n][3])
			flag=0;
			
				if(flag==1)
				{
				gcawa4[count4][0]=p->gcb4[n][0];
				gcawa4[count4][1]=p->gcb4[n][1];
				gcawa4[count4][2]=p->gcb4[n][3];
				++count4;
				}
        }
    }
	
	p->Iresize(gcawa1, gcawa1_count,count1, 4, 4); 
	p->Iresize(gcawa2, gcawa2_count,count2, 4, 4); 
	p->Iresize(gcawa3, gcawa3_count,count3, 4, 4); 
	p->Iresize(gcawa4, gcawa4_count,count4, 4, 4); 
	
	gcawa1_count=count1;
	gcawa2_count=count2;
	gcawa3_count=count3;
	gcawa4_count=count4;
	
	//cout<<p->mpirank<<" GCAWA_COUNT: "<<gcawa4_count<<endl;	
}

void iowave::gen_ini(lexer *p, fdm *a, ghostcell *pgc)
{
    int count1,count2,count3,count4;
	int flag,q;
	
	count=0;
    for(n=0;n<p->gcb4_count;++n)
    {
        if(p->gcb4[n][4]==1||p->gcb4[n][4]==6)
        ++count;
    }
	
	p->Iresize(gcgen1, gcgen1_count,count, 4, 4); 
	p->Iresize(gcgen2, gcgen2_count,count, 4, 4); 
	p->Iresize(gcgen3, gcgen3_count,count, 4, 4); 
	p->Iresize(gcgen4, gcgen4_count,count, 4, 4); 
	gcgen1_count=count;
	gcgen2_count=count;
	gcgen3_count=count;
	gcgen4_count=count;
    
	// 1
    count1=0;
    for(n=0;n<p->gcb1_count;++n)
    {
        if(p->gcb1[n][4]==1||p->gcb1[n][4]==6)
        {
			flag=1;
			for(q=0;q<count1;++q)
			if(gcgen1[q][0]==p->gcb1[n][0] && gcgen1[q][1]==p->gcb1[n][1] && gcgen1[q][2]==p->gcb1[n][3])
			flag=0;
			
				if(flag==1)
				{
				gcgen1[count1][0]=p->gcb1[n][0];
				gcgen1[count1][1]=p->gcb1[n][1];
				gcgen1[count1][2]=p->gcb1[n][3];
				++count1;
				}
        }
    }
	
	// 2
    count2=0;
    for(n=0;n<p->gcb2_count;++n)
    {
        if(p->gcb2[n][4]==1||p->gcb2[n][4]==6)
        {
			flag=1;
			for(q=0;q<count2;++q)
			if(gcgen2[q][0]==p->gcb2[n][0] && gcgen2[q][1]==p->gcb2[n][1] && gcgen2[q][2]==p->gcb2[n][3])
			flag=0;
			
				if(flag==1)
				{
				gcgen2[count2][0]=p->gcb2[n][0];
				gcgen2[count2][1]=p->gcb2[n][1];
				gcgen2[count2][2]=p->gcb2[n][3];
				++count2;
				}
        }
    }
	
	// 3
    count3=0;
    for(n=0;n<p->gcb3_count;++n)
    {
        if(p->gcb3[n][4]==1||p->gcb3[n][4]==6)
        {
			flag=1;
			for(q=0;q<count3;++q)
			if(gcgen3[q][0]==p->gcb3[n][0] && gcgen3[q][1]==p->gcb3[n][1] && gcgen3[q][2]==p->gcb3[n][3])
			flag=0;
			
				if(flag==1)
				{
				gcgen3[count3][0]=p->gcb3[n][0];
				gcgen3[count3][1]=p->gcb3[n][1];
				gcgen3[count3][2]=p->gcb3[n][3];
				++count3;
				}
        }
    }
	
	// 4
    count4=0;
    for(n=0;n<p->gcb4_count;++n)
    {
        if(p->gcb4[n][4]==1||p->gcb4[n][4]==6)
        {
			flag=1;
			for(q=0;q<count4;++q)
			if(gcgen4[q][0]==p->gcb4[n][0] && gcgen4[q][1]==p->gcb4[n][1] && gcgen4[q][2]==p->gcb4[n][3])
			flag=0;
			
				if(flag==1)
				{
				gcgen4[count4][0]=p->gcb4[n][0];
				gcgen4[count4][1]=p->gcb4[n][1];
				gcgen4[count4][2]=p->gcb4[n][3];
				++count4;
				}
        }
    }
	//cout<<p->mpirank<<" GCGEN_COUNT: "<<gcgen4_count<<endl;	
}

void iowave::awa_update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void iowave::inflow_walldist(lexer *p, fdm *a, ghostcell *pgc, convection *pconvec, reini *preini, ioflow *pflow)
{
}

void iowave::veltimesave(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    pvrans->veltimesave(p,a,pgc);
    
}

void iowave::vrans_sed_update(lexer *p,fdm *a,ghostcell *pgc, vrans *pvrans)
{
    pvrans->sed_update(p,a,pgc);
}
