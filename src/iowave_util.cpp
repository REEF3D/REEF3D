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

#include"ioflow_gcio.h"
#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include"patchBC_interface.h"

void iowave::gcio_update(lexer *p, fdm *a, ghostcell *pgc)
{
    ioflow_gcio_lists(p,p->DF,nullptr);

    if(p->I10==1)
    velini(p,a,pgc);

	if(p->B98>=3)
	gen_ini(p,a,pgc);

	if(p->B99==3||p->B99==4||p->B99==5)
	awa_ini(p,a,pgc);

    ioflow_gcio_marks(p,p->DF,pBC);
}


void iowave::gcio_update_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    ioflow_gcio_lists(p,p->DF,nullptr);

    // the legacy inflow ghost cells (gcin) are written along x: an open y edge must be an outflow side in DIVEMesh
    if(p->open_ym==1 || p->open_yp==1)
    {
        int yin=0;
        for(n=0;n<p->gcin_count;++n)
        if((p->open_ym==1 && p->gcin[n][3]==3) || (p->open_yp==1 && p->gcin[n][3]==2))
        ++yin;
        
        if(pgc->globalisum(yin)>0)
        {
            if(p->mpirank==0)
            cout<<endl<<"!!! iowave: a Riemann / Flather edge at y- or y+ needs that side flagged as outflow in DIVEMesh (C 12 / C 13 2) !!!"<<endl<<endl;
            exit(1);
        }
    }
    
    

    //if(p->I10==1)
    //velini(p,a,pgc);
	    
    
    ioflow_gcio_marks(p,p->DF,pBC);
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
