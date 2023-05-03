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

#include"ghostcell.h"
#include"lexer.h"
#include"slice.h"

void ghostcell::gcsl_setbc1(lexer *p)
{
    int cs,bc;
    
    GCSL1LOOP
    {
    i = p->gcbsl1[n][0];
    j = p->gcbsl1[n][1];
    cs = p->gcbsl1[n][3];
    bc = p->gcbsl1[n][4];
    
    if(cs==1 && bc==21 && i+p->origin_i==0)
    p->gcbsl1[n][4]=p->bcside1;
    
    if(cs==4 && bc==21 && i+p->origin_i==p->gknox-2)
    p->gcbsl1[n][4]=p->bcside4;
    
    if(cs==3 && bc==21 && j+p->origin_j==0)
    p->gcbsl1[n][4]=p->bcside3;
    
    if(cs==2 && bc==21 && j+p->origin_j==p->gknoy-1)
    p->gcbsl1[n][4]=p->bcside2;
    }    
}

void ghostcell::gcsl_setbc2(lexer *p)
{
    int cs,bc;
    
    GCSL2LOOP
    {
    i = p->gcbsl2[n][0];
    j = p->gcbsl2[n][1];
    cs = p->gcbsl2[n][3];
    bc = p->gcbsl2[n][4];
    
    if(cs==1 && bc==21 && i+p->origin_i==0)
    p->gcbsl2[n][4]=p->bcside1;
    
    if(cs==4 && bc==21 && i+p->origin_i==p->gknox-1)
    p->gcbsl2[n][4]=p->bcside4;
    
    if(cs==3 && bc==21 && j+p->origin_j==0)
    p->gcbsl2[n][4]=p->bcside3;
    
    if(cs==2 && bc==21 && j+p->origin_j==p->gknoy-2)
    p->gcbsl2[n][4]=p->bcside2;
    }    
}

void ghostcell::gcsl_setbc4(lexer *p)
{
    int cs,bc;
    
    GCSL4LOOP
    {
    i = p->gcbsl4[n][0];
    j = p->gcbsl4[n][1];
    cs = p->gcbsl4[n][3];
    bc = p->gcbsl4[n][4];
    
    if(cs==1 && bc==21 && i+p->origin_i==0)
    p->gcbsl4[n][4]=p->bcside1;
    
    if(cs==4 && bc==21 && i+p->origin_i==p->gknox-1)
    p->gcbsl4[n][4]=p->bcside4;
    
    if(cs==3 && bc==21 && j+p->origin_j==0)
    p->gcbsl4[n][4]=p->bcside3;
    
    if(cs==2 && bc==21 && j+p->origin_j==p->gknoy-1)
    p->gcbsl4[n][4]=p->bcside2;
    }    
}

void ghostcell::gcsl_setbcio(lexer *p)
{
    int cs,bc;

    p->gcslin_count=p->gcslout_count=0;
    GCSL4LOOP
    {
    i = p->gcbsl4[n][0];
    j = p->gcbsl4[n][1];
    cs = p->gcbsl4[n][3];
    bc = p->gcbsl4[n][4];
    
    if(bc==1 || bc==6)
    ++p->gcslin_count;

    if(bc==2 || bc==7 || bc==8)
    ++p->gcslout_count;
    }  
    
    //cout<<p->mpirank<<" "<<p->gcslin_count<<" "<<p->gcslout_count<<endl;
    
    
    p->Iarray(p->gcslin,p->gcslin_count,6);
    p->Iarray(p->gcslout,p->gcslout_count,6);
    
    p->Iarray(p->gcslawa1,p->gcslout_count,6);
    p->Iarray(p->gcslawa2,p->gcslout_count,6);


    int count1=0;
    int count2=0;
    
    GCSL4LOOP
    {
    i = p->gcbsl4[n][0];
    j = p->gcbsl4[n][1];
    cs = p->gcbsl4[n][3];
    bc = p->gcbsl4[n][4];
    
        if(bc==1 || bc==6)
        {
        p->gcslin[count1][0]=i;
        p->gcslin[count1][1]=j;
        p->gcslin[count1][3]=cs;
        p->gcslin[count1][4]=bc;
        p->gcslin[count1][5]=1;
        ++count1;
        }

        if(bc==2 || bc==7 || bc==8)
        {
        p->gcslout[count2][0]=i;
        p->gcslout[count2][1]=j;
        p->gcslout[count2][3]=cs;
        p->gcslout[count2][4]=bc;
        p->gcslout[count2][5]=1;
        ++count2;
        }
    }  
    
    count2=0;
    GCSL1LOOP
    {
    i = p->gcbsl1[n][0];
    j = p->gcbsl1[n][1];
    cs = p->gcbsl1[n][3];
    bc = p->gcbsl1[n][4];
    
        if(bc==2 || bc==7 || bc==8)
        {
        p->gcslawa1[count2][0]=i;
        p->gcslawa1[count2][1]=j;
        p->gcslawa1[count2][2]=cs;
        ++count2;
        }
    } 
    p->gcslawa1_count=count2;
    
    count2=0;
    GCSL2LOOP
    {
    i = p->gcbsl2[n][0];
    j = p->gcbsl2[n][1];
    cs = p->gcbsl2[n][3];
    bc = p->gcbsl2[n][4];
    
        if(bc==2 || bc==7 || bc==8)
        {
        p->gcslawa2[count2][0]=i;
        p->gcslawa2[count2][1]=j;
        p->gcslawa2[count2][2]=cs;
        ++count2;
        }
    } 
    p->gcslawa2_count=count2;
}
