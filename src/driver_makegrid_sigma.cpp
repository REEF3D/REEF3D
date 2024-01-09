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

#include"driver.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"mgc1.h"
#include"mgc2.h"
#include"mgc3.h"
#include"mgc4.h"
#include"mgcslice4.h"

void driver::makegrid_sigma(lexer *p, ghostcell *pgc)
{	
    int q;
    
// flag7
    p->Iarray(p->flag7,p->imax*p->jmax*(p->kmax+2));
    
    for(i=0;i<p->imax*p->jmax*(p->kmax+2);++i)
    p->flag7[i]=-10;
    
    // flag4
    BASELOOP
    p->flag7[FIJK]=p->flag4[IJK];
    
    // add solid structures
    BASELOOP
    if(p->flagslice4[IJ]<0)
    p->flag7[FIJK]=-10;
    
    // add solid structures
    BASELOOP
    if(p->flagslice4[IJ]<0)
    p->flag4[IJK]=-10;

    
    k=p->knoz;
    SLICEBASELOOP
    p->flag7[FIJK] = p->flag7[FIJKm1];
   
// gcx
    int xz,yz;
    int maxnum;
    
    p->Iarray(p->gcx7_count,6);
    p->Iarray(p->gcxco7_count,6);
    
// gcx7 allocation
    xz=yz=0;
    
    xz = (p->knox+2)*(p->knoz+2);
    yz = (p->knoy+2)*(p->knoz+2);
    
    maxnum = MAX(xz,yz);
    
    p->Iarray(p->gcx7,4,maxnum,3);
    
// gcxco7 allocation
    xz=yz=0;
    
    xz = 2*(p->knox+p->knoz+18);
    yz = 2*(p->knoy+p->knoz+18);
    
    maxnum = MAX(xz,yz);
    
    p->Iarray(p->gcxco7,4,maxnum,3);

// ---------------
// gcx7 
    //nb1
    q=0;
    if(p->nb1>=0)
    for(j=0;j<p->knoy;++j)
    for(k=0;k<p->knoz+1;++k)
    {
    p->gcx7[0][q][0] = 0;
    p->gcx7[0][q][1] = j;
    p->gcx7[0][q][2] = k;
    
    ++q;
    }
    p->gcx7_count[0]=q;
    

    //nb4
    q=0;
    if(p->nb4>=0)
    for(j=0;j<p->knoy;++j)
    for(k=0;k<p->knoz+1;++k)
    {
    p->gcx7[3][q][0] = p->knox-1;
    p->gcx7[3][q][1] = j;
    p->gcx7[3][q][2] = k;

    ++q;
    }
    p->gcx7_count[3]=q;
    
    //nb2
    q=0;
    if(p->nb2>=0)
    for(i=0;i<p->knox;++i)
    for(k=0;k<p->knoz+1;++k)
    {
    p->gcx7[1][q][0] = i;
    p->gcx7[1][q][1] = p->knoy-1;
    p->gcx7[1][q][2] = k;
    
    ++q;
    }
    p->gcx7_count[1]=q;
    
    //nb3
    q=0;
    if(p->nb3>=0)
    for(i=0;i<p->knox;++i)
    for(k=0;k<p->knoz+1;++k)
    {
    p->gcx7[2][q][0] = i;
    p->gcx7[2][q][1] = 0;
    p->gcx7[2][q][2] = k;
    
    ++q;
    }
    p->gcx7_count[2]=q;

// ---------------
// gcxco7 
    //nb1
    q=0;
    if(p->nb1>=0)
    {
        for(j=0;j<p->knoy;++j)
        {
        p->gcxco7[0][q][0] = 0;
        p->gcxco7[0][q][1] = j;
        p->gcxco7[0][q][2] = -1;
        ++q;
        }

        for(j=0;j<p->knoy;++j)
        {
        p->gcxco7[0][q][0] = 0;
        p->gcxco7[0][q][1] = j;
        p->gcxco7[0][q][2] = p->knoz+1;
        ++q;
        }
        
        for(k=-1;k<p->knoz+1;++k)
        {
        p->gcxco7[0][q][0] = 0;
        p->gcxco7[0][q][1] = -1;
        p->gcxco7[0][q][2] = k;
        ++q;
        }
        
        for(k=-1;k<p->knoz+1;++k)
        {
        p->gcxco7[0][q][0] = 0;
        p->gcxco7[0][q][1] = p->knoy;
        p->gcxco7[0][q][2] = k;
        ++q;
        }
    }
    
    p->gcxco7_count[0]=q;
    

    //nb2
    q=0;
    if(p->nb2>=0)
    {
        for(i=0;i<p->knox;++i)
        {
        p->gcxco7[1][q][0] = i;
        p->gcxco7[1][q][1] = p->knoy-1;
        p->gcxco7[1][q][2] = -1;
        ++q;
        }

        for(i=0;i<p->knox;++i)
        {
        p->gcxco7[1][q][0] = i;
        p->gcxco7[1][q][1] = p->knoy-1;
        p->gcxco7[1][q][2] = p->knoz+1;
        ++q;
        }
        
        for(k=-1;k<p->knoz+1;++k)
        {
        p->gcxco7[1][q][0] = -1;
        p->gcxco7[1][q][1] = p->knoy-1;
        p->gcxco7[1][q][2] = k;
        ++q;
        }
        
        for(k=-1;k<p->knoz+1;++k)
        {
        p->gcxco7[1][q][0] = p->knox;
        p->gcxco7[1][q][1] = p->knoy-1;
        p->gcxco7[1][q][2] = k;
        ++q;
        }
    }
    
    p->gcxco7_count[1]=q;
    

    //nb3
    q=0;
    if(p->nb3>=0)
    {
        for(i=0;i<p->knox;++i)
        {
        p->gcxco7[2][q][0] = i;
        p->gcxco7[2][q][1] = 0;
        p->gcxco7[2][q][2] = -1;
        ++q;
        }

        for(i=0;i<p->knox;++i)
        {
        p->gcxco7[2][q][0] = i;
        p->gcxco7[2][q][1] = 0;
        p->gcxco7[2][q][2] = p->knoz+1;
        ++q;
        }
        
        for(k=-1;k<p->knoz+1;++k)
        {
        p->gcxco7[2][q][0] = -1;
        p->gcxco7[2][q][1] = 0;
        p->gcxco7[2][q][2] = k;
        ++q;
        }
        
        for(k=-1;k<p->knoz+1;++k)
        {
        p->gcxco7[2][q][0] = p->knox;
        p->gcxco7[2][q][1] = 0;
        p->gcxco7[2][q][2] = k;
        ++q;
        }
    }
    
    p->gcxco7_count[2]=q;
    
    //cout<<p->mpirank<<" q3: "<<q<<endl;

    //nb4
    q=0;
    if(p->nb4>=0)
    {
        for(j=0;j<p->knoy;++j)
        {
        p->gcxco7[3][q][0] = p->knox-1;
        p->gcxco7[3][q][1] = j;
        p->gcxco7[3][q][2] = -1;
        ++q;
        }

        for(j=0;j<p->knoy;++j)
        {
        p->gcxco7[3][q][0] = p->knox-1;
        p->gcxco7[3][q][1] = j;
        p->gcxco7[3][q][2] = p->knoz+1;
        ++q;
        }
        
        for(k=-1;k<p->knoz+1;++k)
        {
        p->gcxco7[3][q][0] = p->knox-1;
        p->gcxco7[3][q][1] = -1;
        p->gcxco7[3][q][2] = k;
        ++q;
        }
        
        for(k=-1;k<p->knoz+1;++k)
        {
        p->gcxco7[3][q][0] = p->knox-1;
        p->gcxco7[3][q][1] = p->knoy;
        p->gcxco7[3][q][2] = k;
        ++q;
        }
    }
    p->gcxco7_count[3]=q;

    // -----
    pgc->flagx7(p,p->flag7);
    
    // ------

    p->vecsize(pgc);

    p->xcoormax=-1.0e9;
    p->xcoormin=1.0e9;
    p->ycoormax=-1.0e9;
    p->ycoormin=1.0e9;
    p->zcoormax=-1.0e9;
    p->zcoormin=1.0e9;

    LOOP
    {
        p->xcoormax = MAX(p->xcoormax,p->XN[IP1]);
        p->xcoormin = MIN(p->xcoormin,p->XN[IP]);
        p->ycoormax = MAX(p->ycoormax,p->YN[JP1]);
        p->ycoormin = MIN(p->ycoormin,p->YN[JP]);
        p->zcoormax = MAX(p->zcoormax,p->ZN[KP1]);
        p->zcoormin = MIN(p->zcoormin,p->ZN[KP]);
     }
     
     p->xcoormax=pgc->globalmax(p->xcoormax);
	 p->ycoormax=pgc->globalmax(p->ycoormax);
	 p->zcoormax=pgc->globalmax(p->zcoormax);
	 
	 p->xcoormin=pgc->globalmin(p->xcoormin);
	 p->ycoormin=pgc->globalmin(p->ycoormin);
	 p->zcoormin=pgc->globalmin(p->zcoormin);

}

void driver::makegrid2D_basic(lexer *p, ghostcell *pgc)
{
    // 2D
    pgc->gcsl_tpflag(p);    
    pgc->gcslflagx(p,p->flagslice4);
    
    mgcslice4 msl4(p);
    
    msl4.makemgc(p);
    msl4.gcb_seed(p);
    msl4.mgcsetup(p);
    msl4.fillmgc(p);
    msl4.gcdirfill(p);
    
    msl4.make_ggc(p);
    msl4.fill_ggc(p);
    
    pgc->gcsl_setbc4(p);
    pgc->gcsl_setbcio(p);
    
	pgc->dgcslini4(p); 
}
	
void driver::makegrid_sigma_cds(lexer *p, ghostcell *pgc)
{	
    p->flagini2D();
    p->gridini2D();	
    
    pgc->sizeS_update(p);
    pgc->gcxslupdate(p);
}
	
