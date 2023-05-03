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

#include"lexer.h"
#include"fdm.h"
#include"sliceint.h"
#include"field.h"
#include"ghostcell.h"

void ghostcell::start1V(lexer *p, double *f, int gcv)
{
    
    //  MPI Boundary Swap
    /*if(p->M10>0)
    {
    starttime=timer();
	gcparax(p,f,1);
	gcparacox(p,f,gcv);
    gcparacox(p,f,gcv);
	gcparacox(p,f,gcv);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    if(p->F10==1)
    nse1(p,a,f,gcv);*/
    
    // solid ghostcells
    /*starttime=timer();
	QQGC1LOOP
	gcdistro1V(p,f,p->gcb1[qq][0], p->gcb1[qq][1], p->gcb1[qq][2], p->gcb1[qq][5], p->gcd1[qq], gcv, p->gcb1[qq][4], p->gcb1[qq][3]);
	endtime=timer();
	p->gctime+=endtime-starttime;*/
    
    LOOP
    {  
        //if(p->B98!=3||bc(i-1,j)==0)
        if(p->flag4[Im1JK]<0)
        {
        f[Im1JK] = f[IJK];
        f[Im2JK] = f[IJK];
        f[Im3JK] = f[IJK];
        }
          
        //if(p->B99!=3||bc(i+1,j)==0)
        if(p->flag4[Ip1JK]<0)
        {
        f[Ip1JK] = f[IJK];
        f[Ip2JK] = f[IJK];
        f[Ip3JK] = f[IJK];
        }
        
        if(p->flag4[IJm1K]<0)
        {
        f[IJm1K] = f[IJK];
        f[IJm2K] = f[IJK];
        f[IJm3K] = f[IJK];
        }
        
        if(p->flag4[IJp1K]<0)
        {
        f[IJp1K] = f[IJK];
        f[IJp2K] = f[IJK];
        f[IJp3K] = f[IJK];
        }
        
        if(p->flag4[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
    }
    
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    
    // periodic ghostcells
    /*gcperiodicx(p,f,1);
    
    if(p->periodic1==1)
    gc_periodic(p, f, 1, 1);
    
    if(p->periodic2==1)
    gc_periodic(p, f, 1, 2);
    
    if(p->periodic3==1)
    gc_periodic(p, f, 1, 3);
    
    if(p->F10==1)
    nse1(p,a,f,gcv);
    
    if(p->Y40==1  || p->Y40==3)
    dgcpol1(p,f,gcv);
        
    if(p->Y40==2  || p->Y40==3)
    f.ggcpol(p);
    
    if(p->M10>0)
	gcparacox(p,f,gcv); */
    
    
    
    /*if(gcv==250)
    fivec(p,x,bc);
    
    if(gcv==150)
    fivec2D(p,x,bc);
    
    if(gcv==210)
    fivec_vel(p,x,bc);
    
    if(gcv==110)
    fivec2D_vel(p,x,bc);*/
  
}

void ghostcell::start2V(lexer *p, double *X, int gcv)
{
    
}

void ghostcell::start3V(lexer *p, double *f, int gcv)
{
    if(p->M10>0)
    {
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    }
    
    LOOP
    {  
        //if(p->B98!=3||bc(i-1,j)==0)
        if(p->flag4[Im1JK]<0)
        {
        f[Im1JK] = f[IJK];
        f[Im2JK] = f[IJK];
        f[Im3JK] = f[IJK];
        }
          
        //if(p->B99!=3||bc(i+1,j)==0)
        if(p->flag4[Ip1JK]<0)
        {
        f[Ip1JK] = f[IJK];
        f[Ip2JK] = f[IJK];
        f[Ip3JK] = f[IJK];
        }
        
        if(p->flag4[IJm1K]<0)
        {
        f[IJm1K] = f[IJK];
        f[IJm2K] = f[IJK];
        f[IJm3K] = f[IJK];
        }
        
        if(p->flag4[IJp1K]<0)
        {
        f[IJp1K] = f[IJK];
        f[IJp2K] = f[IJK];
        f[IJp3K] = f[IJK];
        }
        /*
        if(p->flag4[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }*/
    }
    
}

void ghostcell::start4V(lexer *p, double *f, int gcv)
{
    if(p->M10>0)
    {
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    }
    
    LOOP
    {  
        //if(p->B98!=3||bc(i-1,j)==0)
        if(p->flag4[Im1JK]<0)
        {
        f[Im1JK] = f[IJK];
        f[Im2JK] = f[IJK];
        f[Im3JK] = f[IJK];
        }
          
        //if(p->B99!=3||bc(i+1,j)==0)
        if(p->flag4[Ip1JK]<0)
        {
        f[Ip1JK] = f[IJK];
        f[Ip2JK] = f[IJK];
        f[Ip3JK] = f[IJK];
        }
        
        if(p->flag4[IJm1K]<0)
        {
        f[IJm1K] = f[IJK];
        f[IJm2K] = f[IJK];
        f[IJm3K] = f[IJK];
        }
        
        if(p->flag4[IJp1K]<0)
        {
        f[IJp1K] = f[IJK];
        f[IJp2K] = f[IJK];
        f[IJp3K] = f[IJK];
        }
        
        if(p->flag4[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
    }
    
    /*fivec_vel(p,x,bc);
    
    
    if(gcv==250)
    fivec(p,x,bc);
    
    if(gcv==150)
    fivec2D(p,x,bc);
    
    if(gcv==210)
    fivec_vel(p,x,bc);
    
    if(gcv==110)
    fivec2D_vel(p,x,bc);*/
}

void ghostcell::start7V(lexer *p, double *x, sliceint &bc, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,x,7);
    gcparax7co(p,x,7);
    gcparax7co(p,x,7);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    
    if(gcv==250)
    fivec(p,x,bc);
    
    if(gcv==150)
    fivec2D(p,x,bc);
    
    if(gcv==210)
    fivec_vel(p,x,bc);
    
    if(gcv==110)
    fivec2D_vel(p,x,bc);
}

void ghostcell::start7S(lexer *p, double *x, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,x,7);
    gcparax7co(p,x,7);
    gcparax7co(p,x,7);
    gcparax7co(p,x,7);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
}
