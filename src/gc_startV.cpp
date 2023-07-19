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
#include"fdm_nhf.h"
#include"sliceint.h"
#include"field.h"
#include"ghostcell.h"

void ghostcell::start1V(lexer *p, double *f, int gcv)
{
    //  MPI Boundary Swap
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    
    ULOOP
    {  
        //if(p->B98!=3||bc(i-1,j)==0)
        if(p->flag1[Im1JK]<0 && gcv==10)
        {
        f[Im1JK] = 0.5*fabs(p->W22)*d->ETAs(i,j)*d->ETAs(i,j) + fabs(p->W22)*d->ETAs(i,j)*d->dfx(i,j);
        f[Im2JK] = 0.5*fabs(p->W22)*d->ETAs(i,j)*d->ETAs(i,j) + fabs(p->W22)*d->ETAs(i,j)*d->dfx(i,j);
        f[Im3JK] = 0.5*fabs(p->W22)*d->ETAs(i,j)*d->ETAs(i,j) + fabs(p->W22)*d->ETAs(i,j)*d->dfx(i,j);
        }
        
        if(p->flag1[Im1JK]<0 && gcv==14)
        {
        f[Im1JK] = 0.0;
        f[Im2JK] = 0.0;
        f[Im3JK] = 0.0;
        }
        
        if(p->flag1[Im1JK]<0 && gcv!=10 && gcv!=14)
        {
        f[Im1JK] = 0.0;
        f[Im2JK] = 0.0;
        f[Im3JK] = 0.0;
        }
          
        //if(p->B99!=3||bc(i+1,j)==0)
        if(p->flag1[Ip1JK]<0 && gcv==10)
        {
        f[Ip1JK] = 0.5*fabs(p->W22)*d->ETAn(i,j)*d->ETAn(i,j) + fabs(p->W22)*d->ETAn(i,j)*d->dfx(i,j);
        f[Ip2JK] = 0.5*fabs(p->W22)*d->ETAn(i,j)*d->ETAn(i,j) + fabs(p->W22)*d->ETAn(i,j)*d->dfx(i,j);
        f[Ip3JK] = 0.5*fabs(p->W22)*d->ETAn(i,j)*d->ETAn(i,j) + fabs(p->W22)*d->ETAn(i,j)*d->dfx(i,j);
        }
        
        if(p->flag1[Ip1JK]<0 && gcv==14)
        {
        f[Ip1JK] = 0.0;
        f[Ip2JK] = 0.0;
        f[Ip3JK] = 0.0;
        }
        
        
        if(p->flag1[Ip1JK]<0 && gcv!=10 && gcv!=14)
        {
        f[Ip1JK] = 0.0;
        f[Ip2JK] = 0.0;
        f[Ip3JK] = 0.0;
        }
        
        if(p->flag1[IJm1K]<0 && p->j_dir==1)
        {
        f[IJm1K] = f[IJK];
        f[IJm2K] = f[IJK];
        f[IJm3K] = f[IJK];
        }
        
        if(p->flag1[IJp1K]<0 && p->j_dir==1)
        {
        f[IJp1K] = f[IJK];
        f[IJp2K] = f[IJK];
        f[IJp3K] = f[IJK];
        }
        
        if(p->flag1[IJKm1]<0)
        {
        f[IJKm1] = f[IJK];
        f[IJKm2] = f[IJK];
        f[IJKm3] = f[IJK];
        }
        
        if(p->flag1[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
    }
}

void ghostcell::start2V(lexer *p, double *f, int gcv)
{
    //  MPI Boundary Swap
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    
    VLOOP
    {  
        //if(p->B98!=3||bc(i-1,j)==0)
        if(p->flag2[Im1JK]<0)
        {
        f[Im1JK] = f[IJK];
        f[Im2JK] = f[IJK];
        f[Im3JK] = f[IJK];
        }
          
        //if(p->B99!=3||bc(i+1,j)==0)
        if(p->flag2[Ip1JK]<0)
        {
        f[Ip1JK] = f[IJK];
        f[Ip2JK] = f[IJK];
        f[Ip3JK] = f[IJK];
        }
        
        if(p->flag2[IJm1K]<0 &&  gcv==11 && p->j_dir==1)
        {
        f[IJm1K] = 0.5*fabs(p->W22)*d->ETAe(i,j)*d->ETAe(i,j) + fabs(p->W22)*d->ETAe(i,j)*d->dfy(i,j);
        f[IJm2K] = 0.5*fabs(p->W22)*d->ETAe(i,j)*d->ETAe(i,j) + fabs(p->W22)*d->ETAe(i,j)*d->dfy(i,j);
        f[IJm3K] = 0.5*fabs(p->W22)*d->ETAe(i,j)*d->ETAe(i,j) + fabs(p->W22)*d->ETAe(i,j)*d->dfy(i,j);
        }
        
        if(p->flag2[IJm1K]<0 &&  gcv==14 && p->j_dir==1)
        {
        f[IJm1K] = 0.0;
        f[IJm2K] = 0.0;
        f[IJm3K] = 0.0;
        }
        
        if(p->flag2[IJm1K]<0 && gcv!=11 && gcv!=14 && p->j_dir==1)
        {
        f[IJm1K] = 0.0;
        f[IJm2K] = 0.0;
        f[IJm3K] = 0.0;
        }
        
        if(p->flag2[IJp1K]<0 &&  gcv==11 && p->j_dir==1)
        {
        f[IJp1K] = 0.5*fabs(p->W22)*d->ETAw(i,j)*d->ETAw(i,j) + fabs(p->W22)*d->ETAw(i,j)*d->dfy(i,j);
        f[IJp2K] = 0.5*fabs(p->W22)*d->ETAw(i,j)*d->ETAw(i,j) + fabs(p->W22)*d->ETAw(i,j)*d->dfy(i,j);
        f[IJp3K] = 0.5*fabs(p->W22)*d->ETAw(i,j)*d->ETAw(i,j) + fabs(p->W22)*d->ETAw(i,j)*d->dfy(i,j);
        }
        
        if(p->flag2[IJp1K]<0 &&  gcv==14 && p->j_dir==1)
        {
        f[IJp1K] = 0.0;
        f[IJp2K] = 0.0;
        f[IJp3K] = 0.0;
        }
        
        if(p->flag2[IJp1K]<0 && gcv!=11 && gcv!=14 && p->j_dir==1)
        {
        f[IJp1K] = 0.0;
        f[IJp2K] = 0.0;
        f[IJp3K] = 0.0;
        }
        
        if(p->flag2[IJKm1]<0 && p->j_dir==1)
        {
        f[IJKm1] = f[IJK];
        f[IJKm2] = f[IJK];
        f[IJKm3] = f[IJK];
        }
        
        if(p->flag2[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
    }
}

void ghostcell::start3V(lexer *p, double *f, int gcv)
{
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);

    
    WLOOP
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
     
    }
}

void ghostcell::start4V(lexer *p, double *f, int gcv)
{
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    
    LOOP
    {  
        //if(p->B98!=3||bc(i-1,j)==0)
    // xxxxxxx
        if(p->flag4[Im1JK]<0 && gcv==10)
        {
        f[Im1JK] = f[IJK];
        f[Im2JK] = f[IJK];
        f[Im3JK] = f[IJK];
        }
        
        if(p->flag4[Im1JK]<0 && gcv!=10)
        {
        f[Im1JK] = f[IJK];
        f[Im2JK] = f[IJK];
        f[Im3JK] = f[IJK];
        }
          
        //if(p->B99!=3||bc(i+1,j)==0)
        if(p->flag4[Ip1JK]<0 && gcv==10)
        {
        f[Ip1JK] = 0.0;
        f[Ip2JK] = 0.0;
        f[Ip3JK] = 0.0;
        }
        
        if(p->flag4[Ip1JK]<0 && gcv!=10)
        {
        f[Ip1JK] = f[IJK];
        f[Ip2JK] = f[IJK];
        f[Ip3JK] = f[IJK];
        }
        
        
    // yyyyy
        if(p->flag4[IJm1K]<0 && p->j_dir==1 && gcv==11)
        {
        f[IJm1K] = 0.0;
        f[IJm2K] = 0.0;
        f[IJm3K] = 0.0;
        }
        
        if(p->flag4[IJm1K]<0 && p->j_dir==1 && gcv!=11)
        {
        f[IJm1K] = f[IJK];
        f[IJm2K] = f[IJK];
        f[IJm3K] = f[IJK];
        }
        
        if(p->flag4[IJp1K]<0 && p->j_dir==1 && gcv==11)
        {
        f[IJp1K] = 0.0;
        f[IJp2K] = 0.0;
        f[IJp3K] = 0.0;
        }
        
        if(p->flag4[IJp1K]<0 && p->j_dir==1 && gcv!=11)
        {
        f[IJp1K] = f[IJK];
        f[IJp2K] = f[IJK];
        f[IJp3K] = f[IJK];
        }
        
    // zzzzz
        /*if(p->flag4[IJKp1]<0 && gcv==12)
        {
        f[IJKp1] = 0.0;
        f[IJKp2] = 0.0;
        f[IJKp3] = 0.0;
        }*/
        
        if(p->flag4[IJKp1]<0 && gcv!=12)
        {
        f[IJKp1] = 0.0;
        f[IJKp2] = 0.0;
        f[IJKp3] = 0.0;
        }
        
        /*if(p->flag4[IJKm1]<0 && gcv==12)
        {
        f[IJKm1] = 0.0;
        f[IJKm2] = 0.0;
        f[IJKm3] = 0.0;
        }*/
        
        if(p->flag4[IJKm1]<0 && gcv!=12)
        {
        f[IJKm1] = 0.0;
        f[IJKm2] = 0.0;
        f[IJKm3] = 0.0;
        }
    }
}


void ghostcell::start5V(lexer *p, double *f, int gcv)
{    
    starttime=timer();
	gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
	p->xtime+=endtime-starttime;

}

void ghostcell::start7V(lexer *p, double *f, sliceint &bc, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    if(gcv==250)
    fivec(p,f,bc);
    
    if(gcv==150)
    fivec2D(p,f,bc);
    
    if(gcv==210)
    fivec_vel(p,f,bc);
    
    if(gcv==110)
    fivec2D_vel(p,f,bc);
}

void ghostcell::start7S(lexer *p, double *f, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
}


/*7
void ghostcell::start1V(lexer *p, double *f, int gcv)
{
    //  MPI Boundary Swap
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    
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
        
        if(p->flag4[IJKm1]<0)
        {
        f[IJKm1] = 0.0;
        f[IJKm2] = 0.0;
        f[IJKm3] = 0.0;
        }
        
        if(p->flag4[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
        
    }
}

void ghostcell::start2V(lexer *p, double *f, int gcv)
{
    //  MPI Boundary Swap
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    
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
        
        if(p->flag4[IJKm1]<0)
        {
        f[IJKm1] = 0.0;
        f[IJKm2] = 0.0;
        f[IJKm3] = 0.0;
        }
        
        if(p->flag4[IJKp1]<0)
        {
        f[IJKp1] = f[IJK];
        f[IJKp2] = f[IJK];
        f[IJKp3] = f[IJK];
        }
        
    }
    
}

void ghostcell::start3V(lexer *p, double *f, int gcv)
{
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);

    
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
    }
    
}

void ghostcell::start4V(lexer *p, double *f, int gcv)
{
    gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    
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

}

void ghostcell::start5V(lexer *p, double *f, int gcv)
{    
    starttime=timer();
	gcparaxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
    gcparacoxV(p, f, gcv);
	p->xtime+=endtime-starttime;
}

void ghostcell::start7V(lexer *p, double *f, sliceint &bc, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
    
    if(gcv==250)
    fivec(p,f,bc);
    
    if(gcv==150)
    fivec2D(p,f,bc);
    
    if(gcv==210)
    fivec_vel(p,f,bc);
    
    if(gcv==110)
    fivec2D_vel(p,f,bc);
}

void ghostcell::start7S(lexer *p, double *f, int gcv)
{
    if(p->M10>0)
    {
    starttime=timer();
	gcparax7(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
    gcparax7co(p,f,7);
	endtime=timer();
	p->xtime+=endtime-starttime;
    }
}*/