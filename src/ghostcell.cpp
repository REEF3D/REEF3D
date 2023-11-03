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
#include"fdm.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include"density_f.h"

ghostcell::ghostcell(int& argc, char **argv,lexer *pp):size(15),tag1(1),tag2(2),tag3(3),tag4(4),tag5(5),tag6(6),eps(1.0e-10),
														gcx(1)
{  
    p=pp;
    
	MPI_Init(&argc,&argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD,&p->mpirank);
	MPI_Comm_size(MPI_COMM_WORLD,&p->mpi_size);	
    rank=p->mpirank;
	
	mpi_comm = MPI_COMM_WORLD;
}

ghostcell::~ghostcell()
{
}

void ghostcell::gcini(lexer* p)
{

    margin=p->margin;
	paramargin=p->margin;
    gamma=p->B29;
    Qi=p->W10;
    orderext=2;
    orderext2=2;
    orderpress=2;
    orderdir2=2;
    
    if(p->B23==1)
    orderdir=2;
    
    if(p->B23==2)
    orderdir=3;
	
	imin=p->imin;
    imax=p->imax;
    jmin=p->jmin;
    jmax=p->jmax;
    kmin=p->kmin;
    kmax=p->kmax;
    
    originx=p->originx;
    originy=p->originy;
    originz=p->originz;
	
    for(m=0;m<size;m++)
	{
	y[m]=0.0;
	x[m]=0.0;
	}

    tag=0;    
	

	gcx_count[0] = p->gcpara1_count*paramargin + p->gcparaco1_count*paramargin;
	gcx_count[1] = p->gcpara2_count*paramargin + p->gcparaco2_count*paramargin;
	gcx_count[2] = p->gcpara3_count*paramargin + p->gcparaco3_count*paramargin;
	gcx_count[3] = p->gcpara4_count*paramargin + p->gcparaco4_count*paramargin;
	gcx_count[4] = p->gcpara5_count*paramargin + p->gcparaco5_count*paramargin;
	gcx_count[5] = p->gcpara6_count*paramargin + p->gcparaco6_count*paramargin;
    

	p->Darray(send1,(p->gcpara1_count+p->flast)*paramargin + p->gcparaco1_count*paramargin);
	p->Darray(send2,(p->gcpara2_count+p->flast)*paramargin + p->gcparaco2_count*paramargin);
	p->Darray(send3,(p->gcpara3_count+p->flast)*paramargin + p->gcparaco3_count*paramargin);
	p->Darray(send4,(p->gcpara4_count+p->flast)*paramargin + p->gcparaco4_count*paramargin);
	p->Darray(send5,p->gcpara5_count*paramargin + p->gcparaco5_count*paramargin);
	p->Darray(send6,p->gcpara6_count*paramargin + p->gcparaco6_count*paramargin);
	
	p->Darray(recv1,(p->gcpara1_count+p->flast)*paramargin + p->gcparaco1_count*paramargin);
	p->Darray(recv2,(p->gcpara2_count+p->flast)*paramargin + p->gcparaco2_count*paramargin);
	p->Darray(recv3,(p->gcpara3_count+p->flast)*paramargin + p->gcparaco3_count*paramargin);
	p->Darray(recv4,(p->gcpara4_count+p->flast)*paramargin + p->gcparaco4_count*paramargin);
	p->Darray(recv5,p->gcpara5_count*paramargin + p->gcparaco5_count*paramargin);
	p->Darray(recv6,p->gcpara6_count*paramargin + p->gcparaco6_count*paramargin);
    
	p->Darray(send,6,gcx_count);
	p->Darray(recv,6,gcx_count);
	
	p->Iarray(isend1,(p->gcpara1_count+p->flast)*paramargin + p->gcparaco1_count*paramargin);
	p->Iarray(isend2,(p->gcpara2_count+p->flast)*paramargin + p->gcparaco2_count*paramargin);
	p->Iarray(isend3,(p->gcpara3_count+p->flast)*paramargin + p->gcparaco3_count*paramargin);
	p->Iarray(isend4,(p->gcpara4_count+p->flast)*paramargin + p->gcparaco4_count*paramargin);
	p->Iarray(isend5,p->gcpara5_count*paramargin + p->gcparaco5_count*paramargin);
	p->Iarray(isend6,p->gcpara6_count*paramargin + p->gcparaco6_count*paramargin);
	
	p->Iarray(irecv1,(p->gcpara1_count+p->flast)*paramargin + p->gcparaco1_count*paramargin);
	p->Iarray(irecv2,(p->gcpara2_count+p->flast)*paramargin + p->gcparaco2_count*paramargin);
	p->Iarray(irecv3,(p->gcpara3_count+p->flast)*paramargin + p->gcparaco3_count*paramargin);
	p->Iarray(irecv4,(p->gcpara4_count+p->flast)*paramargin + p->gcparaco4_count*paramargin);
	p->Iarray(irecv5,p->gcpara5_count*paramargin + p->gcparaco5_count*paramargin);
	p->Iarray(irecv6,p->gcpara6_count*paramargin + p->gcparaco6_count*paramargin);
	
	p->dgc1_count=1;
	p->dgc2_count=1;
	p->dgc3_count=1;
	p->dgc4_count=1;
	
	p->Iarray(p->dgc1,p->dgc1_count,8);
	p->Iarray(p->dgc2,p->dgc2_count,8);
	p->Iarray(p->dgc3,p->dgc3_count,8);
	p->Iarray(p->dgc4,p->dgc4_count,8);
    
    p->gcdf4_count=1;
    p->Iarray(p->gcdf4,p->gcdf4_count,6);

    if(p->B20==1)
    {
    gclabel_u=4;
    gclabel_v=4;
    gclabel_w=4;
    
    gclabel_utopo=4;
    gclabel_vtopo=4;
    gclabel_wtopo=4;
    }

    if(p->B20==2)
    {
    gclabel_u=5;
    gclabel_v=5;
    gclabel_w=5;
    
    gclabel_utopo=5;
    gclabel_vtopo=5;
    gclabel_wtopo=5;
    }

    if(p->B20==3)
    {
    gclabel_u=2;
    gclabel_v=2;
    gclabel_w=2;
    
    gclabel_utopo=2;
    gclabel_vtopo=2;
    gclabel_wtopo=2;
    }
    
    if(p->B20==4)
    {
    gclabel_u=5;
    gclabel_v=5;
    gclabel_w=5;
    
    gclabel_utopo=4;
    gclabel_vtopo=4;
    gclabel_wtopo=4;
    }
    

    gclabel_k=4;
    gclabel_e=4;

    gclabel_u_orth=1;
    gclabel_v_orth=1;
    gclabel_w_orth=1;
	gclabel_vel=5;    

    // for reflective BC
    if(p->B23==2)
    {
    gclabel_u=12;
    gclabel_v=12;
    gclabel_w=12;
    
    gclabel_u_orth=11;
    gclabel_v_orth=11;
    gclabel_w_orth=11;
    }
	

	gclabel_lsm=4;
	
	
	awa_lable=0;
	if(p->B99>=3)
	awa_lable=1;

	
	if(p->B75==1)
    {
    gclabel_u_out=4;
    gclabel_v_out=4;
    gclabel_w_out=4;
    }
	
	if(p->B75==2)
    {
    gclabel_u_out=6;
    gclabel_v_out=6;
    gclabel_w_out=6;
    }
    
    if(p->B75==3)
    {
    gclabel_u_out=0;
    gclabel_v_out=6;
    gclabel_w_out=6;
    }
	
	gclabel_outflow=1;
	if(p->B60==3||p->B60==4)
	gclabel_outflow=0;
    
    gclabel_u_in=1;
    gclabel_v_in=1;
    gclabel_w_in=1;
    gclabel_lsm_in=gclabel_lsm;
    
    if(p->I230>0 || p->B98>=3 || p->B60>0)
    {
    gclabel_u_in=0;
    gclabel_v_in=0;
    gclabel_w_in=0;
    gclabel_lsm_in=0;
    }
    
    if(p->B98>=3)
    {
    gclabel_u_in=0;
    gclabel_v_in=0;
    gclabel_w_in=0;
    gclabel_lsm_in=0;
    }
    
    // pressure bc labels
    gclabel_press=4;
    gclabel_press_in=gclabel_press;
    
    if(p->B76==2 || p->B76==3)
	gclabel_press_in=0;
    
    // pressure inflow
    pressin_lable=0;
	if(p->B76!=1)
	pressin_lable=1;
    
    // pressure outflow
    pressout_lable=0;
	if(p->B77==-1 || p->B77==2)
	pressout_lable=1;
    
    
    // sflow slip/no-slip
    if(p->A217==1 && p->A10==2)
    {
    gclabel_u=4;
    gclabel_v=4;
    }

    if(p->A217==2 && p->A10==2)
    {
    gclabel_u=5;
    gclabel_v=5;
    }
    
    
	

	margin=p->margin;
	
	nb[0] = p->nb1;
	nb[1] = p->nb2;
	nb[2] = p->nb3;
	nb[3] = p->nb4;
	nb[4] = p->nb5;
	nb[5] = p->nb6;
	
	stag[0] = 1;
	stag[1] = 2;
	stag[2] = 3;
	stag[3] = 4;
	stag[4] = 5;
	stag[5] = 6;
	
	rtag[0] = 4;
	rtag[1] = 3;
	rtag[2] = 2;
	rtag[3] = 1;
	rtag[4] = 6;
	rtag[5] = 5;
	
	p->Iarray(isend,6,p->gcpara_sum*9);
	p->Iarray(irecv,6,p->gcpara_sum*9);
	
	p->Darray(dsend,6,p->gcpara_sum*9);
	p->Darray(drecv,6,p->gcpara_sum*9);
	
	p->Darray(trecv,p->gcpara_sum*9*6);
	
	
	
	p->Iarray(gcbfb_count,6);
    p->Iarray(gcxfb_count,6);
    p->Iarray(gcbsd_count,6);
    p->Iarray(gcxsd_count,6);
	
	for(n=0;n<6;++n)
	gcbfb_count[n]=1;
	
	p->Iarray(gcbfb,6,gcbfb_count,6);
	
	
	for(n=0;n<6;++n)
	gcxfb_count[n]=1;
	
	p->Iarray(gcxfb,6,gcxfb_count,6);
    
    
    for(n=0;n<6;++n)
	gcxsd_count[n]=1;
	
	p->Iarray(gcxsd,6,gcxsd_count,6);
    
    
    for(n=0;n<6;++n)
	gcbsd_count[n]=1;
	
	p->Iarray(gcbsd,6,gcbsd_count,6);
	
	
	p->colnum = new int[p->M10+1];
	
	pdens = new density_f(p);
}

void ghostcell::fdm_update(fdm *aa)
{
    a=aa;
}

void ghostcell::fdm_fnpf_update(fdm_fnpf *cc)
{
    c=cc;
}

void ghostcell::fdm_nhf_update(fdm_nhf *dd)
{
    d=dd;
}

void ghostcell::final()
{
       MPI_Finalize();
}
