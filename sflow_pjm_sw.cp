/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2017 Hans Bihs

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

#include"sflow_pjm_sw.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver2D.h"
#include"momentum.h"
#include"ioflow.h"

#define HXIJ (fabs(b->hx(i,j))>1.0e-20?b->hx(i,j):1.0e20)
#define HXIMJ (fabs(b->hx(i-1,j))>1.0e-20?b->hx(i-1,j):1.0e20)

#define HYIJ (fabs(b->hy(i,j))>1.0e-20?b->hy(i,j):1.0e20)
#define HYIJM (fabs(b->hy(i,j-1))>1.0e-20?b->hy(i,j-1):1.0e20)
 
sflow_pjm_sw::sflow_pjm_sw(lexer* p, fdm2D *b) : wsn(p), wbn(p)
{
    if(p->B76==0)
    gcval_press=40;  

    if(p->B76==1)
    gcval_press=41;

    if(p->B76==2)
    gcval_press=42;

    if(p->B76==3)
    gcval_press=43;
	
	if(p->B76==4) 
    gcval_press=44;
	
	if(p->B76==5) 
    gcval_press=45;
	
	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
}

sflow_pjm_sw::~sflow_pjm_sw()
{
}

void sflow_pjm_sw::start(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &uvel, slice &vvel, double alpha)
{
    if(p->mpirank==0)
    cout<<".";
	
    rhs(p,b,uvel,vvel,alpha);
	pgc->gcsl_start4(p,b->press,gcval_press);
	pgc->gcsl_start4(p,b->ws,50);
	pgc->gcsl_start4(p,b->wb,50);
	
    poisson(p,b);
	
        starttime=pgc->timer();

    psolv->start(p,b,pgc,b->press,b->xvec,b->rhsvec,4,gcval_press,p->N44);
	
        endtime=pgc->timer();

	pgc->gcsl_start4(p,b->press,gcval_press);
	
	ucorr(p,b,uvel,alpha);
	vcorr(p,b,vvel,alpha);
	wcorr(p,b,alpha,uvel,vvel);
	pgc->gcsl_start4(p,b->ws,50);
	pgc->gcsl_start4(p,b->wb,50);

    p->poissoniter=p->solveriter;

	p->poissontime=endtime-starttime;

	if(p->mpirank==0)
	cout<<"piter: "<<p->solveriter<<"  ptime: "<<setprecision(3)<<p->poissontime<<endl;
}


void sflow_pjm_sw::ucorr(lexer* p, fdm2D* b, slice& uvel,double alpha)
{	
	SLICELOOP1
	uvel(i,j) -= alpha*p->dt*((b->press(i+1,j)*(b->depth(i,j)+b->eta(i+1,j))-b->press(i,j)*(b->depth(i+1,j)+b->eta(i,j)))/(2.0*p->dx*HXIJ)
	
				+ 0.5*(b->press(i+1,j)+b->press(i,j))*((b->eta(i+1,j)-b->eta(i,j))-(b->depth(i+1,j)-b->depth(i,j)))
				/(2.0*p->dx*HXIJ));
}

void sflow_pjm_sw::vcorr(lexer* p, fdm2D* b, slice& vvel,double alpha)
{		
	SLICELOOP2
	vvel(i,j) -= alpha*p->dt*((b->press(i,j+1)*(b->depth(i,j)+b->eta(i,j+1))-b->press(i,j)*(b->depth(i,j+1)+b->eta(i,j)))/(2.0*p->dx*HYIJ)
	
				+ 0.5*(b->press(i,j+1)-b->press(i,j))*((b->eta(i,j+1)-b->eta(i,j))-(b->depth(i,j+1)-b->depth(i,j)))
				/(2.0*p->dx*HYIJ));
}

void sflow_pjm_sw::wcalc(lexer* p, fdm2D* b,double alpha, slice &uvel, slice &vvel)
{	
    SLICELOOP4
    wbn(i,j) = b->wb(i,j);
    
	/*
    SLICELOOP4
    b->wb(i,j) = -0.25*(uvel(i,j)+uvel(i-1,j))*(b->depth(i+1,j)-b->depth(i-1,j))/p->dx
                
                -0.25*(vvel(i,j)+vvel(i,j-1))*(b->depth(i,j+1)-b->depth(i,j-1))/p->dx;
	
	*/
	SLICELOOP4
    b->wb(i,j) = -MAX(0.0,uvel(i-1,j))*(b->depth(i,j)-b->depth(i-1,j))/p->dx
	
				 -MIN(0.0,uvel(i+1,j))*(b->depth(i+1,j)-b->depth(i,j))/p->dx
                
                -MAX(0.0,vvel(i,j-1))*(b->depth(i,j)-b->depth(i,j-1))/p->dx
	
				 -MIN(0.0,vvel(i,j+1))*(b->depth(i,j+1)-b->depth(i,j))/p->dx;
    
    SLICELOOP4
	b->ws(i,j) += - (b->wb(i,j)-wbn(i,j));
}

void sflow_pjm_sw::wcorr(lexer* p, fdm2D* b,double alpha, slice &uvel, slice &vvel)
{	
    SLICELOOP4
	b->ws(i,j) += p->dt*alpha*(2.0*b->press(i,j)/(b->hp(i,j)));
}

void sflow_pjm_sw::rhs(lexer *p, fdm2D* b, slice &u, slice &v, double alpha)
{
	double uvel,vvel,wvel;
	
	
    NSLICELOOP4
	b->rhsvec.V[n]=0.0;
	
	count=0;
    SLICELOOP4
    {
    b->rhsvec.V[count] =   -((u(i,j)*(b->hp(i,j)+b->dpx(i,j)) - u(i-1,j)*(b->hp(i,j)-b->dpx(i,j) ))
						    + (v(i,j)*(b->hp(i,j)+b->dpy(i,j)) - v(i,j-1)*(b->hp(i,j)-b->dpy(i,j))))/(alpha*p->dt*p->dx)
                           
                           -(b->ws(i,j)-b->wb(i,j))/(alpha*p->dt);
    ++count;
    }
}

void sflow_pjm_sw::upgrad(lexer*p, fdm2D* b)
{
	if(p->Z221==1)
    SLICELOOP1
	b->F(i,j)-= fabs(p->W22)*(b->eta(i+1,j)-b->eta(i,j))/(p->dx);

}

void sflow_pjm_sw::vpgrad(lexer*p, fdm2D* b)
{
	if(p->Z221==1)
    SLICELOOP2
	b->G(i,j)-= fabs(p->W22)*(b->eta(i,j+1)-b->eta(i,j))/(p->dx);
}

void sflow_pjm_sw::poisson(lexer*p, fdm2D* b)
{
    sqd = (1.0/(2.0*p->dx*p->dx));
    
	n=0;
    SLICELOOP4
	{
	b->M.p[n]  =  ((b->eta(i,j)+b->depth(i+1,j))/HXIJ) *(b->hp(i,j)+b->dpx(i,j))*sqd*p->x_dir
				+ ((b->eta(i,j)+b->depth(i-1,j))/HXIMJ)*(b->hp(i,j)-b->dpx(i,j))*sqd*p->x_dir
				+ ((b->eta(i,j)+b->depth(i,j+1))/HYIJ) *(b->hp(i,j)+b->dpy(i,j))*sqd*p->y_dir
				+ ((b->eta(i,j)+b->depth(i,j-1))/HYIJM)*(b->hp(i,j)-b->dpy(i,j))*sqd*p->y_dir
               +   2.0/b->hp(i,j);

   	b->M.n[n] = -((b->eta(i+1,j)+b->depth(i,j))/HXIJ) *(b->hp(i,j)+b->dpx(i,j))*sqd*p->x_dir;
	b->M.s[n] = -((b->eta(i-1,j)+b->depth(i,j))/HXIMJ)*(b->hp(i,j)-b->dpx(i,j))*sqd*p->x_dir;

	b->M.w[n] = -((b->eta(i,j+1)+b->depth(i,j))/HYIJ) *(b->hp(i,j)+b->dpy(i,j))*sqd*p->y_dir;
	b->M.e[n] = -((b->eta(i,j-1)+b->depth(i,j))/HYIJM)*(b->hp(i,j)-b->dpy(i,j))*sqd*p->y_dir;

	++n;
	}

    n=0;
	SLICELOOP4
	{
		if(p->flagslice4[Im1J]<0)
		{
		b->rhsvec.V[n] -= b->M.s[n]*b->press(i-1,j);
		b->M.s[n] = 0.0;
		}
		
		if(p->flagslice4[Ip1J]<0)
		{
		b->rhsvec.V[n] -= b->M.n[n]*b->press(i+1,j);
		b->M.n[n] = 0.0;
		}
		
		if(p->flagslice4[IJm1]<0)
		{
		b->rhsvec.V[n] -= b->M.e[n]*b->press(i,j-1);
		b->M.e[n] = 0.0;
		}
		
		if(p->flagslice4[IJp1]<0)
		{
		b->rhsvec.V[n] -= b->M.w[n]*b->press(i,j+1);
		b->M.w[n] = 0.0;
		}
		
	++n;
	}
}






