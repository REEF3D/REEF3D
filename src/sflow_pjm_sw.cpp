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

#include"sflow_pjm_sw.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver2D.h"
#include"momentum.h"
#include"ioflow.h"
#include"sflow_weno_hj.h"
#include"patchBC_interface.h"

#define HXIJ (fabs(b->hx(i,j))>1.0e-20?b->hx(i,j):1.0e20)
#define HXIMJ (fabs(b->hx(i-1,j))>1.0e-20?b->hx(i-1,j):1.0e20)

#define HXP (0.5*(HXIJ + HXIMJ))
#define HYP (0.5*(HYIJ + HYIJM))

#define HYIJ (fabs(b->hy(i,j))>1.0e-20?b->hy(i,j):1.0e20)
#define HYIJM (fabs(b->hy(i,j-1))>1.0e-20?b->hy(i,j-1):1.0e20)

#define HP (fabs(b->hp(i,j))>1.0e-20?b->hp(i,j):1.0e20)
#define HPI (fabs(b->hp(i+1,j))>1.0e-20?b->hp(i+1,j):1.0e20)
#define HPJ (fabs(b->hp(i,j+1))>1.0e-20?b->hp(i,j+1):1.0e20)

#define HPX (0.5*(HP + HPI))
#define HPY (0.5*(HP + HPJ))
 
sflow_pjm_sw::sflow_pjm_sw(lexer* p, fdm2D *b, patchBC_interface *ppBC) : wb(p), wsn(p), wbn(p)
{
    pBC = ppBC;
    
    gcval_press=40;  

	gcval_u=10;
	gcval_v=11;
	gcval_w=12;
	
    
    
    wd_criterion=0.00005;
    
    if(p->A244==1)
    wd_criterion=p->A244_val;
    
    if(p->A245==1)
    wd_criterion=p->A245_val*p->DXM;
	
}

sflow_pjm_sw::~sflow_pjm_sw()
{
}

void sflow_pjm_sw::start(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, ioflow *pflow, slice &P, slice &Q, slice &Pn, slice &Qn, slice &ws, slice &eta, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
	
	starttime=pgc->timer();
	
    pgc->gcsl_start4(p,ws,12);
	pgc->gcsl_start4(p,wb,12);

    rhs(p,b,P,Q,ws,alpha);
	pgc->gcsl_start4(p,b->press,gcval_press);
	
    poisson(p,b);
	
        solvtime=pgc->timer();

    psolv->start(p,pgc,b->press,b->M,b->xvec,b->rhsvec,4);
	
        p->poissontime=pgc->timer()-solvtime;
    
    SLICELOOP4
    if(b->hp(i,j)<wd_criterion)
    b->press(i,j)=0.0;
  
    pflow->pm_relax(p,pgc,b->press);
	pgc->gcsl_start4(p,b->press,gcval_press);
  
	ucorr(p,b,P,eta,alpha);
	vcorr(p,b,Q,eta,alpha);
    wcorr(p,b, alpha, P, Q, ws);
	pgc->gcsl_start4(p,ws,12);
	pgc->gcsl_start4(p,wb,12);

    p->poissoniter=p->solveriter;

	ptime=pgc->timer()-starttime;

	if(p->mpirank==0 && (p->count%p->P12==0) && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  solvtime: "<<setprecision(3)<<p->poissontime<<"  ptime: "<<setprecision(3)<<ptime<<endl;
    
}

void sflow_pjm_sw::ucorr(lexer* p, fdm2D* b, slice& P, slice &eta, double alpha)
{	
	SLICELOOP1
    if(b->breaking(i,j)==0 && b->breaking(i+1,j)==0)
	P(i,j) -= alpha*p->dt*(((b->press(i+1,j)-b->press(i,j))/(2.0*p->DXM)));

                
    SLICELOOP1
    if(b->breaking(i,j)==0 && b->breaking(i+1,j)==0)
	P(i,j) -= alpha*p->dt*(0.5*(b->press(i+1,j)+b->press(i,j))*((eta(i+1,j)-eta(i,j))-(b->depth(i+1,j)-b->depth(i,j)))
				             /(2.0*p->DXM*HXIJ ));
}

void sflow_pjm_sw::vcorr(lexer* p, fdm2D* b, slice& Q, slice &eta, double alpha)
{	
	SLICELOOP2
    if(b->breaking(i,j)==0 && b->breaking(i,j+1)==0)
	Q(i,j) -= alpha*p->dt*(((b->press(i,j+1)-b->press(i,j))/(2.0*p->DXM)));
                
    SLICELOOP2
    if(b->breaking(i,j)==0 && b->breaking(i,j+1)==0)
	Q(i,j) -= alpha*p->dt*(0.5*(b->press(i,j+1)+b->press(i,j))*((eta(i,j+1)-eta(i,j))-(b->depth(i,j+1)-b->depth(i,j)))
                            /(2.0*p->DXM*HYIJ ));
}

void sflow_pjm_sw::wcorr(lexer* p, fdm2D* b,double alpha, slice &P, slice &Q, slice &ws)
{	    
    SLICELOOP4
    if(b->breaking(i,j)==0)
	ws(i,j) += p->dt*alpha*(2.0*b->press(i,j)/HP);
}

void sflow_pjm_sw::wcalc(lexer* p, fdm2D* b,double alpha, slice &P, slice &Q, slice &ws)
{	
    SLICELOOP4
    wbn(i,j) = wb(i,j);

	
	SLICELOOP4
    wb(i,j) = -0.25*(P(i,j)+P(i-1,j))*(b->depth(i+1,j)-b->depth(i-1,j))/p->DXM
                
                -0.25*(Q(i,j)+Q(i,j-1))*(b->depth(i,j+1)-b->depth(i,j-1))/p->DXM;
	
    SLICELOOP4
	ws(i,j) += -(wb(i,j)-wbn(i,j));
}

void sflow_pjm_sw::rhs(lexer *p, fdm2D* b, slice &u, slice &v, slice &ws, double alpha)
{
    NSLICELOOP4
	b->rhsvec.V[n]=0.0;
    
	count=0;
    SLICELOOP4
    {
    b->rhsvec.V[count] =   -((u(i,j) - u(i-1,j))*(b->hp(i,j))
                           + (v(i,j) - v(i,j-1))*(b->hp(i,j)))/(alpha*p->dt*p->DXM)
                           
                           -(ws(i,j)-wb(i,j))/(alpha*p->dt);
    ++count;
    }
}

void sflow_pjm_sw::upgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{
    SLICELOOP1
    WETDRY1
	b->F(i,j) -= fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j) )/(p->DXM); 
    
    pBC->patchBC_pressure2D_ugrad(p,b,eta,eta_n);
}

void sflow_pjm_sw::vpgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{
    SLICELOOP2
    WETDRY2
	b->G(i,j) -= fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j) )/(p->DXM); 
    
    pBC->patchBC_pressure2D_vgrad(p,b,eta,eta_n);
}

void sflow_pjm_sw::poisson(lexer*p, fdm2D* b)
{
    sqd = (1.0/(2.0*p->DXM*p->DXM));

    n=0;
    SLICELOOP4
	{
	b->M.p[n]  =  (b->hp(i,j)*sqd + b->hp(i,j)*sqd)*p->x_dir
               +  (b->hp(i,j)*sqd + b->hp(i,j)*sqd)*p->y_dir + 2.0/HP;

   	b->M.n[n] = -b->hp(i,j)*sqd*p->x_dir;
	b->M.s[n] = -b->hp(i,j)*sqd*p->x_dir;

	b->M.w[n] = -b->hp(i,j)*sqd*p->y_dir;
	b->M.e[n] = -b->hp(i,j)*sqd*p->y_dir;

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

void sflow_pjm_sw::wpgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{	    
}

