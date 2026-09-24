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

#include"sflow_pjm_lin.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"solver2D.h"
#include"ioflow.h"
#include"patchBC_interface.h"

#define HP (WL(i,j)>1.0e-20?WL(i,j):1.0e20)

sflow_pjm_lin::sflow_pjm_lin(lexer* p, fdm2D *b, patchBC_interface *ppBC) : cb(2.0)
{
    pBC = ppBC;
    
    gcval_press=40;  
}

sflow_pjm_lin::~sflow_pjm_lin()
{
}

void sflow_pjm_lin::start(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, ioflow *pflow, 
                          slice &UH, slice &VH, slice &WH, slice &WL, slice &Un, slice &Vn, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
    
	starttime=pgc->timer();

    rhs(p,b,WL,alpha);
    poisson(p,b,WL,alpha);

        solvtime=pgc->timer();

    psolv->start(p,pgc,b->press,b->M,b->xvec,b->rhsvec,4);

        p->poissontime=pgc->timer()-solvtime;
    
    pflow->pm_relax(p,pgc,b->press);
	pgc->gcsl_start4(p,b->press,gcval_press);

	ucorr(p,b,UH,WL,alpha);
	vcorr(p,b,VH,WL,alpha);
    wcorr(p,b,WH,WL,alpha);

    p->poissoniter=p->solveriter;

	ptime=pgc->timer()-starttime;

	if(p->mpirank==0 && (p->count%p->P12==0))
	cout<<"piter: "<<p->solveriter<<"  solvtime: "<<setprecision(3)<<p->poissontime<<"  ptime: "<<setprecision(3)<<ptime<<endl;
}

int sflow_pjm_lin::active(lexer *p, fdm2D *b)
{
    // non-hydrostatic pressure is solved in wet, deep, non-breaking cells
    if(p->wet[IJ]==0 || p->deep[IJ]==0 || b->breaking(i,j)>0)
    return 0;
    
    if(p->wet[Im1J]==0 || p->wet[Ip1J]==0)
    return 0;
    
    if(p->j_dir==1 && (p->wet[IJm1]==0 || p->wet[IJp1]==0))
    return 0;
    
    return 1;
}

void sflow_pjm_lin::ucorr(lexer* p, fdm2D* b, slice &UH, slice &WL, double alpha)
{	
    double dx;
    
	SLICELOOP4
    if(active(p,b)==1)
    {
    dx = p->DXP[IP] + p->DXP[IM1];
    
	UH(i,j) -= alpha*p->dt*(1.0/p->W1)*((WL(i+1,j)*b->press(i+1,j) - WL(i-1,j)*b->press(i-1,j))/dx
                                    - cb*b->press(i,j)*(b->depth(i+1,j) - b->depth(i-1,j))/dx);
    }
}

void sflow_pjm_lin::vcorr(lexer* p, fdm2D* b, slice &VH, slice &WL, double alpha)
{	
    double dy;
    
    if(p->j_dir==1)
	SLICELOOP4
    if(active(p,b)==1)
    {
    dy = p->DYP[JP] + p->DYP[JM1];
    
	VH(i,j) -= alpha*p->dt*(1.0/p->W1)*((WL(i,j+1)*b->press(i,j+1) - WL(i,j-1)*b->press(i,j-1))/dy
                                    - cb*b->press(i,j)*(b->depth(i,j+1) - b->depth(i,j-1))/dy);
    }
}

void sflow_pjm_lin::wcorr(lexer* p, fdm2D* b, slice &WH, slice &WL, double alpha)
{	    
    SLICELOOP4
    if(active(p,b)==1)
	WH(i,j) += alpha*p->dt*(1.0/p->W1)*cb*b->press(i,j);
}

void sflow_pjm_lin::rhs(lexer *p, fdm2D* b, slice &WL, double alpha)
{
    double dx,dy,dudx,dvdy,dddx,dddy;
    
    b->rhsvec.reset();
    
    n=0;
    SLICELOOP4
    {
        if(active(p,b)==1)
        {
        dx = p->DXP[IP] + p->DXP[IM1];
        dy = p->DYP[JP] + p->DYP[JM1];
        
        dudx = (b->U(i+1,j) - b->U(i-1,j))/dx;
        dvdy = (b->V(i,j+1) - b->V(i,j-1))/dy*p->y_dir;
        
        dddx = (b->depth(i+1,j) - b->depth(i-1,j))/dx;
        dddy = (b->depth(i,j+1) - b->depth(i,j-1))/dy*p->y_dir;
        
        b->rhsvec.V[n] = -(WL(i,j)*(dudx + dvdy) 
                          + 2.0*(b->W(i,j) + b->U(i,j)*dddx + b->V(i,j)*dddy))/(alpha*p->dt);
        }
    ++n;
    }
}

void sflow_pjm_lin::poisson(lexer*p, fdm2D* b, slice &WL, double alpha)
{
    n=0;
    SLICELOOP4
	{
	b->M.p[n] =  WL(i,j)/(p->W1*p->DXP[IP]*p->DXN[IP])
               + WL(i,j)/(p->W1*p->DXP[IM1]*p->DXN[IP])
               + WL(i,j)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir
               + WL(i,j)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir
               + 2.0*cb/(HP*p->W1);
    
   	b->M.n[n] = -WL(i,j)/(p->W1*p->DXP[IP]*p->DXN[IP]);
	b->M.s[n] = -WL(i,j)/(p->W1*p->DXP[IM1]*p->DXN[IP]);
	b->M.w[n] = -WL(i,j)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir;
	b->M.e[n] = -WL(i,j)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir;
	++n;
	}
    
    n=0;
	SLICELOOP4
	{
        // inflow / outflow: q=0
		if(p->flagslice4[Im1J]<0 && p->IOSL[Im1J]==1)
		b->M.s[n] = 0.0;
        
        if(p->flagslice4[Ip1J]<0 && p->IOSL[Ip1J]==1)
		b->M.n[n] = 0.0;
        
        // walls: dq/dn=0
        if(p->flagslice4[Im1J]<0 && p->IOSL[Im1J]==0)
		{
		b->M.p[n] += b->M.s[n];
		b->M.s[n] = 0.0;
		}

        if(p->flagslice4[Ip1J]<0 && p->IOSL[Ip1J]==0)
		{
		b->M.p[n] += b->M.n[n];
		b->M.n[n] = 0.0;
		}

		if(p->flagslice4[IJm1]<0)
		{
		b->M.p[n] += b->M.e[n];
		b->M.e[n] = 0.0;
		}

		if(p->flagslice4[IJp1]<0)
		{
		b->M.p[n] += b->M.w[n];
		b->M.w[n] = 0.0;
		}
        
        // hydrostatic cells
        if(active(p,b)==0)
        {
        b->M.p[n] = 1.0;
        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;
        b->M.w[n] = 0.0;
        b->M.e[n] = 0.0;
        b->rhsvec.V[n] = 0.0;
        b->press(i,j) = 0.0;
        }
	++n;
	}
}

void sflow_pjm_lin::upgrad(lexer*p, fdm2D* b, slice &eta)
{
    SLICELOOP4
    WETDRY
    b->F(i,j) += fabs(p->W22)*eta(i,j)*(b->dfx(i,j) - b->dfx(i-1,j))/p->DXN[IP];
}

void sflow_pjm_lin::vpgrad(lexer*p, fdm2D* b, slice &eta)
{
    SLICELOOP4
    WETDRY
    b->G(i,j) += fabs(p->W22)*eta(i,j)*(b->dfy(i,j) - b->dfy(i,j-1))/p->DYN[JP];
}

void sflow_pjm_lin::wpgrad(lexer*p, fdm2D* b, slice &eta)
{
}
