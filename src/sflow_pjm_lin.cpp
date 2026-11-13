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
#include"sflow_amr.h"
#include"sflow_curv.h"
#include"gcsl_yskip.h"

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
    poisson(p,b,WL,cb);

        solvtime=pgc->timer();

    // mesh refinement: one pressure for level 0 and the patches
    if(amr!=nullptr && amr->nh_patches())
    amr->nh_solve(p,b,pgc,UH,VH,WH,WL,alpha);
    
    else
    // G 7 1: level 0 alone, the cells under the patches fixed at their restricted pressure
    if(amr!=nullptr && amr->nh_level0())
    amr->nh_solve0(p,b,pgc);
    
    else
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

void sflow_pjm_lin::assemble(lexer *p, fdm2D *b, ghostcell *pgc, slice &UH, slice &VH, slice &WL, slice &Un, slice &Vn, double alpha)
{
    rhs(p,b,WL,alpha);
    poisson(p,b,WL,cb);
}

void sflow_pjm_lin::correct(lexer *p, fdm2D *b, slice &UH, slice &VH, slice &WH, slice &WL, double alpha)
{
    ucorr(p,b,UH,WL,alpha);
	vcorr(p,b,VH,WL,alpha);
    wcorr(p,b,WH,WL,alpha);
}

void sflow_pjm_lin::ucorr(lexer* p, fdm2D* b, slice &UH, slice &WL, double alpha)
{	
    double dx;
    
    // curvilinear grid: central gradients by the chain rule
    if(p->geo_curv==1)
    {
    double gx,gy,dx_,dy_;

	SLICELOOP4
    if(active(p,b)==1)
    {
    grad_q(p,b,WL,gx,gy);
    sflow_curv::grad(p,b->depth,i,j,dx_,dy_);
    
	UH(i,j) -= alpha*p->dt*(1.0/p->W1)*(gx - cb*b->press(i,j)*dx_);
    }
    return;
    }

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
    
    // curvilinear grid: central gradients by the chain rule
    if(p->geo_curv==1)
    {
    double gx,gy,dx_,dy_;

	SLICELOOP4
    if(active(p,b)==1)
    {
    grad_q(p,b,WL,gx,gy);
    sflow_curv::grad(p,b->depth,i,j,dx_,dy_);
    
	VH(i,j) -= alpha*p->dt*(1.0/p->W1)*(gy - cb*b->press(i,j)*dy_);
    }
    return;
    }

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

    // curvilinear grid: chain-rule gradients, the non-orthogonal part of the operator on the
    // previous pressure (deferred correction)
    if(p->geo_curv==1)
    {
    double ux,uy,vx,vy;

    n=0;
    SLICELOOP4
    {
        if(active(p,b)==1)
        {
        sflow_curv::grad(p,b->U,i,j,ux,uy);
        sflow_curv::grad(p,b->V,i,j,vx,vy);
        sflow_curv::grad(p,b->depth,i,j,dddx,dddy);

        b->rhsvec.V[n] = -(WL(i,j)*(ux + vy)
                          + 2.0*(b->W(i,j) + b->U(i,j)*dddx + b->V(i,j)*dddy))/(alpha*p->dt);

        const int si=i, sj=j;
        const double kf = sflow_curv::kflux(p,b->press,si,sj,[&](int a, int c)
                          {
                              if(p->flagslice4[(a-p->imin)*p->jmax+(c-p->jmin)]<0) return false;
                              i=a; j=c; const int r=active(p,b); i=si; j=sj;
                              return r==1;
                          });
        b->rhsvec.V[n] += WL(i,j)/p->W1*kf;
        }
    ++n;
    }
    return;
    }

    const bool ysk = gcsl_yskip(p);
    n=0;
    SLICELOOP4
    {
        if(active(p,b)==1)
        {
        dx = p->DXP[IP] + p->DXP[IM1];
        dy = p->DYP[JP] + p->DYP[JM1];
        
        dudx = (b->U(i+1,j) - b->U(i-1,j))/dx;
        // (one cell in y: the y derivatives are 0, as from the ghost cells; gcsl_yskip leaves those unset)
        dvdy = ysk ? 0.0 : (b->V(i,j+1) - b->V(i,j-1))/dy*p->y_dir;
        
        dddx = (b->depth(i+1,j) - b->depth(i-1,j))/dx;
        dddy = ysk ? 0.0 : (b->depth(i,j+1) - b->depth(i,j-1))/dy*p->y_dir;
        
        b->rhsvec.V[n] = -(WL(i,j)*(dudx + dvdy) 
                          + 2.0*(b->W(i,j) + b->U(i,j)*dddx + b->V(i,j)*dddy))/(alpha*p->dt);
        }
    ++n;
    }
}

void sflow_pjm_lin::upgrad(lexer*p, fdm2D* b, slice &eta)
{
    bedslope_x(p,b,eta);
}

void sflow_pjm_lin::vpgrad(lexer*p, fdm2D* b, slice &eta)
{
    bedslope_y(p,b,eta);
}

void sflow_pjm_lin::wpgrad(lexer*p, fdm2D* b, slice &eta)
{
}

// curvilinear grid: gradient of h q at cell (i,j)
void sflow_pjm_lin::grad_q(lexer *p, fdm2D *b, slice &WL, double &gx, double &gy)
{
    sflow_curv::grad5(p,i,j,WL(i+1,j)*b->press(i+1,j),WL(i-1,j)*b->press(i-1,j),
                            WL(i,j+1)*b->press(i,j+1),WL(i,j-1)*b->press(i,j-1),gx,gy);
}
