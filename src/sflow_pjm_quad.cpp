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

#include"sflow_pjm_quad.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"solver2D.h"
#include"ioflow.h"
#include"patchBC_interface.h"

#define HP (WL(i,j)>1.0e-20?WL(i,j):1.0e20)

sflow_pjm_quad::sflow_pjm_quad(lexer* p, fdm2D *b, ghostcell *ppgc, patchBC_interface *ppBC) : cb(1.5), phi(p), Uest(p), Vest(p), Ld(p), Gx(p), Gy(p)
{
    pgc = ppgc;
    
    // dispersion parameter
    adisp = 1.0;
    
    if(p->A220==3)
    adisp = MAX(p->A224,1.0);
    
    cw = cb/adisp;
    Bdisp = (adisp-1.0)/3.0;
    
    if(p->mpirank==0 && p->A220==3)
    cout<<"SFLOW quadratic non-hydrostatic pressure with improved dispersion, alpha: "<<adisp<<endl;

    pBC = ppBC;
    
    gcval_press=40;  
}

sflow_pjm_quad::~sflow_pjm_quad()
{
}

void sflow_pjm_quad::start(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, ioflow *pflow, 
                           slice &UH, slice &VH, slice &WH, slice &WL, slice &Un, slice &Vn, double alpha)
{
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<".";
    
	starttime=pgc->timer();
    
    quad_calc(p,b,pgc,UH,VH,WL,Un,Vn,alpha);
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

int sflow_pjm_quad::active(lexer *p, fdm2D *b)
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

void sflow_pjm_quad::quad_calc(lexer *p, fdm2D *b, ghostcell *pgc, slice &UH, slice &VH, slice &WL, slice &Un, slice &Vn, double alpha)
{
    double dx,dy;
    double dddx,dddy,dddxx,dddyy,dddxy;
    double dudt,dvdt;
    
    // estimate of the corrected velocities with the previous pressure
    SLICELOOP4
    {
    Uest(i,j) = b->U(i,j);
    Vest(i,j) = b->V(i,j);
    }
    
    SLICELOOP4
    if(active(p,b)==1)
    {
    dx = p->DXP[IP] + p->DXP[IM1];
    dy = p->DYP[JP] + p->DYP[JM1];
    
    Uest(i,j) = (UH(i,j) - alpha*p->dt*(1.0/p->W1)*((WL(i+1,j)*b->press(i+1,j) - WL(i-1,j)*b->press(i-1,j))/dx
                                               - cb*b->press(i,j)*(b->depth(i+1,j) - b->depth(i-1,j))/dx))/HP;
    
    if(p->j_dir==1)
    Vest(i,j) = (VH(i,j) - alpha*p->dt*(1.0/p->W1)*((WL(i,j+1)*b->press(i,j+1) - WL(i,j-1)*b->press(i,j-1))/dy
                                               - cb*b->press(i,j)*(b->depth(i,j+1) - b->depth(i,j-1))/dy))/HP;
    }
    
    // ghost cells: boundary velocities
    GCSL4LOOP
    {
    i = p->gcbsl4[n][0];
    j = p->gcbsl4[n][1];
    
        for(q=1;q<=3;++q)
        {
        if(p->gcbsl4[n][3]==1)
        {
        Uest(i-q,j) = b->U(i-q,j);
        Vest(i-q,j) = b->V(i-q,j);
        }
        
        if(p->gcbsl4[n][3]==4)
        {
        Uest(i+q,j) = b->U(i+q,j);
        Vest(i+q,j) = b->V(i+q,j);
        }
        
        if(p->gcbsl4[n][3]==2)
        {
        Uest(i,j+q) = b->U(i,j+q);
        Vest(i,j+q) = b->V(i,j+q);
        }
        
        if(p->gcbsl4[n][3]==3)
        {
        Uest(i,j-q) = b->U(i,j-q);
        Vest(i,j-q) = b->V(i,j-q);
        }
        }
    }
    
    pgc->gcsl_start4(p,Uest,10);
    pgc->gcsl_start4(p,Vest,10);
    
    // bed acceleration phi = Dw_b/Dt,  w_b = -u.grad(d)
    SLICELOOP4
    {
    phi(i,j) = 0.0;
    
        if(active(p,b)==1)
        {
        dx = p->DXP[IP] + p->DXP[IM1];
        dy = p->DYP[JP] + p->DYP[JM1];
        
        dddx  = (b->depth(i+1,j) - b->depth(i-1,j))/dx;
        dddxx = ((b->depth(i+1,j) - b->depth(i,j))/p->DXP[IP] - (b->depth(i,j) - b->depth(i-1,j))/p->DXP[IM1])/p->DXN[IP];
        
        dudt = (Uest(i,j) - Un(i,j))/(alpha*p->dt)
             + Uest(i,j)*(Uest(i+1,j) - Uest(i-1,j))/dx
             + Vest(i,j)*(Uest(i,j+1) - Uest(i,j-1))/dy*p->y_dir;
        
        phi(i,j) = -dudt*dddx - Uest(i,j)*Uest(i,j)*dddxx;
        
            if(p->j_dir==1)
            {
            dddy  = (b->depth(i,j+1) - b->depth(i,j-1))/dy;
            dddyy = ((b->depth(i,j+1) - b->depth(i,j))/p->DYP[JP] - (b->depth(i,j) - b->depth(i,j-1))/p->DYP[JM1])/p->DYN[JP];
            
            dddxy = 0.0;
            
            if(p->flagslice4[Ip1Jp1]>0 && p->flagslice4[Ip1Jm1]>0 && p->flagslice4[Im1Jp1]>0 && p->flagslice4[Im1Jm1]>0)
            dddxy = (b->depth(i+1,j+1) - b->depth(i+1,j-1) - b->depth(i-1,j+1) + b->depth(i-1,j-1))/(dx*dy);
            
            dvdt = (Vest(i,j) - Vn(i,j))/(alpha*p->dt)
                 + Uest(i,j)*(Vest(i+1,j) - Vest(i-1,j))/dx
                 + Vest(i,j)*(Vest(i,j+1) - Vest(i,j-1))/dy;
            
            phi(i,j) += -dvdt*dddy - Vest(i,j)*Vest(i,j)*dddyy - 2.0*Uest(i,j)*Vest(i,j)*dddxy;
            }
        }
    }
}

void sflow_pjm_quad::ucorr(lexer* p, fdm2D* b, slice &UH, slice &WL, double alpha)
{	
    double dx,qb;
    
	SLICELOOP4
    if(active(p,b)==1)
    {
    dx = p->DXP[IP] + p->DXP[IM1];
    
    qb = cb*b->press(i,j) + adisp*0.25*p->W1*WL(i,j)*phi(i,j);
    
	UH(i,j) -= alpha*p->dt*(1.0/p->W1)*((WL(i+1,j)*b->press(i+1,j) - WL(i-1,j)*b->press(i-1,j))/dx
                                    - qb*(b->depth(i+1,j) - b->depth(i-1,j))/dx);
    }
}

void sflow_pjm_quad::vcorr(lexer* p, fdm2D* b, slice &VH, slice &WL, double alpha)
{	
    double dy,qb;
    
    if(p->j_dir==1)
	SLICELOOP4
    if(active(p,b)==1)
    {
    dy = p->DYP[JP] + p->DYP[JM1];
    
    qb = cb*b->press(i,j) + adisp*0.25*p->W1*WL(i,j)*phi(i,j);
    
	VH(i,j) -= alpha*p->dt*(1.0/p->W1)*((WL(i,j+1)*b->press(i,j+1) - WL(i,j-1)*b->press(i,j-1))/dy
                                    - qb*(b->depth(i,j+1) - b->depth(i,j-1))/dy);
    }
}

void sflow_pjm_quad::wcorr(lexer* p, fdm2D* b, slice &WH, slice &WL, double alpha)
{	    
    SLICELOOP4
    if(active(p,b)==1)
	WH(i,j) += alpha*p->dt*((1.0/p->W1)*cw*b->press(i,j) + 0.25*WL(i,j)*phi(i,j));
}

void sflow_pjm_quad::rhs(lexer *p, fdm2D* b, slice &WL, double alpha)
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
                          + 2.0*(b->W(i,j) + b->U(i,j)*dddx + b->V(i,j)*dddy))/(alpha*p->dt)
                         - 0.5*phi(i,j);
        }
    ++n;
    }
}

void sflow_pjm_quad::poisson(lexer*p, fdm2D* b, slice &WL, double alpha)
{
    n=0;
    SLICELOOP4
	{
	b->M.p[n] =  WL(i,j)/(p->W1*p->DXP[IP]*p->DXN[IP])
               + WL(i,j)/(p->W1*p->DXP[IM1]*p->DXN[IP])
               + WL(i,j)/(p->W1*p->DYP[JP]*p->DYN[JP])*p->y_dir
               + WL(i,j)/(p->W1*p->DYP[JM1]*p->DYN[JP])*p->y_dir
               + 2.0*cw/(HP*p->W1);
    
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

void sflow_pjm_quad::upgrad(lexer*p, fdm2D* b, slice &eta)
{
    SLICELOOP4
    WETDRY
    b->F(i,j) += fabs(p->W22)*eta(i,j)*(b->dfx(i,j) - b->dfx(i-1,j))/p->DXN[IP];
    
    // dispersion correction (A 220 3): B g grad(h^3 div(grad(eta)))
    // all derivatives central: the collocated projection only acts through central
    // gradients, a compact Laplacian here would limit the time step to dt ~ dx^2 much earlier
    if(Bdisp>0.0)
    {
    double h;
    
        SLICELOOP4
        {
        Gx(i,j) = (eta(i+1,j) - eta(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
        Gy(i,j) = (eta(i,j+1) - eta(i,j-1))/(p->DYP[JP] + p->DYP[JM1])*p->y_dir;
        }
        
        pgc->gcsl_start4(p,Gx,1);
        pgc->gcsl_start4(p,Gy,1);
        
        SLICELOOP4
        {
        Ld(i,j) = 0.0;
        
            if(active(p,b)==1)
            {
            h = MAX(eta(i,j) + b->depth(i,j),0.0);
            
            Ld(i,j) = h*h*h*((Gx(i+1,j) - Gx(i-1,j))/(p->DXP[IP] + p->DXP[IM1])
                           + (Gy(i,j+1) - Gy(i,j-1))/(p->DYP[JP] + p->DYP[JM1])*p->y_dir);
            }
        }
        
        pgc->gcsl_start4(p,Ld,1);
        
        SLICELOOP4
        if(active(p,b)==1)
        b->F(i,j) += Bdisp*fabs(p->W22)*(Ld(i+1,j) - Ld(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
    }
}

void sflow_pjm_quad::vpgrad(lexer*p, fdm2D* b, slice &eta)
{
    SLICELOOP4
    WETDRY
    b->G(i,j) += fabs(p->W22)*eta(i,j)*(b->dfy(i,j) - b->dfy(i,j-1))/p->DYN[JP];
    
    // dispersion correction (A 220 3), Ld from upgrad
    if(Bdisp>0.0 && p->j_dir==1)
    SLICELOOP4
    if(active(p,b)==1)
    b->G(i,j) += Bdisp*fabs(p->W22)*(Ld(i,j+1) - Ld(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
}

void sflow_pjm_quad::wpgrad(lexer*p, fdm2D* b, slice &eta)
{
}
