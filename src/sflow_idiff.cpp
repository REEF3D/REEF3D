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

#include"sflow_idiff.h"
#include"sflow_diffusion_flux.h"
#include"lexer.h"
#include"ghostcell.h"
#include"fdm2D.h"
#include"solver2D.h"

sflow_idiff::sflow_idiff(lexer* p) : f(p)
{
    gcval_u = 10;
    gcval_v = 11;
}

sflow_idiff::~sflow_idiff()
{
}

void sflow_idiff::bc(lexer *p, fdm2D *b, slice &u)
{
    // boundary values from the ghost cells of the velocity
    n=0;
    SLICELOOP4
    {
        if(p->flagslice4[Im1J]<0)
		{
		b->rhsvec.V[n] -= b->M.s[n]*u(i-1,j);
		b->M.s[n] = 0.0;
		}
		
		if(p->flagslice4[Ip1J]<0)
		{
		b->rhsvec.V[n] -= b->M.n[n]*u(i+1,j);
		b->M.n[n] = 0.0;
		}
		
		if(p->flagslice4[IJm1]<0)
		{
		b->rhsvec.V[n] -= b->M.e[n]*u(i,j-1);
		b->M.e[n] = 0.0;
		}
		
		if(p->flagslice4[IJp1]<0)
		{
		b->rhsvec.V[n] -= b->M.w[n]*u(i,j+1);
		b->M.w[n] = 0.0;
		}
        
        // dry neighbours inside the domain: free slip, zero normal gradient (the coupling goes to the diagonal);
        // the dry cell itself is solved as f = 0 and would otherwise act as a no-slip wall at the shoreline
        if(p->wet[IJ]==1)
        {
            if(p->flagslice4[Im1J]>0 && p->wet[Im1J]==0)
            {
            b->M.p[n] += b->M.s[n];
            b->M.s[n] = 0.0;
            }
            
            if(p->flagslice4[Ip1J]>0 && p->wet[Ip1J]==0)
            {
            b->M.p[n] += b->M.n[n];
            b->M.n[n] = 0.0;
            }
            
            if(p->flagslice4[IJm1]>0 && p->wet[IJm1]==0)
            {
            b->M.p[n] += b->M.e[n];
            b->M.e[n] = 0.0;
            }
            
            if(p->flagslice4[IJp1]>0 && p->wet[IJp1]==0)
            {
            b->M.p[n] += b->M.w[n];
            b->M.w[n] = 0.0;
            }
        }
        
        if(p->wet[IJ]==0)
        {
        b->M.p[n] = 1.0;
        b->M.n[n] = 0.0;
        b->M.s[n] = 0.0;
        b->M.e[n] = 0.0;
        b->M.w[n] = 0.0;
        b->rhsvec.V[n] = 0.0;
        }
	++n;
	}
}

void sflow_idiff::diff_u(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &UHdiff, slice &UH, slice &U, slice &V, slice &WL, double alpha)
{
    double cxp,cxm,cyp,cym,dxy;
    
    starttime=pgc->timer();
    
    // flux form (sflow_diffusion_flux.h), divided by the depth of the cell: h_P du/dt = sum_f (h nu)_f du/dn_f
    n=0;
    SLICELOOP4
    {
    cxp=cxm=cyp=cym=dxy=0.0;

        if(p->wet[IJ]==1 && WL(i,j)>0.0)
        {
        const double ih = 1.0/WL(i,j);

        cxp = 2.0*sflow_dflux::hnu(p,b,WL,i,j,i+1,j)*ih/(p->DXP[IP]*p->DXN[IP]);
        cxm = 2.0*sflow_dflux::hnu(p,b,WL,i,j,i-1,j)*ih/(p->DXP[IM1]*p->DXN[IP]);
        cyp = sflow_dflux::hnu(p,b,WL,i,j,i,j+1)*ih/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
        cym = sflow_dflux::hnu(p,b,WL,i,j,i,j-1)*ih/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;

        dxy = sflow_dflux::cross(p,b,WL,V,i,j,1)*ih;
        }
    
	b->M.p[n] = cxp + cxm + cyp + cym + 1.0/(alpha*p->dt);
    b->M.n[n] = -cxp;
    b->M.s[n] = -cxm;
    b->M.w[n] = -cyp;
    b->M.e[n] = -cym;
    
	b->rhsvec.V[n] = dxy + U(i,j)/(alpha*p->dt);
    
    f(i,j) = U(i,j);
	++n;
	}
    
    bc(p,b,U);
    
	psolv->start(p,pgc,f,b->M,b->xvec,b->rhsvec,4);
    
    SLICELOOP4
    UHdiff(i,j) = (p->wet[IJ]==1)?UH(i,j) + WL(i,j)*(f(i,j) - U(i,j)):UH(i,j);
    
	time=pgc->timer()-starttime;
	p->uiter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<"udiffiter: "<<p->uiter<<"  udifftime: "<<setprecision(3)<<time<<endl;
}

void sflow_idiff::diff_v(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &VHdiff, slice &VH, slice &U, slice &V, slice &WL, double alpha)
{
    double cxp,cxm,cyp,cym,dxy;
    
    if(p->j_dir==0)
    {
    SLICELOOP4
    VHdiff(i,j) = VH(i,j);
    
    return;
    }
    
    starttime=pgc->timer();
    
    n=0;
    SLICELOOP4
    {
    cxp=cxm=cyp=cym=dxy=0.0;

        if(p->wet[IJ]==1 && WL(i,j)>0.0)
        {
        const double ih = 1.0/WL(i,j);

        cxp = sflow_dflux::hnu(p,b,WL,i,j,i+1,j)*ih/(p->DXP[IP]*p->DXN[IP]);
        cxm = sflow_dflux::hnu(p,b,WL,i,j,i-1,j)*ih/(p->DXP[IM1]*p->DXN[IP]);
        cyp = 2.0*sflow_dflux::hnu(p,b,WL,i,j,i,j+1)*ih/(p->DYP[JP]*p->DYN[JP]);
        cym = 2.0*sflow_dflux::hnu(p,b,WL,i,j,i,j-1)*ih/(p->DYP[JM1]*p->DYN[JP]);

        dxy = sflow_dflux::cross(p,b,WL,U,i,j,2)*ih;
        }
    
	b->M.p[n] = cxp + cxm + cyp + cym + 1.0/(alpha*p->dt);
    b->M.n[n] = -cxp;
    b->M.s[n] = -cxm;
    b->M.w[n] = -cyp;
    b->M.e[n] = -cym;
    
	b->rhsvec.V[n] = dxy + V(i,j)/(alpha*p->dt);
    
    f(i,j) = V(i,j);
	++n;
	}
    
    bc(p,b,V);
    
	psolv->start(p,pgc,f,b->M,b->xvec,b->rhsvec,4);
    
    SLICELOOP4
    VHdiff(i,j) = (p->wet[IJ]==1)?VH(i,j) + WL(i,j)*(f(i,j) - V(i,j)):VH(i,j);
    
	time=pgc->timer()-starttime;
	p->viter=p->solveriter;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<"vdiffiter: "<<p->viter<<"  vdifftime: "<<setprecision(3)<<time<<endl;
}

void sflow_idiff::diff_w(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &WHdiff, slice &WH, slice &U, slice &V, slice &W, slice &WL, double alpha)
{
    double cxp,cxm,cyp,cym;
    
    starttime=pgc->timer();
    
    n=0;
    SLICELOOP4
    {
    cxp=cxm=cyp=cym=0.0;

        if(p->wet[IJ]==1 && WL(i,j)>0.0)
        {
        const double ih = 1.0/WL(i,j);

        cxp = sflow_dflux::hnu(p,b,WL,i,j,i+1,j)*ih/(p->DXP[IP]*p->DXN[IP]);
        cxm = sflow_dflux::hnu(p,b,WL,i,j,i-1,j)*ih/(p->DXP[IM1]*p->DXN[IP]);
        cyp = sflow_dflux::hnu(p,b,WL,i,j,i,j+1)*ih/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
        cym = sflow_dflux::hnu(p,b,WL,i,j,i,j-1)*ih/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;
        }
    
	b->M.p[n] = cxp + cxm + cyp + cym + 1.0/(alpha*p->dt);
    b->M.n[n] = -cxp;
    b->M.s[n] = -cxm;
    b->M.w[n] = -cyp;
    b->M.e[n] = -cym;
    
	b->rhsvec.V[n] = W(i,j)/(alpha*p->dt);
    
    f(i,j) = W(i,j);
	++n;
	}
    
    bc(p,b,W);
    
	psolv->start(p,pgc,f,b->M,b->xvec,b->rhsvec,4);
    
    SLICELOOP4
    WHdiff(i,j) = (p->wet[IJ]==1)?WL(i,j)*f(i,j):WH(i,j);
    
	time=pgc->timer()-starttime;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<"wdiffiter: "<<p->solveriter<<"  wdifftime: "<<setprecision(3)<<time<<endl;
}

void sflow_idiff::diff_scalar(lexer* p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &f, double sig, double alpha)
{
    count=0;

    sqd = (1.0/(p->DXM*p->DXM));


	SLICELOOP4
	{
	ev_ij=b->eddyv(i,j);
	visc_ij=p->W2;

	b->M.p[count]  +=   0.5*sqd*(visc_ij+b->eddyv(i+1,j)/sig + visc_ij+ev_ij/sig)
    
					+   0.5*sqd*(visc_ij+ev_ij/sig + visc_ij+b->eddyv(i-1,j)/sig)
                    
					+   0.5*sqd*(visc_ij+b->eddyv(i,j+1)/sig + visc_ij+ev_ij/sig)*p->y_dir
                    
					+   0.5*sqd*(visc_ij+ev_ij/sig + visc_ij+b->eddyv(i,j-1)/sig)*p->y_dir;

	 
	 b->M.s[count] -= 0.5*sqd*(visc_ij+ev_ij/sig + visc_ij+b->eddyv(i-1,j)/sig);
	 b->M.n[count] -= 0.5*sqd*(visc_ij+b->eddyv(i+1,j)/sig + visc_ij+ev_ij/sig);
	 
	 b->M.e[count] -= 0.5*sqd*(visc_ij+ev_ij/sig + visc_ij+b->eddyv(i,j-1)/sig)*p->y_dir;
	 b->M.w[count] -= 0.5*sqd*(visc_ij+b->eddyv(i,j+1)/sig + visc_ij+ev_ij/sig)*p->y_dir;

	 ++count;
	}
}
