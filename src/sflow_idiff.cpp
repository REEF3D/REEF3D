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
#include"sflow_curv.h"

sflow_idiff::sflow_idiff(lexer* p) : f(p)
{
    gcval_u = 10;
    gcval_v = 11;
}

sflow_idiff::~sflow_idiff()
{
    delete own_solver;
}

// the linear solve of a component: level 0 with the solver of the model; a patch of the mesh refinement with its own
// solver, the ghost cells of f with the values of the field (the cells around the patch are Dirichlet values of the
// solver), the reductions local to the patch
void sflow_idiff::solve(lexer *p, ghostcell *pgc, solver2D *psolv, fdm2D *b, slice &u)
{
    if(own_solver==nullptr)
    {
    psolv->start(p,pgc,f,b->M,b->xvec,b->rhsvec,4);
    return;
    }

    for(int a=p->imin; a<p->imin+p->imax; ++a)
    for(int d=p->jmin; d<p->jmin+p->jmax; ++d)
    if(a<0 || a>=p->knox || d<0 || d>=p->knoy)
    f(a,d) = u(a,d);

    const bool ol = pgc->set_local(true);
    own_solver->start(p,pgc,f,b->M,b->xvec,b->rhsvec,4);
    pgc->set_local(ol);
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

    // curvilinear grids (G 50): the stress of sflow_ediff::curv_stress, its orthogonal part implicit
    if(p->geo_curv==1)
    {
    curv_diff(p,b,pgc,psolv,UHdiff,UH,U,V,WL,alpha,1);
    p->uiter=p->solveriter;
    return;
    }
    
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
    
	solve(p,pgc,psolv,b,U);
    
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

    // curvilinear grids (G 50): the stress of sflow_ediff::curv_stress, its orthogonal part implicit
    if(p->geo_curv==1)
    {
    curv_diff(p,b,pgc,psolv,VHdiff,VH,U,V,WL,alpha,2);
    p->viter=p->solveriter;
    return;
    }
    
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
    
	solve(p,pgc,psolv,b,V);
    
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

        if(p->wet[IJ]==1 && WL(i,j)>0.0 && p->geo_curv==1)
        {
        // curvilinear grids: the orthogonal part of the Laplacian (sflow_curv::lapcoef, L/(A dn) per face) with the
        // (h nu) of the face; no stress through a face to a wall, an open boundary or a dry cell
        double cE,cW,cN,cS;
        sflow_curv::lapcoef(p,i,j,cE,cW,cN,cS);
        const double ih = 1.0/WL(i,j);
        auto fl = [&](int a, int d) { return sflow_dflux::fluid(p,a,d); };

        cxp = fl(i+1,j) ? sflow_dflux::hnu(p,b,WL,i,j,i+1,j)*ih*cE : 0.0;
        cxm = fl(i-1,j) ? sflow_dflux::hnu(p,b,WL,i,j,i-1,j)*ih*cW : 0.0;
        cyp = fl(i,j+1) ? sflow_dflux::hnu(p,b,WL,i,j,i,j+1)*ih*cN : 0.0;
        cym = fl(i,j-1) ? sflow_dflux::hnu(p,b,WL,i,j,i,j-1)*ih*cS : 0.0;
        }

        else
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
    
	solve(p,pgc,psolv,b,W);
    
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

// ---------------------------------------------------------------------------------------------------------------------
// curvilinear grids (G 50): the depth-integrated stress tau = nu h (grad u + grad u^T) of sflow_ediff::curv_stress
// through the faces of the cell over its area, divided by the depth of the cell, with the orthogonal part of the
// normal stress implicit.  Per face (face vector S = n L, cell centres P and N at the distance dn along n):
//   implicit   (h nu)_f L/dn (1 + n_c^2) (f_N - f_P),   n_c = n_x for the x-momentum (comp 1), n_y for comp 2
//   explicit   the rest of the full face flux of curv_stress, evaluated with the velocities of the step:
//              (h nu)_f ( (f_N - f_P) L/dn + k.grad f + transposed part ) - implicit part (deferred correction)
// On a Cartesian grid the implicit part is that of the Cartesian operator (2 (h nu)_f du/dx on the x-faces, 1 on the
// y-faces).  The explicit rest is at most the transposed part of the normal stress, half of the implicit coefficient,
// so the step is stable at any time step (amplification of the shortest mode <= 1/2).  A face to a wall, an open
// boundary or a dry cell carries no stress (free slip), as curv_stress.
void sflow_idiff::curv_diff(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, slice &Qdiff, slice &Q, slice &U, slice &V, slice &WL, double alpha, int comp)
{
    auto fluid = [&](int a, int d) { const int q=(a-p->imin)*p->jmax+(d-p->jmin); return p->flagslice4[q]>0 && p->wet[q]==1; };

    auto cgrad = [&](slice &g, int a, int d, double &gx, double &gy)
    {
        const double gc = g(a,d);
        sflow_curv::grad5(p,a,d, fluid(a+1,d) ? g(a+1,d) : gc, fluid(a-1,d) ? g(a-1,d) : gc,
                                  fluid(a,d+1) ? g(a,d+1) : gc, fluid(a,d-1) ? g(a,d-1) : gc, gx, gy);
    };

    slice &u = (comp==1) ? U : V;
    const int di[4]={1,-1,0,0}, dj[4]={0,0,1,-1};

    n=0;
    SLICELOOP4
    {
    double cf[4] = {0.0,0.0,0.0,0.0};
    double rhs = 0.0;

        if(p->wet[IJ]==1 && WL(i,j)>0.0 && p->flagslice4[IJ]>0)
        {
        const double ihA = 1.0/(WL(i,j)*sflow_curv::area(p,i,j));
        double uxP, uyP, vxP, vyP;
        cgrad(U,i,j,uxP,uyP);
        cgrad(V,i,j,vxP,vyP);

            for(int d=0; d<4; ++d)
            {
                const int a=i+di[d], e=j+dj[d];

                if(!fluid(a,e))
                continue;

                const int dir = (d<2) ? 1 : 2;
                const int fi = (d==1) ? i-1 : i, fj = (d==3) ? j-1 : j;
                const double sg = (d==0 || d==2) ? 1.0 : -1.0;

                double nx, ny, len, kx, ky, uxN, uyN, vxN, vyN;
                sflow_curv::normal(p,dir,fi,fj,nx,ny,len);
                sflow_curv::kvec(p,dir,fi,fj,kx,ky);
                const double dn = sflow_curv::dnorm(p,dir,fi,fj);
                cgrad(U,a,e,uxN,uyN);
                cgrad(V,a,e,vxN,vyN);

                const double Sx = sg*nx*len, Sy = sg*ny*len;
                kx *= sg; ky *= sg;

                const double gx = 0.5*((comp==1 ? uxP : vxP) + (comp==1 ? uxN : vxN));
                const double gy = 0.5*((comp==1 ? uyP : vyP) + (comp==1 ? uyN : vyN));

                const double tr = (comp==1) ? 0.5*(uxP+uxN)*Sx + 0.5*(vxP+vxN)*Sy
                                            : 0.5*(uyP+uyN)*Sx + 0.5*(vyP+vyN)*Sy;

                // (h nu) of the face: the mean of the two viscosities times the smaller depth (sflow_diffusion_flux.h)
                const double nuh = 0.5*(sflow_dflux::nu(p,b,i,j) + sflow_dflux::nu(p,b,a,e))*MIN(MAX(WL(i,j),0.0),MAX(WL(a,e),0.0));

                const double nc = (comp==1) ? nx : ny;
                const double cimp = nuh*len/dn*(1.0 + nc*nc);
                const double du = u(a,e) - u(i,j);
                const double full = nuh*(du*len/dn + gx*kx + gy*ky + tr);

                cf[d] = cimp*ihA;
                rhs += (full - cimp*du)*ihA;
            }
        }

	b->M.p[n] = cf[0] + cf[1] + cf[2] + cf[3] + 1.0/(alpha*p->dt);
    b->M.n[n] = -cf[0];
    b->M.s[n] = -cf[1];
    b->M.w[n] = -cf[2];
    b->M.e[n] = -cf[3];

	b->rhsvec.V[n] = rhs + u(i,j)/(alpha*p->dt);

    f(i,j) = u(i,j);
	++n;
	}

    bc(p,b,u);

	solve(p,pgc,psolv,b,u);

    SLICELOOP4
    Qdiff(i,j) = (p->wet[IJ]==1)?Q(i,j) + WL(i,j)*(f(i,j) - u(i,j)):Q(i,j);

	time=pgc->timer()-starttime;
	if(p->mpirank==0 && p->D21==1 && p->count%p->P12==0)
	cout<<(comp==1 ? "u" : "v")<<"diffiter (curvilinear): "<<p->solveriter<<"  difftime: "<<setprecision(3)<<time<<endl;
}
