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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"momentum_IMEX.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"bcmom.h"
#include"convection.h"
#include"diffusion.h"
#include"idiff_IMEX.h"
#include"idiff_IMEX_2D.h"
#include"pjm_IMEX.h"
#include"poisson.h"
#include"ioflow.h"
#include"turbulence.h"
#include"solver.h"
#include"fluid_update_rheology.h"
#include"fluid_update_void.h"

momentum_IMEX::momentum_IMEX(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson, turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow):bcmom(p),un(p),vn(p),wn(p),F0_ex(p),G0_ex(p),H0_ex(p),F1_ex(p),G1_ex(p),H1_ex(p),F2_ex(p),G2_ex(p),H2_ex(p),F1_im(p),G1_im(p),H1_im(p),F2_im(p),G2_im(p),H2_im(p)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
	
    // 2nd order scheme
    gamma = 1.0/6.0*(3.0 + sqrt(3.0));
    a11 = 0.5;
    a21 = 0.0;
    a22 = 0.0;
    ahat10 = 0.5;
    ahat20 = 0.0;
    ahat21 = 0.0;
    b1 = 1.0;
    b2 = 0.0;
    bhat1 = 1.0;
    bhat2 = 0.0;

    // 3rd order scheme
    gamma = 1.0/6.0*(3.0 + sqrt(3.0));
    a11 = gamma;
    a21 = 1.0 - 2.0*gamma;
    a22 = gamma;
    ahat10 = gamma;
    ahat20 = gamma - 1.0;
    ahat21 = 2.0*(1.0 - gamma);
    b1 = 0.5;
    b2 = 0.5;
    bhat1 = 0.5;
    bhat2 = 0.5;


	pconvec=pconvection;
	ppress=ppressure;
    pdiff=pdiffusion;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;
    
    if(p->W90==0  || p->F300>0)
	pupdate = new fluid_update_void();
    
    if(p->W90==1 && p->F300==0)
	pupdate = new fluid_update_rheology(p,a);
}

momentum_IMEX::~momentum_IMEX()
{
}

void momentum_IMEX::start(lexer *p, fdm* a, ghostcell* pgc, vrans *pvrans)
{
    // Apply in- and outflow conditions
	pflow->discharge(p,a,pgc);
	pflow->inflow(p,a,pgc,a->u,a->v,a->w);


    // Store old velocities
	ULOOP
	un(i,j,k)=a->u(i,j,k);
	
    VLOOP
	vn(i,j,k)=a->v(i,j,k);
    
    WLOOP
	wn(i,j,k)=a->w(i,j,k);

	pgc->start1(p,un,gcval_u);
    pgc->start2(p,vn,gcval_v);
	pgc->start3(p,wn,gcval_w);
	
    
    //- First RK step

    // Calculate explicit terms
    eval_ex_F(p, a, pgc, pvrans, F0_ex, a->u, a->v, a->w);
    eval_ex_G(p, a, pgc, pvrans, G0_ex, a->u, a->v, a->w);
    eval_ex_H(p, a, pgc, pvrans, H0_ex, a->u, a->v, a->w);
    
    // Calculate constant rhs for implicit step

    ULOOP
    {
        a->F(i,j,k) = ahat10/a11*F0_ex(i,j,k) + a->gi + (CPOR1*un(i,j,k))/(a11*p->dt);
    }
    VLOOP
    {
        a->G(i,j,k) = ahat10/a11*G0_ex(i,j,k) + a->gj + (CPOR2*vn(i,j,k))/(a11*p->dt);
    }
    WLOOP
    {
        a->H(i,j,k) = ahat10/a11*H0_ex(i,j,k) + a->gk + (CPOR3*wn(i,j,k))/(a11*p->dt);
    }

    for (int it_rk = 0; it_rk < 2; it_rk++)
    {
        // Evaluate momentum predictor step in a->u,v,w 
        pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,a11);                  
        pgc->start1(p,a->u,gcval_u);
        pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,a11);
        pgc->start2(p,a->v,gcval_v);
        pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,a11);
        pgc->start3(p,a->w,gcval_w);

        // Evaluate pressure in a->press
        pflow->pressure_io(p,a,pgc);
	    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,a->u,a->v,a->w,a11);
	    pflow->p_relax(p,a,pgc,a->press);
    }

    // Divergence free velocity field
    ppress->ucorr(p,a,a->u,a11);
    ppress->vcorr(p,a,a->v,a11);
    ppress->wcorr(p,a,a->w,a11);
    
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	
	pgc->start1(p,a->u,gcval_u);
	pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);
        
    pupdate->start(p,a,pgc);


    // Calculate explicit terms
    eval_ex_F(p, a, pgc, pvrans, F1_ex, a->u, a->v, a->w);
    eval_ex_G(p, a, pgc, pvrans, G1_ex, a->u, a->v, a->w);
    eval_ex_H(p, a, pgc, pvrans, H1_ex, a->u, a->v, a->w);

    // Calculate implicit terms
	ULOOP
    {
        F1_im(i,j,k) = (a->u(i,j,k) - un(i,j,k))/(a11*p->dt) - ahat10/a11*F0_ex(i,j,k);
    }
	VLOOP
    {
        G1_im(i,j,k) = (a->v(i,j,k) - vn(i,j,k))/(a11*p->dt) - ahat10/a11*G0_ex(i,j,k); 
    }
	WLOOP
    {
        H1_im(i,j,k) = (a->w(i,j,k) - wn(i,j,k))/(a11*p->dt) - ahat10/a11*H0_ex(i,j,k); 
    }
    
    
    // Recombination step 2nd order
    ULOOP
    {
        a->u(i,j,k) = un(i,j,k) + p->dt*(b1*F1_im(i,j,k) + bhat1*F1_ex(i,j,k));
    }
    VLOOP
    {
        a->v(i,j,k) = vn(i,j,k) + p->dt*(b1*G1_im(i,j,k) + bhat1*G1_ex(i,j,k));
    }
    WLOOP
    {
        a->w(i,j,k) = wn(i,j,k) + p->dt*(b1*H1_im(i,j,k) + bhat1*H1_ex(i,j,k));
    }

    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);

    pflow->pressure_io(p,a,pgc);
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,a->u,a->v,a->w,b1);
    pflow->p_relax(p,a,pgc,a->press);
    ppress->ucorr(p,a,a->u,b1);
    ppress->vcorr(p,a,a->v,b1);
    ppress->wcorr(p,a,a->w,b1);

    pflow->u_relax(p,a,pgc,a->u);
    pflow->v_relax(p,a,pgc,a->v);
    pflow->w_relax(p,a,pgc,a->w);
    pflow->p_relax(p,a,pgc,a->press);

    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);
        
    pupdate->start(p,a,pgc);


    // 3rd-order version
    if (a22 != 0.0)
    {
        // Calculate constant rhs for implicit step
        ULOOP
        {
            a->F(i,j,k) = a21/a22*F1_im(i,j,k) + ahat20/a22*F0_ex(i,j,k) + ahat21/a22*F1_ex(i,j,k) + a->gi + (CPOR1*un(i,j,k))/(a22*p->dt);
        }
        VLOOP
        {
            a->G(i,j,k) = a21/a22*G1_im(i,j,k) + ahat20/a22*G0_ex(i,j,k) + ahat21/a22*G1_ex(i,j,k) + a->gj + (CPOR2*vn(i,j,k))/(a22*p->dt);
        }
        WLOOP
        {
            a->H(i,j,k) = a21/a22*H1_im(i,j,k) + ahat20/a22*H0_ex(i,j,k) + ahat21/a22*H1_ex(i,j,k) + a->gk + (CPOR3*wn(i,j,k))/(a22*p->dt);
        }
       
        for (int it_rk = 0; it_rk < 2; it_rk++)
        {
            // Evaluate momentum predictor step in a->u,v,w 
            pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,a22);     
            pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,a22);
            pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,a22);

            // Evaluate pressure in a->press
            pflow->pressure_io(p,a,pgc);
            ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,a->u,a->v,a->w,a22);
            pflow->p_relax(p,a,pgc,a->press);
        }
        
        // Divergence free velocity field
        ppress->ucorr(p,a,a->u,a22);
        ppress->vcorr(p,a,a->v,a22);
        ppress->wcorr(p,a,a->w,a22);

        pflow->u_relax(p,a,pgc,a->u);
        pflow->v_relax(p,a,pgc,a->v);
        pflow->w_relax(p,a,pgc,a->w);
        
        pgc->start1(p,a->u,gcval_u);
        pgc->start2(p,a->v,gcval_v);
        pgc->start3(p,a->w,gcval_w);
        
        pupdate->start(p,a,pgc);

        // Calculate explicit terms
        eval_ex_F(p, a, pgc, pvrans, F2_ex, a->u, a->v, a->w);
        eval_ex_G(p, a, pgc, pvrans, G2_ex, a->u, a->v, a->w);
        eval_ex_H(p, a, pgc, pvrans, H2_ex, a->u, a->v, a->w);

        // Calculate implicit terms
        ULOOP
        {
            F2_im(i,j,k) = (a->u(i,j,k) - un(i,j,k))/(a22*p->dt) - ahat20/a22*F0_ex(i,j,k) - ahat21/a22*F1_ex(i,j,k) - a21/a22*F1_im(i,j,k);
        }
        VLOOP
        {
            G2_im(i,j,k) = (a->v(i,j,k) - vn(i,j,k))/(a22*p->dt) - ahat20/a22*G0_ex(i,j,k) - ahat21/a22*G1_ex(i,j,k) - a21/a22*G1_im(i,j,k); 
        }
        WLOOP
        {
            H2_im(i,j,k) = (a->w(i,j,k) - wn(i,j,k))/(a22*p->dt) - ahat20/a22*H0_ex(i,j,k) - ahat21/a22*H1_ex(i,j,k) - a21/a22*H1_im(i,j,k); 
        }


        //- Recombination step

        // Calculate velocity
        ULOOP
        {
            a->u(i,j,k) = un(i,j,k) + p->dt*(b1*F1_im(i,j,k) + b2*F2_im(i,j,k) + bhat1*F1_ex(i,j,k) + bhat2*F2_ex(i,j,k)); 
        }
        VLOOP
        {
            a->v(i,j,k) = vn(i,j,k) + p->dt*(b1*G1_im(i,j,k) + b2*G2_im(i,j,k) + bhat1*G1_ex(i,j,k) + bhat2*G2_ex(i,j,k)); 
        }
        WLOOP
        {
            a->w(i,j,k) = wn(i,j,k) + p->dt*(b1*H1_im(i,j,k) + b2*H2_im(i,j,k) + bhat1*H1_ex(i,j,k) + bhat2*H2_ex(i,j,k)); 
        }
        
        pgc->start1(p,a->u,gcval_u);
        pgc->start2(p,a->v,gcval_v);
        pgc->start3(p,a->w,gcval_w);

        // Calculate pressure
        pflow->pressure_io(p,a,pgc);
        ppress->start(a,p,ppois,ppoissonsolv,pgc,pflow,a->u,a->v,a->w,b2);
        pflow->p_relax(p,a,pgc,a->press);
        ppress->ucorr(p,a,a->u,b2);
        ppress->vcorr(p,a,a->v,b2);
        ppress->wcorr(p,a,a->w,b2);
        
        pflow->u_relax(p,a,pgc,a->u);
        pflow->v_relax(p,a,pgc,a->v);
        pflow->w_relax(p,a,pgc,a->w);
        pflow->p_relax(p,a,pgc,a->press);

        pgc->start1(p,a->u,gcval_u);
        pgc->start2(p,a->v,gcval_v);
        pgc->start3(p,a->w,gcval_w);
    
        pupdate->start(p,a,pgc);
    }
}

void momentum_IMEX::eval_ex_F(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, field &f, field &uvel, field &vvel, field &wvel)
{
    // Evaluate explicit term in F and store in f
	pturb->isource(p,a);
	pflow->isource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,uvel,gcval_u);
	pconvec->start(p,a,uvel,1,uvel,vvel,wvel);

    ULOOP
    {
        f(i,j,k) = a->F(i,j,k);
        
        a->F(i,j,k) = 0.0;
    }
}

void momentum_IMEX::eval_ex_G(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, field &f, field &uvel, field &vvel, field &wvel)
{
    // Evaluate explicit term in G and store in f
	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,vvel,gcval_v);
	pconvec->start(p,a,vvel,2,uvel,vvel,wvel);

    VLOOP
    {
        f(i,j,k) = a->G(i,j,k);
        
        a->G(i,j,k) = 0.0;
    }
}

void momentum_IMEX::eval_ex_H(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans, field &f, field &uvel, field &vvel, field &wvel)
{
    // Evaluate explicit term in H and store in f
	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc,pvrans);
	bcmom_start(a,p,pgc,pturb,wvel,gcval_w);
	pconvec->start(p,a,wvel,3,uvel,vvel,wvel);

    WLOOP
    {
        f(i,j,k) = a->H(i,j,k) ;
        
        a->H(i,j,k) = 0.0;
    }
}

void momentum_IMEX::utimesave(lexer *p, fdm *a, ghostcell* pgc)
{
}

void momentum_IMEX::vtimesave(lexer *p, fdm *a, ghostcell* pgc)
{
}

void momentum_IMEX::wtimesave(lexer *p, fdm *a, ghostcell* pgc)
{
}
