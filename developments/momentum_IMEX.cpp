/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

momentum_IMEX::momentum_IMEX(lexer *p, fdm *a, convection *pconvection, diffusion *pdiffusion, pressure* ppressure, poisson* ppoisson,
                                                    turbulence *pturbulence, solver *psolver, solver *ppoissonsolver, ioflow *pioflow)
                                                    :bcmom(p),un(p),vn(p),wn(p),F0_ex(p),G0_ex(p),H0_ex(p),F1_ex(p),G1_ex(p),H1_ex(p),F2_ex(p),G2_ex(p),H2_ex(p),F1_im(p),G1_im(p),H1_im(p),F2_im(p),G2_im(p),H2_im(p), u_theo(p), Fu_theo(p), v_theo(p), Fv_theo(p), p_theo(p)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
	
    gcval_urk=20;
	gcval_vrk=21;
	gcval_wrk=22;

    // 1st order scheme
/*    gamma = 1.0/6.0*(3.0 + sqrt(3.0));
    a11 = 1.0;
    a21 = 0.0;
    a22 = 0.0;
    ahat10 = 1.0;
    ahat20 = 0.0;
    ahat21 = 0.0;
    b1 = 1.0;
    b2 = 0.0;
    bhat1 = 1.0;
    bhat2 = 0.0;
*/

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

/*
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
*/
	pconvec=pconvection;
	pdiff=pdiffusion;
	ppois=ppoisson;
	pturb=pturbulence;
	psolv=psolver;
    ppoissonsolv=ppoissonsolver;
	pflow=pioflow;

	ppress = new pjm_IMEX(p,a);
	
    if(p->j_dir==1)
    {
        twoD = 1.0;
        pdiff = new idiff_IMEX(p);
    }
    
    if(p->j_dir==0)
    {
        twoD = 0.0;
        pdiff = new idiff_IMEX_2D(p);
    }
}

momentum_IMEX::~momentum_IMEX()
{
}

void momentum_IMEX::start(lexer *p, fdm* a, ghostcell* pgc, momentum *pmom)
{
    // Apply in- and outflow conditions
	pflow->discharge(p,a,pgc);
	pflow->inflow(p,a,pgc,a->u,a->v,a->w);


// Taylor Green vortex
/*double t = 0.0;
ULOOP
{
     a->u(i,j,k) = exp(-2.0*p->W2*t)*cos(p->pos1_x())*sin(p->pos1_z());
}
WLOOP
{
     a->w(i,j,k) = -exp(-2.0*p->W2*t)*sin(p->pos3_x())*cos(p->pos3_z());
}
LOOP
{
     a->press(i,j,k) = -0.25*exp(-2.0*p->W2*t)*exp(-2.0*p->W2*t)*(cos(2.0*p->pos_x()) + cos(2.0*p->pos_z()));
}

p->dt = 0.001;
int nIt = 10;
for (int tt=0; tt<nIt; tt++)
{    
t += p->dt;
*/



// Analytical solution
//double pi = PI;
//double nu = p->W2;

//double t = 0.0;
//p->dt = 0.001; //0.0625;
//double nIt = 10;// 1.0/p->dt;

//for (int tt=0; tt<nIt; tt++)
//{    
//t += p->dt;
//ULOOP
//{
//    double x = p->pos1_x();
//    double y = p->pos1_z();
///*
//    Fu_theo(i,j,k) = pi * cos(t) * sin((2 * pi * y)) * pow(sin((pi * x)), 0.2e1) + 0.2e1 * pow((double) pi, (double) 3) * pow(sin(t), 0.2e1) * pow(sin((double) (2 * pi * y)), 0.2e1) * pow(sin((double) (pi * x)), 0.3e1) * cos((double) (pi * x)) - 0.2e1 * (double) pow((double) pi, (double) 3) * pow(sin(t), 0.2e1) * sin((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1) * cos((double) (2 * pi * y)) * pow(sin((double) (pi * x)), 0.2e1) - 0.2e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * y)) * pow(cos((double) (pi * x)), 0.2e1) + 0.2e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * y)) * pow(sin((double) (pi * x)), 0.2e1) - sin(t) * (double) pi * sin((double) (pi * x)) * sin((double) (pi * y));
//*/
///*
//    // Stokes problem
//    Fu_theo(i,j,k) = (double) pi * cos(t) * sin((double) (2 * pi * y)) * pow(sin((double) (pi * x)), 0.2e1) - 0.2e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * y)) * pow(cos((double) (pi * x)), 0.2e1) + 0.2e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * y)) * pow(sin((double) (pi * x)), 0.2e1) - sin(t) * (double) pi * sin((double) (pi * x)) * sin((double) (pi * y));
//*/
//    Fu_theo(i,j,k) = (double) pi * cos(t) * sin((double) (2 * pi * y)) * pow(sin((double) (pi * x)), 0.2e1) - 0.4e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * y)) * pow(cos((double) (pi * x)), 0.2e1) + 0.4e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * y)) * pow(sin((double) (pi * x)), 0.2e1) - nu * (-0.4e1 * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * y)) * pow(sin((double) (pi * x)), 0.2e1) - 0.4e1 * (double) pow((double) pi, (double) 3) * sin(t) * cos((double) (2 * pi * x)) * sin((double) (pi * y)) * cos((double) (pi * y))) - sin(t) * (double) pi * sin((double) (pi * x)) * sin((double) (pi * y));
//
//}
//
//WLOOP
//{
//    double x = p->pos3_x();
//    double y = p->pos3_z();
///*
//    Fv_theo(i,j,k) = -(double) pi * cos(t) * sin((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1) - 0.2e1 * (double) pow((double) pi, (double) 3) * pow(sin(t), 0.2e1) * sin((double) (2 * pi * y)) * pow(sin((double) (pi * x)), 0.2e1) * cos((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1) + 0.2e1 * (double) pow((double) pi, (double) 3) * pow(sin(t), 0.2e1) * pow(sin((double) (2 * pi * x)), 0.2e1) * pow(sin((double) (pi * y)), 0.3e1) * cos((double) (pi * y)) + 0.2e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * x)) * pow(cos((double) (pi * y)), 0.2e1) - 0.2e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1) + sin(t) * cos((double) (pi * x)) * (double) pi * cos((double) (pi * y));
//  */  
//   /* 
//    // Stokes problem
//    Fv_theo(i,j,k) = -(double) pi * cos(t) * sin((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1) + 0.2e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * x)) * pow(cos((double) (pi * y)), 0.2e1) - 0.2e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1) + sin(t) * cos((double) (pi * x)) * (double) pi * cos((double) (pi * y));
//*/
//    Fv_theo(i,j,k) = -(double) pi * cos(t) * sin((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1) + 0.4e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * x)) * pow(cos((double) (pi * y)), 0.2e1) - 0.4e1 * nu * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1) - nu * (0.4e1 * (double) pow((double) pi, (double) 3) * sin(t) * cos((double) (2 * pi * y)) * sin((double) (pi * x)) * cos((double) (pi * x)) + 0.4e1 * (double) pow((double) pi, (double) 3) * sin(t) * sin((double) (2 * pi * x)) * pow(sin((double) (pi * y)), 0.2e1)) + sin(t) * cos((double) (pi * x)) * (double) pi * cos((double) (pi * y));
//}
//    


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
	
	//starttime=pgc->timer();


	
    //- First RK step

    // Calculate explicit terms
    eval_ex_F(p, a, pgc, F0_ex, a->u, a->v, a->w);
    eval_ex_G(p, a, pgc, G0_ex, a->u, a->v, a->w);
    eval_ex_H(p, a, pgc, H0_ex, a->u, a->v, a->w);
    
    // Calculate constant rhs for implicit step

    ULOOP
    {
        a->F(i,j,k) = ahat10/a11*F0_ex(i,j,k) + a->gi + (CPOR1*un(i,j,k))/(a11*p->dt) ;//+ Fu_theo(i,j,k);
    }
    VLOOP
    {
        a->G(i,j,k) = ahat10/a11*G0_ex(i,j,k) + a->gj + (CPOR2*vn(i,j,k))/(a11*p->dt);
    }
    WLOOP
    {
        a->H(i,j,k) = ahat10/a11*H0_ex(i,j,k) + a->gk + (CPOR3*wn(i,j,k))/(a11*p->dt) ;//+ Fv_theo(i,j,k);
    }

    for (int it_rk = 0; it_rk < 10; it_rk++)
    {
        // Evaluate momentum predictor step in a->u,v,w 
        pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,a11);                  
        pgc->start1(p,a->u,gcval_urk);
        pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,a11);
        pgc->start2(p,a->v,gcval_vrk);
        pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,a11);
        pgc->start3(p,a->w,gcval_wrk);

        // Evaluate pressure in a->press
        pflow->pressure_io(p,a,pgc);
	    ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow,a->u,a->v,a->w,a11);
	    pflow->p_relax(p,a,pgc,a->press);
    }

    // Divergence free velocity field
    ppress->ucorr(a,p,a->u,a11);
    ppress->vcorr(a,p,a->v,a11);
    ppress->wcorr(a,p,a->w,a11);
    
	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	
	pgc->start1(p,a->u,gcval_urk);
	pgc->start2(p,a->v,gcval_vrk);
	pgc->start3(p,a->w,gcval_wrk);


    // Calculate explicit terms
    eval_ex_F(p, a, pgc, F1_ex, a->u, a->v, a->w);
    eval_ex_G(p, a, pgc, G1_ex, a->u, a->v, a->w);
    eval_ex_H(p, a, pgc, H1_ex, a->u, a->v, a->w);

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
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow,a->u,a->v,a->w,b1);
    pflow->p_relax(p,a,pgc,a->press);
    ppress->ucorr(a,p,a->u,b1);
    ppress->vcorr(a,p,a->v,b1);
    ppress->wcorr(a,p,a->w,b1);

    pflow->u_relax(p,a,pgc,a->u);
    pflow->v_relax(p,a,pgc,a->v);
    pflow->w_relax(p,a,pgc,a->w);
    pflow->p_relax(p,a,pgc,a->press);

    pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
    pgc->start3(p,a->w,gcval_w);


/*    
    // Calculate constant rhs for implicit step
    ULOOP
    {
        a->F(i,j,k) = a21/a22*F1_im(i,j,k) + ahat20/a22*F0_ex(i,j,k) + ahat21/a22*F1_ex(i,j,k) + a->gi + (CPOR1*un(i,j,k))/(a22*p->dt);// + Fu_theo(i,j,k);
    }
    VLOOP
    {
        a->G(i,j,k) = a21/a22*G1_im(i,j,k) + ahat20/a22*G0_ex(i,j,k) + ahat21/a22*G1_ex(i,j,k) + a->gj + (CPOR2*vn(i,j,k))/(a22*p->dt);
    }
    WLOOP
    {
        a->H(i,j,k) = a21/a22*H1_im(i,j,k) + ahat20/a22*H0_ex(i,j,k) + ahat21/a22*H1_ex(i,j,k) + a->gk + (CPOR3*wn(i,j,k))/(a22*p->dt);// + Fv_theo(i,j,k);
    }
   
    for (int it_rk = 0; it_rk < 10; it_rk++)
    {
        // Evaluate momentum predictor step in a->u,v,w 
        pdiff->diff_u(p,a,pgc,psolv,a->u,a->v,a->w,a22);     
        pdiff->diff_v(p,a,pgc,psolv,a->u,a->v,a->w,a22);
        pdiff->diff_w(p,a,pgc,psolv,a->u,a->v,a->w,a22);

        // Evaluate pressure in a->press
        pflow->pressure_io(p,a,pgc);
	    ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow,a->u,a->v,a->w,a22);
	    pflow->p_relax(p,a,pgc,a->press);
    }
    
    // Divergence free velocity field
    ppress->ucorr(a,p,a->u,a22);
    ppress->vcorr(a,p,a->v,a22);
    ppress->wcorr(a,p,a->w,a22);

	pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	
	pgc->start1(p,a->u,gcval_urk);
	pgc->start2(p,a->v,gcval_vrk);
	pgc->start3(p,a->w,gcval_wrk);

    // Calculate explicit terms
    eval_ex_F(p, a, pgc, F2_ex, a->u, a->v, a->w);
    eval_ex_G(p, a, pgc, G2_ex, a->u, a->v, a->w);
    eval_ex_H(p, a, pgc, H2_ex, a->u, a->v, a->w);

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
    ppress->start(a,p,ppois,ppoissonsolv,pgc,pmom,pflow,a->u,a->v,a->w,b2);
    pflow->p_relax(p,a,pgc,a->press);
    ppress->ucorr(a,p,a->u,b2);
    ppress->vcorr(a,p,a->v,b2);
    ppress->wcorr(a,p,a->w,b2);
	
    pflow->u_relax(p,a,pgc,a->u);
	pflow->v_relax(p,a,pgc,a->v);
	pflow->w_relax(p,a,pgc,a->w);
	pflow->p_relax(p,a,pgc,a->press);

	pgc->start1(p,a->u,gcval_u);
    pgc->start2(p,a->v,gcval_v);
	pgc->start3(p,a->w,gcval_w);

	pflow->periodic(a->u,p);
	pflow->periodic(a->v,p);
	pflow->periodic(a->w,p);
*/


    

//Exact solution vortex
/*ULOOP
{
     u_theo(i,j,k) = exp(-2.0*p->W2*t)*cos(p->pos1_x())*sin(p->pos1_z());
}
WLOOP
{
     v_theo(i,j,k) = -exp(-2.0*p->W2*t)*sin(p->pos3_x())*cos(p->pos3_z());
}
LOOP
{
     p_theo(i,j,k) = -0.25*exp(-2.0*p->W2*t)*exp(-2.0*p->W2*t)*(cos(2.0*p->pos_x()) + cos(2.0*p->pos_z()));
     a->test(i,j,k) = p_theo(i,j,k);
}*/

/*
ULOOP
{
    double x = p->pos1_x();
    double y = p->pos1_z();

    u_theo(i,j,k) = (double) pi * sin(t) * sin((double) (2 * pi * y)) * sin((double) (pi * x)) * sin((double) (pi * x));
}
WLOOP
{
    double x = p->pos3_x();
    double y = p->pos3_z();
    
    v_theo(i,j,k) = (-0.1e1) * (double) pi * sin(t) * sin((double) (2 * pi * x)) * sin((double) (pi * y)) * sin((double) (pi * y));
}
LOOP
{
    double x = p->pos_x();
    double y = p->pos_z();
    p_theo(i,j,k) = sin(t) * cos(pi * x) * sin(pi * y);

    double u = (double) pi * sin(t) * sin((double) (2 * pi * y)) * sin((double) (pi * x)) * sin((double) (pi * x));
    double v = (-0.1e1) * (double) pi * sin(t) * sin((double) (2 * pi * x)) * sin((double) (pi * y)) * sin((double) (pi * y));
    
    a->test(i,j,k) = p_theo(i,j,k);
}*/
/*
}

        ofstream result;
        result.open("./REEF3D-Analytical-U.csv", ios::binary);
        ULOOP
        {
            result<<p->pos1_x()<<","<<p->pos1_y()<<","<<p->pos1_z()<<","<<a->u(i,j,k)<<","<<u_theo(i,j,k)<<endl;
        }
        result.close();
        result.open("./REEF3D-Analytical-V.csv", ios::binary);
        WLOOP
        {
            result<<p->pos3_x()<<","<<p->pos3_y()<<","<<p->pos3_z()<<","<<a->w(i,j,k)<<","<<v_theo(i,j,k)<<endl;
        }
        result.close();
        result.open("./REEF3D-Analytical-P.csv", ios::binary);
        LOOP
        {
            result<<p->pos_x()<<","<<p->pos_y()<<","<<p->pos_z()<<","<<a->press(i,j,k)<<","<<p_theo(i,j,k)<<endl;
        }
        result.close();
*/

}

void momentum_IMEX::eval_ex_F(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel)
{
    // Evaluate explicit term in F and store in f
	pturb->isource(p,a);
	pflow->isource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,uvel,gcval_u);
	pconvec->start(p,a,uvel,1,uvel,vvel,wvel);

    ULOOP
    {
        f(i,j,k) = a->F(i,j,k) ;//+ Fu_theo(i,j,k);
        
        a->F(i,j,k) = 0.0;
    }
}

void momentum_IMEX::eval_ex_G(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel)
{
    // Evaluate explicit term in G and store in f
	pturb->jsource(p,a);
	pflow->jsource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,vvel,gcval_v);
	pconvec->start(p,a,vvel,2,uvel,vvel,wvel);

    VLOOP
    {
        f(i,j,k) = a->G(i,j,k);
        
        a->G(i,j,k) = 0.0;
    }
}

void momentum_IMEX::eval_ex_H(lexer *p, fdm *a, ghostcell *pgc, field &f, field &uvel, field &vvel, field &wvel)
{
    // Evaluate explicit term in H and store in f
	pturb->ksource(p,a);
	pflow->ksource(p,a,pgc);
	bcmom_start(a,p,pgc,pturb,wvel,gcval_w);
	pconvec->start(p,a,wvel,3,uvel,vvel,wvel);

    WLOOP
    {
        f(i,j,k) = a->H(i,j,k) ;//+ Fv_theo(i,j,k);
        
        a->H(i,j,k) = 0.0;
    }
}

void momentum_IMEX::irhs(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
}

void momentum_IMEX::jrhs(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
}

void momentum_IMEX::krhs(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
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

void momentum_IMEX::fillaij1(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

void momentum_IMEX::fillaij2(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}

void momentum_IMEX::fillaij3(lexer *p, fdm *a, ghostcell* pgc, solver *psolv)
{
}


void::momentum_IMEX::eval_im_F
(
    lexer *p, 
    fdm *a, 
    ghostcell *pgc, 
    field& u, field& v, field& w,
    field& f
)
{
}


void::momentum_IMEX::eval_im_G
(
    lexer *p, 
    fdm *a, 
    ghostcell *pgc, 
    field& u, field& v, field& w,
    field& f
)
{
}


void::momentum_IMEX::eval_im_H
(
    lexer *p, 
    fdm *a, 
    ghostcell *pgc, 
    field& u, field& v, field& w,
    field& f
)
{
}

