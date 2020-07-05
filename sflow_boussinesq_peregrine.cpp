/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"sflow_boussinesq_peregrine.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"slice1.h"
#include"slice2.h"

#define DI (fabs(b->depth(i,j))>1.0e-20?b->depth(i,j):1.0e20)
#define DJ (fabs(b->depth(i,j))>1.0e-20?b->depth(i,j):1.0e20)

#define DIP (fabs(b->depth(i+1,j))>1.0e-20?b->depth(i+1,j):1.0e20)
#define DIPP (fabs(b->depth(i+2,j))>1.0e-20?b->depth(i+2,j):1.0e20)
#define DJP (fabs(b->depth(i,j+1))>1.0e-20?b->depth(i,j+1):1.0e20)
#define DJPP (fabs(b->depth(i,j+2))>1.0e-20?b->depth(i,j+2):1.0e20)

#define DIM (fabs(b->depth(i-1,j))>1.0e-20?b->depth(i-1,j):1.0e20)
#define DIMM (fabs(b->depth(i-2,j))>1.0e-20?b->depth(i-2,j):1.0e20)
#define DJM (fabs(b->depth(i,j-1))>1.0e-20?b->depth(i,j-1):1.0e20)
#define DJMM (fabs(b->depth(i,j-2))>1.0e-20?b->depth(i,j-2):1.0e20)

#define DIPJP (fabs(b->depth(i+1,j+1))>1.0e-20?b->depth(i+1,j+1):1.0e20)
#define DIPJM (fabs(b->depth(i+1,j-1))>1.0e-20?b->depth(i+1,j-1):1.0e20)
#define DIMJP (fabs(b->depth(i-1,j+1))>1.0e-20?b->depth(i-1,j+1):1.0e20)
#define DIMJM (fabs(b->depth(i-1,j-1))>1.0e-20?b->depth(i-1,j-1):1.0e20)

#define DX (0.5*(DI + DIP))
#define DY (0.5*(DJ + DJP))

#define DXP (0.5*(DIP + DIPP))
#define DYP (0.5*(DJP + DJPP))

//#define DXM (0.5*(DIM + DIMM))
#define DYM (0.5*(DJM + DJMM))

#define DXY (0.25*(DI + DIP + DJP + DIPJP))

#define DXYM  (0.25*(DI + DIP + DJM + DIPJM))

 
sflow_boussinesq_peregrine::sflow_boussinesq_peregrine(lexer *p, fdm2D *b) : Pxx(p),Qxy(p),Pxx_n(p),Qxy_n(p),
                                                                            Qyy(p),Pxy(p),Qyy_n(p),Pxy_n(p),
                                                                            P1x(p),Q1x(p),Q1y(p),P1x_n(p),Q1x_n(p),Q1y_n(p),
                                                                            P2x(p),P2y(p),Q2x(p),Q2y(p),P2x_n(p),P2y_n(p),Q2x_n(p),Q2y_n(p)
{

}

sflow_boussinesq_peregrine::~sflow_boussinesq_peregrine()
{
}

void sflow_boussinesq_peregrine::ini(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q)
{
    psi1_calc(p,b,pgc,P,Q,b->eta,0.0);
    psi2_calc(p,b,pgc,P,Q,b->eta,0.0);
}

void sflow_boussinesq_peregrine::psi1(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    if(p->count==1)
    ini(p,b,pgc,P,Q);
    
    SLICELOOP1
    {
    Pxx_n(i,j) = Pxx(i,j);
    Qxy_n(i,j) = Qxy(i,j);
    
    P1x_n(i,j) = P1x(i,j);
    Q1x_n(i,j) = Q1x(i,j);
    
    Q1y_n(i,j) = Q1y(i,j);
    }
    
    psi1_calc(p,b,pgc,P,Q,eta,alpha);
    
    // add to source
    SLICELOOP1
    {
    d = 0.5*(b->depth(i+1,j) + b->depth(i,j));
    
    b->F(i,j) -= (1.0/6.0)*pow(d,3.0)*((Pxx(i,j) - Pxx_n(i,j))/(alpha*p->dt) + (Qxy(i,j) - Qxy_n(i,j))/(alpha*p->dt))
    
               - (1.0/2.0)*pow(d,3.0)*((Pxx(i,j) - Pxx_n(i,j))/(alpha*p->dt) + (Qxy(i,j) - Qxy_n(i,j))/(alpha*p->dt));
    }
}

void sflow_boussinesq_peregrine::psi2(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    SLICELOOP2
    {
    Qyy_n(i,j) = Qyy(i,j);
    Pxy_n(i,j) = Pxy(i,j);
    
    Qyy_n(i,j) = Qyy(i,j);
    Pxy_n(i,j) = Pxy(i,j);
    }
    
    psi2_calc(p,b,pgc,P,Q,eta,alpha);
    
    // add to source
    SLICELOOP2
    {
    d = 0.5*(b->depth(i,j) + b->depth(i,j+1));
    
    b->G(i,j) -= (1.0/6.0)*pow(d,3.0)*((Qyy(i,j) - Qyy_n(i,j))/(alpha*p->dt) + (Pxy(i,j) - Pxy_n(i,j))/(alpha*p->dt))
    
                -(1.0/2.0)*pow(d,2.0)*((Qyy(i,j) - Qyy_n(i,j))/(alpha*p->dt) + (Pxy(i,j) - Pxy_n(i,j))/(alpha*p->dt));
    }
}

void sflow_boussinesq_peregrine::psi1_calc(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    SLICELOOP1
    {
    Pxx(i,j) = (P(i+1,j) - 2.0*P(i,j) + P(i-1,j))/(p->DXM*p->DXM);
    Qxy(i,j) = (Q(i+1,j) - Qxy(i,j) - Q(i+1,j-1) + Qxy(i,j-1)) /(p->DXM*p->DXM);
    
    //cout<<Pxx(i,j)<<endl;
    
    P1x(i,j) = (P(i+1,j) - P(i-1,j))/(2.0*p->DXM);
    Q1x(i,j) = (0.5*(Q(i+1,j)+Q(i+1,j-1)) - 0.5*(Qxy(i,j)  + Qxy(i,j-1)))/p->DXM;
    Q1y(i,j) = (0.5*(Q(i+1,j)+Qxy(i,j)) - 0.5*(Q(i+1,j-1) + Qxy(i,j-1)))/p->DXM;
    }
    
    
     
}

void sflow_boussinesq_peregrine::psi2_calc(lexer *p, fdm2D *b, ghostcell *pgc, slice &P, slice &Q, slice &eta, double alpha)
{
    SLICELOOP2
    {

    
    Qyy(i,j) = (Q(i,j+1) - 2.0*Q(i,j) + P(i,j-1))/(p->DXM*p->DXM);
    Pxy(i,j) = (P(i+1,j) - Pxy(i,j) - P(i+1,j-1) + Pxy(i,j-1)) /(p->DXM*p->DXM);
    }
}