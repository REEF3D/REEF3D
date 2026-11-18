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


#include"rheology_cfd_wall.h"
#include"lexer.h"
#include"fdm.h"
#include"field.h"
#include<cmath>

rheology_cfd_wall::rheology_cfd_wall(lexer *p) : coef(2.0*tan(p->W135_dw*M_PI/180.0)/(p->W1*p->W135_b)),
                                                 lambda(p->W126==1 ? p->W127_lambda : 0.0)
{
    if(p->mpirank==0)
    cout<<"CFD rheology core: side-wall friction delta_w "<<p->W135_dw<<" deg, channel width "<<p->W135_b<<" m (W 135)"<<endl;
}

double rheology_cfd_wall::rate(lexer *p, fdm *a, int c, field &U, field &V, field &W, double &other)
{
    double pr,ph,uu,vv,ww;
    
    if(c==0)
    {
    pr = 0.5*(a->press(i,j,k) + a->press(i+1,j,k));
    ph = 0.5*(a->phi(i,j,k) + a->phi(i+1,j,k));
    uu = U(i,j,k);
    vv = p->j_dir==1 ? 0.25*(V(i,j,k) + V(i+1,j,k) + V(i,j-1,k) + V(i+1,j-1,k)) : 0.0;
    ww = 0.25*(W(i,j,k) + W(i+1,j,k) + W(i,j,k-1) + W(i+1,j,k-1));
    }
    
    else
    {
    pr = 0.5*(a->press(i,j,k) + a->press(i,j,k+1));
    ph = 0.5*(a->phi(i,j,k) + a->phi(i,j,k+1));
    uu = 0.25*(U(i,j,k) + U(i,j,k+1) + U(i-1,j,k) + U(i-1,j,k+1));
    vv = p->j_dir==1 ? 0.25*(V(i,j,k) + V(i,j,k+1) + V(i,j-1,k) + V(i,j-1,k+1)) : 0.0;
    ww = W(i,j,k);
    }
    
    const double peff = (1.0 - lambda)*MAX(0.0, pr - p->pressgage);
    
    if(peff<=0.0 || ph<=-p->psi)
    return 0.0;
    
    double H = 1.0;
    
    if(ph<p->psi)
    H = 0.5*(1.0 + ph/p->psi + sin(M_PI*ph/p->psi)/M_PI);
    
    // the in-plane velocity component that is not updated here (the walls are normal to y)
    other = (c==0) ? ww : uu;
    (void)vv;
    
    return H*coef*peff;
}

void rheology_cfd_wall::implicit(lexer *p, fdm *a, int c, double w, field &f, field &U, field &V, field &W)
{
    // Coulomb step: |u*| reduced by w dt F, stopped where |u*| <= w dt F (the friction holds, no reversal)
    double other,F,speed;
    
    if(c==0)
    ULOOP
    {
    F = rate(p,a,0,U,V,W,other);
    speed = sqrt(f(i,j,k)*f(i,j,k) + other*other);
    
    if(F>0.0)
    f(i,j,k) *= MAX(0.0, 1.0 - w*p->dt*CPOR1*F/MAX(speed,1.0e-12));
    }
    
    if(c==2)
    WLOOP
    {
    F = rate(p,a,2,U,V,W,other);
    speed = sqrt(f(i,j,k)*f(i,j,k) + other*other);
    
    if(F>0.0)
    f(i,j,k) *= MAX(0.0, 1.0 - w*p->dt*CPOR3*F/MAX(speed,1.0e-12));
    }
}
