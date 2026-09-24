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

#include"wind_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"slice4.h"

wind_f::wind_f(lexer *p) 
{
    Uref = 31.5; // Zijlema et al. (2012) reference velocity, must be set before Cd is computed
    Cd = 0.0;
    href = 0.0;
    href_set = 0;
    Sw = nullptr;
    Stmp = nullptr;
    Pa = nullptr;
    
    xs = -1.0e10;
    xe =  1.0e10;
    ys = -1.0e10;
    ye =  1.0e10;
    
    if(p->A10==3)
    {
        wind_forcing_drag_coeff_fnpf(p);
        
        if(p->A373==3 || p->A373==4)
        {
        Sw = new slice4(p);
        Stmp = new slice4(p);
        Pa = new slice4(p);
        }
        
        cosa = cos(p->A371_dir*(PI/180.0));
        sina = sin(p->A371_dir*(PI/180.0));
        
        if(p->A373==2 && p->mpirank==0)
        cout<<"A 373 2 is deprecated for FNPF, using A 373 1 (wind setup)"<<endl;
        
        if(p->A372==1)
        {
        xs = p->A372_xs;
        xe = p->A372_xe;
        ys = p->A372_ys;
        ye = p->A372_ye;
        }
    }
    
    
    if(p->A10==5)
    {
        wind_forcing_drag_coeff_nhflow(p);
        
        if(p->A573==3 || p->A573==4)
        {
        Sw = new slice4(p);
        Stmp = new slice4(p);
        Pa = new slice4(p);
        }
        
        cosa = cos(p->A571_dir*(PI/180.0));
        sina = sin(p->A571_dir*(PI/180.0));
        
        if(p->A573==2 && p->mpirank==0)
        cout<<"A 573 2 is deprecated for NHFLOW, using A 573 1 (surface stress)"<<endl;
        
        if(p->A572==1)
        {
        xs = p->A572_xs;
        xe = p->A572_xe;
        ys = p->A572_ys;
        ye = p->A572_ye;
        }
    }

}

wind_f::~wind_f()
{
    delete Sw;
    delete Stmp;
    delete Pa;
}

void wind_f::wind_forcing_ini(lexer *p, ghostcell *pgc)
{
}

// Wind surface pressure for the wave-growth modes, shared by FNPF and NHFLOW. Fills *Pa [Pa].
//
// mode 3 : Miles/Plant, zero-mean pressure in phase with the slope
//          p_a = beta*rho_a*u*^2*deta/dn,  u*^2 = Cd*U10^2,  n = wind direction.
//          Energy growth gamma/omega = beta*(rho_a/rho_w)*(u*/c)^2 for every wave component,
//          beta ~ 32 reproduces Plant (1982), gamma/omega = 0.04*(u*/c)^2.
// mode 4 : modified Jeffreys (Touboul et al. 2006, Kharif et al. 2008):
//          p_a = s*rho_a*(U10-c)^2*deta/dn, only on steep faces |deta/dn| > critical slope sc,
//          switched on smoothly over 0.8..1.0 of sc. c <= 0: celerity from wave generation.
// passes : 1-2-1 low-pass passes on the slope. The input grows like (u*/c)^2, so without the
//          filter unresolved grid-scale waves receive the strongest forcing and grow without bound.
//          Response per pass cos^2(pi*dx/L): 40 cells per wave, 4 passes -> 0.976; 4 cells -> 0.06.
// decay  : 1 (with a forcing area) cos^2 downwind decay over [xs,xe].
void wind_f::wind_pressure(lexer *p, ghostcell *pgc, slice &eta, int mode, double U, double beta, 
                           double s, double sc, double cph, int passes, int decay, int area)
{
    double P,w,psi,dU,slope;
    
    // pressure amplitude per unit slope, p_a = P*deta/dn
    P = 0.0;
    
    if(mode==3)
    P = beta*p->W3*Cd*U*U;
    
    if(mode==4)
    {
    cph = (cph>0.0) ? cph : p->wC;
    dU = MAX(U - cph, 0.0);
    P = s*p->W3*dU*dU;
    }
    
    sc = MAX(sc, 1.0e-6);
    
    // along-wind slope, central differences
    SLICELOOP4
    {
    (*Sw)(i,j) = 0.0;
    
    if(p->wet[IJ]==1)
    {
    (*Sw)(i,j) = cosa*(eta(i+1,j) - eta(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
    
    if(p->j_dir==1)
    (*Sw)(i,j) += sina*(eta(i,j+1) - eta(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
    }
    }
    
    // low-pass filter
    for(int q=0; q<passes; ++q)
    {
        pgc->gcsl_start4(p,*Sw,1);
        
        SLICELOOP4
        (*Stmp)(i,j) = 0.25*(*Sw)(i-1,j) + 0.5*(*Sw)(i,j) + 0.25*(*Sw)(i+1,j);
        
        SLICELOOP4
        (*Sw)(i,j) = (p->wet[IJ]==1) ? (*Stmp)(i,j) : 0.0;
        
        if(p->j_dir==1)
        {
        pgc->gcsl_start4(p,*Sw,1);
        
        SLICELOOP4
        (*Stmp)(i,j) = 0.25*(*Sw)(i,j-1) + 0.5*(*Sw)(i,j) + 0.25*(*Sw)(i,j+1);
        
        SLICELOOP4
        (*Sw)(i,j) = (p->wet[IJ]==1) ? (*Stmp)(i,j) : 0.0;
        }
    }
    
    SLICELOOP4
    {
    (*Pa)(i,j) = 0.0;
    
    if(P>0.0)
    if(p->wet[IJ]==1)
    if( p->XP[IP]>xs && p->XP[IP]<xe)
    if((p->YP[JP]>ys && p->YP[JP]<ye) || p->j_dir==0)
    {
    slope = (*Sw)(i,j);
    
    // Jeffreys: steep faces only, smooth switch-on over 0.8..1.0 of the critical slope
    w = 1.0;
    
    if(mode==4)
    w = MIN(MAX((fabs(slope) - 0.8*sc)/(0.2*sc), 0.0), 1.0);
    
    psi = 1.0;
    
    if(decay==1 && area==1)
    {
    psi = cos(0.5*PI*(p->XP[IP]-xs)/(xe-xs));
    psi = psi*psi;
    }
    
    (*Pa)(i,j) = psi*w*P*slope;
    }
    }
}


