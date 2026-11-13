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

#include"rheology_cfd_core.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cmath>

rheology_cfd_core::rheology_cfd_core(lexer *p) : strain(p), rho1(p->W1), gx(p->W20), gy(p->W21), gz(p->W22)
{
    using namespace rheo;
    
    // W 130 1 HB, 2 HBC, 3 mu(I), 4 Coulomb-viscoplastic -> law3D
    const int map[5] = {LAW_NEWTONIAN,LAW_HB,LAW_HBC,LAW_MUI,LAW_COULOMB_VP};
    q.law = (p->W130>=1 && p->W130<=4) ? map[p->W130] : LAW_NEWTONIAN;
    
    q.mu0 = p->W1*p->W2;
    q.tauy = p->W125_tauy;
    q.K = p->W125_K;
    q.n = p->W125_n;
    q.tanphi = tan(deg2rad(p->W121));
    q.m = p->W132;
    q.mu_min = p->W133_min>=0.0 ? p->W133_min : p->W1*p->W2;
    q.mu_max = p->W133_max;
    q.mus = tan(deg2rad(p->W124_ds));
    q.mu2 = tan(deg2rad(p->W124_d2));
    q.I0 = p->W124_I0;
    q.d = p->W124_d;
    q.rho_s = p->W1/MAX(p->W124_phis,1.0e-6);      // grain density from the bulk density W 1
    q.alpha = p->W134;
    
    double I2N = 0.0;
    
    if(q.law==LAW_MUI && q.alpha>0.0)
    I2N = mu_I_setup(q);
    
    if(p->mpirank==0)
    {
    const char *name[5] = {"Newtonian","Herschel-Bulkley","Herschel-Bulkley-Coulomb","mu(I)","Coulomb-viscoplastic"};
    cout<<"CFD rheology core: "<<name[q.law]<<", effective pressure "<<(p->W131==2 ? "solved" : "lithostatic")
        <<", mu_eff "<<q.mu_min<<" - "<<q.mu_max<<" Pa s";
    if(p->W126==1)
    cout<<", pore pressure ratio "<<p->W127_lambda;
    cout<<endl;
    
    if(q.law==LAW_MUI && q.alpha>0.0)
    cout<<"CFD rheology core: mu(I) partially regularised below I_N = "<<q.I_N<<" (alpha "<<q.alpha<<"), ill posed above I = "<<I2N<<endl;
    
    if(p->W126==2)
    cout<<"CFD rheology core: W 126 2 (advected pore pressure) is not available in CFD, no pore pressure"<<endl;
    }
}

double rheology_cfd_core::pressure(lexer *p, fdm *a)
{
    if(p->W131==2)
    return a->press(i,j,k) - p->pressgage;
    
    // lithostatic: rho_1 (g . n) phi, n the unit normal of the level set into phase 1
    const double phival = a->phi(i,j,k);
    
    if(phival<=0.0)
    return 0.0;
    
    double nx = (a->phi(i+1,j,k) - a->phi(i-1,j,k))/(p->DXP[IM1] + p->DXP[IP]);
    double ny = p->j_dir==1 ? (a->phi(i,j+1,k) - a->phi(i,j-1,k))/(p->DYP[JM1] + p->DYP[JP]) : 0.0;
    double nz = (a->phi(i,j,k+1) - a->phi(i,j,k-1))/(p->DZP[KM1] + p->DZP[KP]);
    const double nn = sqrt(nx*nx + ny*ny + nz*nz);
    
    // grad(phi) points into phase 1 (away from the surface): gravity along it
    double gn = fabs(gz);
    
    if(nn>1.0e-10)
    gn = MAX(0.0,(gx*nx + gy*ny + gz*nz)/nn);
    
    return rho1*gn*phival;
}

double rheology_cfd_core::viscosity(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    const double gamma = strainterm(p,u,v,w);
    
    double peff = pressure(p,a);
    
    if(p->W126==1)
    peff *= 1.0 - p->W127_lambda;
    
    return rheo::viscosity(q,gamma,peff)/rho1;
}
