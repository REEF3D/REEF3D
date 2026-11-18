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
#include"rheology_cfd_3phase.h"
#include"interface_width.h"
#include<cmath>

rheology_cfd_core::rheology_cfd_core(lexer *p) : strain(p), rho1(p->W1), gx(p->W20), gy(p->W21), gz(p->W22)
{
    using namespace rheo;
    
    // W 130 1 HB, 2 HBC, 3 mu(I), 4 Coulomb-viscoplastic, 5 mu(J) -> law3D
    const int map[6] = {LAW_NEWTONIAN,LAW_HB,LAW_HBC,LAW_MUI,LAW_COULOMB_VP,LAW_MUJ};
    q.law = (p->W130>=1 && p->W130<=5) ? map[p->W130] : LAW_NEWTONIAN;
    
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
    
    // three-phase model: the grain density W 142; mu(J): the interstitial fluid is the water (three-phase) or phase 2
    if(p->W140>0)
    q.rho_s = p->W142;
    
    q.J0 = p->W146_J0;
    q.aK = p->W146_aK;
    q.phic = p->W145_phic;
    q.eta_f = p->W140>0 ? p->W1*p->W2 : p->W3*p->W4;
    
    double I2N = 0.0;
    
    if(q.law==LAW_MUI && q.alpha>0.0)
    I2N = mu_I_setup(q);
    
    if(p->mpirank==0)
    {
    const char *name[6] = {"Newtonian","Herschel-Bulkley","Herschel-Bulkley-Coulomb","mu(I)","Coulomb-viscoplastic","mu(J)"};
    cout<<"CFD rheology core: "<<name[q.law]<<", effective pressure "<<(p->W131==2 ? "solved" : (p->W131==3 ? "solved, at least the weight of the interface band" : "lithostatic"))
        <<", mu_eff "<<q.mu_min<<" - "<<q.mu_max<<" Pa s";
    if(p->W126==1)
    cout<<", pore pressure ratio "<<p->W127_lambda;
    cout<<endl;
    
    if(q.law==LAW_MUJ)
    cout<<"CFD rheology core: mu(J) J0 "<<q.J0<<", K = J + "<<q.aK<<" I^2, eta_f "<<q.eta_f<<" Pa s"<<endl;
    
    if(q.law==LAW_MUI && q.alpha>0.0)
    {
    if(q.I_N>0.0)
    cout<<"CFD rheology core: mu(I) partially regularised below I_N = "<<q.I_N<<" (alpha "<<q.alpha<<")"<<endl;
    else
    cout<<"CFD rheology core: mu(I) well posed at low I, no regularisation (alpha "<<q.alpha<<")"<<endl;
    
    if(I2N>0.0)
    cout<<"CFD rheology core: mu(I) ill posed above I = "<<I2N<<endl;
    }
    
    if(p->W126==2)
    cout<<"CFD rheology core: W 126 2 (advected pore pressure) is not available in CFD, no pore pressure"<<endl;
    }
}

// gravity normal to the level set, into phase 1: (g . n), n = grad(phi)/|grad(phi)|
double rheology_cfd_core::gravity_normal(lexer *p, fdm *a)
{
    double nx = (a->phi(i+1,j,k) - a->phi(i-1,j,k))/(p->DXP[IM1] + p->DXP[IP]);
    double ny = p->j_dir==1 ? (a->phi(i,j+1,k) - a->phi(i,j-1,k))/(p->DYP[JM1] + p->DYP[JP]) : 0.0;
    double nz = (a->phi(i,j,k+1) - a->phi(i,j,k-1))/(p->DZP[KM1] + p->DZP[KP]);
    const double nn = sqrt(nx*nx + ny*ny + nz*nz);
    
    // grad(phi) points into phase 1 (away from the surface): gravity along it
    double gn = sqrt(gx*gx + gy*gy + gz*gz);
    
    if(nn>1.0e-10)
    gn = MAX(0.0,(gx*nx + gy*ny + gz*nz)/nn);
    
    return gn;
}

double rheology_cfd_core::pressure(lexer *p, fdm *a)
{
    const double phival = a->phi(i,j,k);
    
    if(p->W131==2 || p->W131==3)
    {
    const double ps = a->press(i,j,k) - p->pressgage;
    
    if(p->W131==2)
    return ps;
    
    // W 131 3: the solved pressure, at least the weight of the smoothed interface band above the cell. In the band
    // (|phi| < e) the material is a mixture with air whose solved pressure is about that of the air, and the
    // strength mu p' of the band was about zero: the band of an arrested deposit kept flowing along the inclined
    // surface (driven by the horizontal density gradient of the smoothed interface) and the deposit slumped. The
    // floor is the hydrostatic pressure of the mixture below a flat surface, rho_1 (g.n) int_{-e}^{phi} H(s) ds, at
    // most rho_1 (g.n) e (reached at phi = e); below the band the solved pressure is larger
    const double e = interface_width(p,a->phi,i,j,k);
    
    if(phival<=-e)
    return ps;
    
    const double s = MIN(phival,e);
    const double hint = 0.5*(s + s*s/(2.0*e) - e*cos(PI*s/e)/(PI*PI) + 0.5*e - e/(PI*PI));
    
    return MAX(ps, rho1*gravity_normal(p,a)*hint);
    }
    
    // lithostatic: rho_1 (g . n) phi, n the unit normal of the level set into phase 1
    if(phival<=0.0)
    return 0.0;
    
    return rho1*gravity_normal(p,a)*phival;
}

// the cell lies in a direct-forcing solid, the topography or a 6DOF body
bool rheology_cfd_core::in_solid(lexer *p, fdm *a, int ii, int jj, int kk)
{
    return (p->solidread==1 && a->solid(ii,jj,kk)<0.0)
        || (p->toporead==1 && a->topo(ii,jj,kk)<0.0)
        || (p->X10>0 && a->fb(ii,jj,kk)<0.0);
}

double rheology_cfd_core::local(lexer *p, fdm *a, field &u, field &v, field &w)
{
    const double gamma = strainterm(p,u,v,w);
    
    double peff = pressure(p,a);
    
    if(p->W126==1)
    peff *= 1.0 - p->W127_lambda;
    
    return rheo::viscosity(q,gamma,peff)/rho1;
}

double rheology_cfd_core::viscosity(lexer *p, fdm *a, ghostcell *pgc, field &u, field &v, field &w)
{
    if(!in_solid(p,a,i,j,k))
    return local(p,a,u,v,w);
    
    // inside a solid or a body: its velocities are forced, its pressure is not the pressure of the material (the
    // projection solves through the solid). The law with that pressure gave the solid cells next to a granular
    // material a viscosity that grew with the pressure there, the pressure with it, and the run stopped (45 deg flume,
    // Validation/rheology 22). The solid cells take the viscosity of a neighbouring cell of the material (the
    // largest), so that the wall stress of the face is that of the material; cells without such a neighbour mu_min.
    const int i0=i, j0=j, k0=k;
    const int di[6] = {1,-1,0,0,0,0};
    const int dj[6] = {0,0,1,-1,0,0};
    const int dk[6] = {0,0,0,0,1,-1};
    
    double nu = q.mu_min/rho1;
    
    for(int n=0; n<6; ++n)
    {
        if(p->j_dir==0 && dj[n]!=0)
        continue;
        
        const int ii=i0+di[n], jj=j0+dj[n], kk=k0+dk[n];
        
        if(p->flag4[(ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin]<0 || in_solid(p,a,ii,jj,kk))
        continue;
        
        i=ii; j=jj; k=kk;
        nu = MAX(nu,local(p,a,u,v,w));
        i=i0; j=j0; k=k0;
    }
    
    return nu;
}

double rheology_cfd_core::mu_eff(lexer *p, fdm *a, field &u, field &v, field &w)
{
    const double gamma = strainterm(p,u,v,w);
    
    rheology_cfd_3phase *s3 = rheology_cfd_3phase::self;
    
    if(s3==nullptr)
    {
        double peff = pressure(p,a);
        
        if(p->W126==1)
        peff *= 1.0 - p->W127_lambda;
        
        return rheo::viscosity(q,gamma,peff);
    }
    
    return rheo::viscosity(q,gamma,s3->peff(p,a),s3->tanpsi());
}
