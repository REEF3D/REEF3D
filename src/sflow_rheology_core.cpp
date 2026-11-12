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

#include"sflow_rheology_core.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<cmath>

sflow_rheology_core::sflow_rheology_core(lexer *p) : lambda(p), lnew(p)
{
    using namespace rheo;
    
    q.law = p->W120;
    q.tandelta = tan(deg2rad(p->W121));
    q.xi = p->W122;
    q.tz1 = tan(deg2rad(p->W123_z1));
    q.tz2 = tan(deg2rad(p->W123_z2));
    q.tz3 = tan(deg2rad(p->W123_z3));
    q.beta = p->W123_beta;
    q.L = p->W123_L;
    q.mus = tan(deg2rad(p->W124_ds));
    q.mu2 = tan(deg2rad(p->W124_d2));
    q.I0 = p->W124_I0;
    q.d = p->W124_d;
    q.phis = p->W124_phis;
    q.tauy = p->W125_tauy;
    q.K = p->W125_K;
    q.n = p->W125_n;
    
    if(p->mpirank==0)
    {
    const char *name[7] = {"off","Coulomb","Voellmy","Pouliquen-Forterre","mu(I)","Herschel-Bulkley","Herschel-Bulkley-Coulomb"};
    cout<<"SFLOW rheology core: basal law "<<name[q.law>=0 && q.law<=6 ? q.law : 0];
    if(p->W126>0)
    cout<<", pore pressure "<<(p->W126==1 ? "constant" : "advected/dissipating");
    if(p->W20!=0.0 || p->W21!=0.0)
    cout<<", along-slope gravity "<<p->W20<<" "<<p->W21;
    cout<<endl;
    }
}

void sflow_rheology_core::ini(lexer *p, fdm2D *b, ghostcell *pgc)
{
    SLICELOOP4
    lambda(i,j) = p->W126>0 ? p->W127_lambda : 0.0;
    
    pgc->gcsl_start4(p,lambda,50);
}

void sflow_rheology_core::u_source(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    // along-slope gravity (inclined reference plane)
    if(p->W20!=0.0)
    SLICELOOP4
    WETDRY
    b->F(i,j) += p->W20*WL(i,j);
}

void sflow_rheology_core::v_source(lexer *p, fdm2D *b, slice &U, slice &V, slice &WL)
{
    if(p->W21!=0.0)
    SLICELOOP4
    WETDRY
    b->G(i,j) += p->W21*WL(i,j);
}

double sflow_rheology_core::gn(lexer *p, fdm2D *b)
{
    // gravity normal to the bed; W 128: cos^2(theta) of the bed slope (vertical depth -> normal stress)
    double g = fabs(p->W22);
    
    if(p->W128==1)
    {
    const double bx = (b->bed(i+1,j) - b->bed(i-1,j))/(p->DXP[IM1] + p->DXP[IP]);
    const double by = p->j_dir==1 ? (b->bed(i,j+1) - b->bed(i,j-1))/(p->DYP[JM1] + p->DYP[JP]) : 0.0;
    g /= 1.0 + bx*bx + by*by;
    }
    
    return g;
}

void sflow_rheology_core::implicit(lexer *p, fdm2D *b, ghostcell *pgc, slice &UH, slice &VH, slice &WL, double alpha, bool finalize)
{
    const double dt = alpha*p->dt;
    rheo::basal_state s;
    s.rho = p->W1;
    s.g = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    
    SLICELOOP4
    {
        const double h = WL(i,j);
        
        if(p->wet[IJ]==0 || h<=p->A244)
        continue;
        
        const double uh = UH(i,j);
        const double vh = p->j_dir==1 ? VH(i,j) : 0.0;
        const double Ustar = sqrt(uh*uh + vh*vh)/h;
        
        if(Ustar<=0.0 && q.law!=rheo::BASAL_OFF)
        continue;
        
        s.h = h;
        s.gn = gn(p,b);
        s.sigma_n = s.rho*s.gn*h;
        s.lambda = 0.0;
        
        if(p->W126==1)
        s.lambda = lambda(i,j);
        
        // advected ratio: at the middle of the time step (lambda is updated once per step)
        if(p->W126==2)
        s.lambda = rheo::pore_relax(lambda(i,j),p->W127_lambdah,p->W127_D,h,0.5*p->dt);
        
        const double sp = rheo::implicit_speed(q,s,Ustar,dt);
        const double f = sp/Ustar;
        
        UH(i,j) *= f;
        VH(i,j) *= f;
    }
    
    if(finalize && p->W126==2)
    pore_update(p,b,pgc,UH,VH,WL);
}

void sflow_rheology_core::pore_update(lexer *p, fdm2D *b, ghostcell *pgc, slice &UH, slice &VH, slice &WL)
{
    // first-order upwind advection of the ratio with the depth-averaged velocity, then relaxation
    pgc->gcsl_start4(p,lambda,50);
    
    SLICELOOP4
    {
        lnew(i,j) = lambda(i,j);
        const double h = WL(i,j);
        
        if(p->wet[IJ]==0 || h<=p->A244)
        continue;
        
        const double u = UH(i,j)/h;
        const double v = p->j_dir==1 ? VH(i,j)/h : 0.0;
        double adv = 0.0;
        
        if(u>0.0)
        adv += u*(lambda(i,j) - lambda(i-1,j))/p->DXP[IM1];
        else
        adv += u*(lambda(i+1,j) - lambda(i,j))/p->DXP[IP];
        
        if(p->j_dir==1)
        {
        if(v>0.0)
        adv += v*(lambda(i,j) - lambda(i,j-1))/p->DYP[JM1];
        else
        adv += v*(lambda(i,j+1) - lambda(i,j))/p->DYP[JP];
        }
        
        lnew(i,j) = rheo::pore_relax(lambda(i,j) - p->dt*adv,p->W127_lambdah,p->W127_D,h,p->dt);
    }
    
    SLICELOOP4
    lambda(i,j) = lnew(i,j);
    
    pgc->gcsl_start4(p,lambda,50);
}
