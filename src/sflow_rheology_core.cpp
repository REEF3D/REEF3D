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

double sflow_rheology_core::gn_at(lexer *p, fdm2D *b, int ii, int jj)
{
    const int si = i, sj = j;
    i = ii;
    j = jj;
    const double r = gn(p,b);
    i = si;
    j = sj;
    return r;
}

double sflow_rheology_core::lambda_mid(lexer *p, double lam, double h)
{
    // pore pressure ratio of the friction: constant (W 126 1), or the advected ratio at the middle of the time step
    // (W 126 2; lambda is updated once per step)
    if(p->W126==1)
    return lam;
    
    if(p->W126==2)
    return rheo::pore_relax(lam,p->W127_lambdah,p->W127_D,h,0.5*p->dt);
    
    return 0.0;
}

double sflow_rheology_core::si_factor(lexer *p, fdm2D *b, int ii, int jj, int dir, double h, double Ustar, double dt)
{
    // semi-implicit step: the speed s of the backward-Euler friction step s + dt tau_b(s)/(rho h) = Ustar, as the
    // factor s/Ustar = 1/(1 + dt tau_b(s)/(rho h s)) on the momentum of the cell or the face (with the implicit
    // gravity term); 0 when the driving force of the step is below the static threshold. At a face the normal gravity
    // and the pore pressure ratio are the means of the two cells
    if(h<=p->A244)
    return 1.0;
    
    if(Ustar<=0.0)
    return 0.0;
    
    const int i2 = dir==1 ? ii+1 : ii;
    const int j2 = dir==2 ? jj+1 : jj;
    
    rheo::basal_state s;
    s.rho = p->W1;
    s.g = sqrt(p->W20*p->W20 + p->W21*p->W21 + p->W22*p->W22);
    s.h = h;
    s.gn = dir==0 ? gn_at(p,b,ii,jj) : 0.5*(gn_at(p,b,ii,jj) + gn_at(p,b,i2,j2));
    s.sigma_n = s.rho*s.gn*h;
    s.lambda = lambda_mid(p,dir==0 ? lambda(ii,jj) : 0.5*(lambda(ii,jj) + lambda(i2,j2)),h);
    
    return rheo::implicit_speed(q,s,Ustar,dt)/Ustar;
}

void sflow_rheology_core::si_end(lexer *p, fdm2D *b, ghostcell *pgc, slice &UH, slice &VH, slice &WL)
{
    if(p->W126==2)
    pore_update(p,b,pgc,UH,VH,WL);
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
        s.lambda = lambda_mid(p,lambda(i,j),h);
        
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
