/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"bc_ikepsilon.h"
#include"fdm.h"
#include"lexer.h"

bc_ikepsilon::bc_ikepsilon(lexer* p):roughness(p),kappa(0.4)
{
    // Initialize wall function constants
    E = 9.8;                 // Wall constant for smooth walls
    yPlusLam = 11.0;         // Transition point between viscous sublayer and log region
    ksPlusSmooth = 2.25;     // Hydraulically smooth limit (OpenFOAM value)
    ksPlusRough = 90.0;      // Fully rough limit (OpenFOAM value)
}

bc_ikepsilon::~bc_ikepsilon()
{
}

void bc_ikepsilon::bckeps_start(fdm* a,lexer* p,field& kin,field& eps,int gcval)
{
    int q;

    if(gcval==20)
    {
        QGC4LOOP
        if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21  || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43)
        wall_law_kin(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
        
    n=0;
    LOOP
    {
        if(p->flag4[Im1JK]<0)
        {
        a->rhsvec.V[n] -= a->M.s[n]*kin(i-1,j,k);
        a->M.s[n] = 0.0;
        }
        
        if(p->flag4[Ip1JK]<0)
        {
        a->rhsvec.V[n] -= a->M.n[n]*kin(i+1,j,k);
        a->M.n[n] = 0.0;
        }
        
        if(p->flag4[IJm1K]<0 && p->j_dir==1)
        {
        a->rhsvec.V[n] -= a->M.e[n]*kin(i,j-1,k);
        a->M.e[n] = 0.0;
        }
        
        if(p->flag4[IJp1K]<0 && p->j_dir==1)
        {
        a->rhsvec.V[n] -= a->M.w[n]*kin(i,j+1,k);
        a->M.w[n] = 0.0;
        }
        
        if(p->flag4[IJKm1]<0)
        {
        a->rhsvec.V[n] -= a->M.b[n]*kin(i,j,k-1);
        a->M.b[n] = 0.0;
        }
        
        if(p->flag4[IJKp1]<0)
        {
        a->rhsvec.V[n] -= a->M.t[n]*kin(i,j,k+1);
        a->M.t[n] = 0.0;
        }

    ++n;
    }
    }

    if(gcval==30)
    {
        QGC4LOOP
        if(p->gcb4[q][4]==5 || p->gcb4[q][4]==21 || p->gcb4[q][4]==22 || p->gcb4[q][4]==41 || p->gcb4[q][4]==42 || p->gcb4[q][4]==43  || (p->gcb4[q][4]==3 && p->gcb4[q][4]==6))
        wall_law_eps(a,p,kin,eps,p->gcb4[q][0], p->gcb4[q][1], p->gcb4[q][2], p->gcb4[q][3], p->gcb4[q][4], p->gcb4[q][5],  p->gcd4[q]);
        
    n=0;
    LOOP
    {
        if(p->flag4[Im1JK]<0)
        {
        a->rhsvec.V[n] -= a->M.s[n]*eps(i-1,j,k);
        a->M.s[n] = 0.0;
        }
        
        if(p->flag4[Ip1JK]<0)
        {
        a->rhsvec.V[n] -= a->M.n[n]*eps(i+1,j,k);
        a->M.n[n] = 0.0;
        }
        
        if(p->flag4[IJm1K]<0 && p->j_dir==1)
        {
        a->rhsvec.V[n] -= a->M.e[n]*eps(i,j-1,k);
        a->M.e[n] = 0.0;
        }
        
        if(p->flag4[IJp1K]<0 && p->j_dir==1)
        {
        a->rhsvec.V[n] -= a->M.w[n]*eps(i,j+1,k);
        a->M.w[n] = 0.0;
        }
        
        if(p->flag4[IJKm1]<0)
        {
        a->rhsvec.V[n] -= a->M.b[n]*eps(i,j,k-1);
        a->M.b[n] = 0.0;
        }
        
        if(p->flag4[IJKp1]<0)
        {
        a->rhsvec.V[n] -= a->M.t[n]*eps(i,j,k+1);
        a->M.t[n] = 0.0;
        }

    ++n;
    }
    }


}
// ****************************
// WALL KIN
// ****************************
void bc_ikepsilon::wall_law_kin(fdm* a,lexer* p,field& kin,field& eps,int ii,int jj,int kk,int cs,int bc, int id, double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Calculate distance to wall
    if(cs==1 || cs==4)
    dist = 0.5*p->DXN[IP];
    
    if(cs==2 || cs==3)
    dist = 0.5*p->DYN[JP];
    
    if(cs==5 || cs==6)
    dist = 0.5*p->DZN[KP];
    
    // Get roughness height
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Ensure minimum distance for numerical stability
    if(dist < 1.0e-10)
    dist = 1.0e-10;

    // Calculate velocity at the cell center
    double uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
    double vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
    double wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
    
    // Velocity magnitude
    u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
    
    // Avoid division by zero
    if(u_abs < 1.0e-10)
    {
        return; // No wall stress for essentially zero velocity
    }
    
    // Get local viscosity
    double nu = a->visc(i,j,k);
    
    // Improved iterative calculation of friction velocity
    // Initial guess based on log law
    double uTau = u_abs * kappa / log(E * dist/MAX(ks, 1.0e-10));
    double uTau_old = 0.0;
    int maxIter = 10;      // Maximum number of iterations
    double tolerance = 1e-4; // Convergence tolerance
    
    // Iterative refinement of uTau
    for(int iter = 0; iter < maxIter; ++iter)
    {
        // Calculate y+ and ks+ with current uTau estimate
        yPlus = dist * uTau / MAX(nu, 1.0e-10);
        ks_plus = ks * uTau / MAX(nu, 1.0e-10);
        
        // Store previous uTau for convergence check
        uTau_old = uTau;
        
        // Calculate new uTau based on appropriate wall model
        if(yPlus < yPlusLam)
        {
            // Viscous sublayer - linear profile
            uTau = sqrt(nu * u_abs / dist);
        }
        else
        {
            // Log law region with roughness effects
            if(ks_plus < ksPlusSmooth)
            {
                // Hydraulically smooth regime
                uTau = u_abs * kappa / log(E * yPlus);
            }
            else if(ks_plus < ksPlusRough)
            {
                // Transitionally rough regime
                // Use Cebeci-Chang roughness model
                double delta_B = (1.0/kappa) * log(ks_plus) * sin(0.4258 * (log(ks_plus) - 0.811));
                uTau = u_abs * kappa / (log(E * yPlus) - delta_B);
            }
            else
            {
                // Fully rough regime
                uTau = u_abs * kappa / log(30.0 * dist/ks);
            }
        }
        
        // Check for convergence
        if(fabs(uTau - uTau_old) < tolerance * uTau)
            break;
        
        // Under-relaxation to improve stability
        uTau = 0.7 * uTau + 0.3 * uTau_old;
    }
    
    // Wall shear stress
    tau = uTau * uTau;
    
    // Calculate k wall value based on friction velocity
    double k_wall = tau / sqrt(p->cmu);
    
    // Apply wall function boundary condition for k
    a->M.p[id] += (pow(p->cmu,0.75)*sqrt(MAX(kin(i,j,k),0.0))*uplus)/dist;
    a->rhsvec.V[id] += (k_wall)/dist;
}

void bc_ikepsilon::wall_law_eps(fdm* a,lexer* p,field& kin,field& eps,int ii,int jj,int kk,int cs,int bc, int id, double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Calculate distance to wall
    if(cs==1 || cs==4)
    dist = 0.5*p->DXN[IP];
    
    if(cs==2 || cs==3)
    dist = 0.5*p->DYN[JP];
    
    if(cs==5 || cs==6)
    dist = 0.5*p->DZN[KP];
    
    // Ensure minimum distance for numerical stability
    if(dist < 1.0e-10)
    dist = 1.0e-10;
    
    // Get roughness height
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Calculate velocity at the cell center
    double uvel=0.5*(a->u(i,j,k)+a->u(i-1,j,k));
    double vvel=0.5*(a->v(i,j,k)+a->v(i,j-1,k));
    double wvel=0.5*(a->w(i,j,k)+a->w(i,j,k-1));
    
    // Velocity magnitude
    u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
    
    // Avoid division by zero
    if(u_abs < 1.0e-10)
    {
        return; // No wall stress for essentially zero velocity
    }
    
    // Get local viscosity
    double nu = a->visc(i,j,k);
    
    // Calculate friction velocity using the same method as in wall_law_kin
    double uTau = u_abs * kappa / log(E * dist/MAX(ks, 1.0e-10));
    double uTau_old = 0.0;
    int maxIter = 10;      // Maximum number of iterations
    double tolerance = 1e-4; // Convergence tolerance
    
    // Iterative refinement of uTau (same as in wall_law_kin)
    for(int iter = 0; iter < maxIter; ++iter)
    {
        // Calculate y+ and ks+ with current uTau estimate
        yPlus = dist * uTau / MAX(nu, 1.0e-10);
        ks_plus = ks * uTau / MAX(nu, 1.0e-10);
        
        // Store previous uTau for convergence check
        uTau_old = uTau;
        
        // Calculate new uTau based on appropriate wall model
        if(yPlus < yPlusLam)
        {
            // Viscous sublayer - linear profile
            uTau = sqrt(nu * u_abs / dist);
        }
        else
        {
            // Log law region with roughness effects
            if(ks_plus < ksPlusSmooth)
            {
                // Hydraulically smooth regime
                uTau = u_abs * kappa / log(E * yPlus);
            }
            else if(ks_plus < ksPlusRough)
            {
                // Transitionally rough regime
                // Use Cebeci-Chang roughness model
                double delta_B = (1.0/kappa) * log(ks_plus) * sin(0.4258 * (log(ks_plus) - 0.811));
                uTau = u_abs * kappa / (log(E * yPlus) - delta_B);
            }
            else
            {
                // Fully rough regime
                uTau = u_abs * kappa / log(30.0 * dist/ks);
            }
        }
        
        // Check for convergence
        if(fabs(uTau - uTau_old) < tolerance * uTau)
            break;
        
        // Under-relaxation to improve stability
        uTau = 0.7 * uTau + 0.3 * uTau_old;
    }
    
    // Calculate epsilon at the wall based on friction velocity
    // Using standard wall function for epsilon
    eps_star = pow(uTau, 3.0) / (kappa * dist);
    
    // Set epsilon value at the wall adjacent cell
    eps(i,j,k) = eps_star;
}



