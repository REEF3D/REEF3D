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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"bcmom.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<math.h>

bcmom::bcmom(lexer* p):surftens(p),roughness(p),kappa(0.4)
{
    if(p->F50==1)
    gcval_phi=51;

    if(p->F50==2)
    gcval_phi=52;

    if(p->F50==3)
    gcval_phi=53;

    if(p->F50==4)
    gcval_phi=54;


    bckin=0;
    if(p->T10>0 || p->T10<20)
    bckin=1;
    
    // Initialize wall function constants
    E = 9.8;                 // Wall constant for smooth walls
    yPlusViscous = 11.0;     // Transition point between viscous sublayer and log region
    ksPlusSmooth = 2.25;     // Hydraulically smooth limit 
    ksPlusRough = 90.0;      // Fully rough limit 
    
}

bcmom::~bcmom()
{
}

void bcmom::bcmom_start(fdm* a, lexer* p,ghostcell *pgc, turbulence *pturb,field& b,int gcval)
{
    int q;

    if(gcval==10 && p->B10!=0)
    {
        QGC1LOOP
        if((p->gcb1[q][4]==5 || p->gcb1[q][4]==21 || p->gcb1[q][4]==22 || p->gcb1[q][4]==41 || p->gcb1[q][4]==42 || p->gcb1[q][4]==43) && p->gcb1[q][3]!=1 && p->gcb1[q][3]!=4)
        wall_law_u(a,p,pturb,b,p->gcb1[q][0], p->gcb1[q][1], p->gcb1[q][2], p->gcb1[q][3], p->gcb1[q][4], p->gcd1[q]);
        
        QGCDF1LOOP
        wall_law_u(a,p,pturb,b,p->gcdf1[q][0], p->gcdf1[q][1], p->gcdf1[q][2], p->gcdf1[q][3], p->gcdf1[q][4],  0.5*p->DXM);
    }

    if(gcval==11 && p->B10!=0 && p->j_dir==1)
    {
        QGC2LOOP
        if((p->gcb2[q][4]==5 || p->gcb2[q][4]==21 || p->gcb2[q][4]==22 || p->gcb2[q][4]==41 || p->gcb2[q][4]==42 || p->gcb2[q][4]==43) && p->gcb2[q][3]!=2 && p->gcb2[q][3]!=3)
        wall_law_v(a,p,pturb,b,p->gcb2[q][0], p->gcb2[q][1], p->gcb2[q][2], p->gcb2[q][3], p->gcb2[q][4], p->gcd2[q]);
        
        QGCDF2LOOP
        wall_law_v(a,p,pturb,b,p->gcdf2[q][0], p->gcdf2[q][1], p->gcdf2[q][2], p->gcdf2[q][3], p->gcdf2[q][4],  0.5*p->DXM);
    }

    if(gcval==12 && p->B10!=0)
    {
        QGC3LOOP
        if((p->gcb3[q][4]==5 || p->gcb3[q][4]==21 || p->gcb3[q][4]==22 || p->gcb3[q][4]==41 || p->gcb3[q][4]==42 || p->gcb3[q][4]==43) && p->gcb3[q][3]!=5 && p->gcb3[q][3]!=6)
        wall_law_w(a,p,pturb,b,p->gcb3[q][0], p->gcb3[q][1], p->gcb3[q][2], p->gcb3[q][3], p->gcb3[q][4], p->gcd3[q]);
        
        QGCDF3LOOP
        wall_law_w(a,p,pturb,b,p->gcdf3[q][0], p->gcdf3[q][1], p->gcdf3[q][2], p->gcdf3[q][3], p->gcdf3[q][4],  0.5*p->DXM);

    }
    surface_tension(a,p,a->phi,gcval);
}

void bcmom::wall_law_u(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Adjust distance based on cell orientation
    if(cs==2 || cs==3)
    dist=p->DYN[JP];
    
    if(cs==5 || cs==6)
    dist=p->DZN[KP];
    
    // Get roughness height
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Ensure minimum distance for numerical stability
    if(dist < 1.0e-10)
    dist = 1.0e-10;
    
    // Get local properties
    double nu = a->visc(i,j,k);      // Kinematic viscosity
    double rho = a->ro(i,j,k);       // Density
    double velocity = a->u(i,j,k);   // Local velocity
    double velMag = fabs(velocity);  // Velocity magnitude
    
    // Avoid division by zero
    if(velMag < 1.0e-10)
    {
        return; // No wall stress for essentially zero velocity
    }
    
    double wallStress; // Will store the calculated wall stress
    
    // Initial guess for friction velocity using log law
    double uTau_guess = velMag * kappa / log(E * dist/max(ks, 1.0e-10));
    
    // Calculate y+
    double yPlus = dist * uTau_guess / max(nu, 1.0e-10);
    
    // Calculate ks+
    double ksPlus = ks * uTau_guess / max(nu, 1.0e-10);
    
    // Calculate wall stress based on appropriate wall model
    if(yPlus < yPlusViscous)
    {
        // Viscous sublayer - linear profile
        // tau_w = mu * dU/dy = mu * U/y
        wallStress = rho * nu * velMag / dist;
    }
    else
    {
        // Log law region with roughness effects
        double uTau; // Friction velocity
        
        if(ksPlus < ksPlusSmooth)
        {
            // Hydraulically smooth regime
            uTau = velMag * kappa / log(E * yPlus);
        }
        else if(ksPlus < ksPlusRough)
        {
            // Transitionally rough regime
            // Use Cebeci-Chang roughness model
            double delta_B = (1.0/kappa) * log(ksPlus) * sin(0.4258 * (log(ksPlus) - 0.811));
            uTau = velMag * kappa / (log(E * yPlus) - delta_B);
        }
        else
        {
            // Fully rough regime
            uTau = velMag * kappa / log(30.0 * dist/ks);
        }
        
        // Wall stress from friction velocity
        wallStress = rho * uTau * uTau;
    }
    
    // Apply wall stress with correct sign
    a->F(i,j,k) -= (wallStress * (velocity/velMag));
}

void bcmom::wall_law_v(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Adjust distance based on cell orientation
    if(cs==1 || cs==4)
    dist=p->DXN[IP];
    
    if(cs==5 || cs==6)
    dist=p->DZN[KP];
    
    // Get roughness height
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Ensure minimum distance for numerical stability
    if(dist < 1.0e-10)
    dist = 1.0e-10;
    
    // Get local properties
    double nu = a->visc(i,j,k);      // Kinematic viscosity
    double rho = a->ro(i,j,k);       // Density
    double velocity = a->v(i,j,k);   // Local velocity
    double velMag = fabs(velocity);  // Velocity magnitude
    
    // Avoid division by zero
    if(velMag < 1.0e-10)
    {
        return; // No wall stress for essentially zero velocity
    }
    
    double wallStress; // Will store the calculated wall stress
    
    // Initial guess for friction velocity using log law
    double uTau_guess = velMag * kappa / log(E * dist/max(ks, 1.0e-10));
    
    // Calculate y+
    double yPlus = dist * uTau_guess / max(nu, 1.0e-10);
    
    // Calculate ks+
    double ksPlus = ks * uTau_guess / max(nu, 1.0e-10);
    
    // Calculate wall stress based on appropriate wall model
    if(yPlus < yPlusViscous)
    {
        // Viscous sublayer - linear profile
        // tau_w = mu * dU/dy = mu * U/y
        wallStress = rho * nu * velMag / dist;
    }
    else
    {
        // Log law region with roughness effects
        double uTau; // Friction velocity
        
        if(ksPlus < ksPlusSmooth)
        {
            // Hydraulically smooth regime
            uTau = velMag * kappa / log(E * yPlus);
        }
        else if(ksPlus < ksPlusRough)
        {
            // Transitionally rough regime
            // Use Cebeci-Chang roughness model
            double delta_B = (1.0/kappa) * log(ksPlus) * sin(0.4258 * (log(ksPlus) - 0.811));
            uTau = velMag * kappa / (log(E * yPlus) - delta_B);
        }
        else
        {
            // Fully rough regime
            uTau = velMag * kappa / log(30.0 * dist/ks);
        }
        
        // Wall stress from friction velocity
        wallStress = rho * uTau * uTau;
    }
    
    // Apply wall stress with correct sign
    a->G(i,j,k) -= (wallStress * (velocity/velMag));
}

void bcmom::wall_law_w(fdm* a,lexer* p, turbulence *pturb,field& b,int ii,int jj,int kk,int cs,int bc,double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Adjust distance based on cell orientation
    if(cs==1 || cs==4)
    dist=p->DXN[IP];
    
    if(cs==2 || cs==3)
    dist=p->DYN[JP];
    
    // Get roughness height
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Ensure minimum distance for numerical stability
    if(dist < 1.0e-10)
    dist = 1.0e-10;
    
    // Get local properties
    double nu = a->visc(i,j,k);      // Kinematic viscosity
    double rho = a->ro(i,j,k);       // Density
    double velocity = a->w(i,j,k);   // Local velocity
    double velMag = fabs(velocity);  // Velocity magnitude
    
    // Avoid division by zero
    if(velMag < 1.0e-10)
    {
        return; // No wall stress for essentially zero velocity
    }
    
    double wallStress; // Will store the calculated wall stress
    
    // Initial guess for friction velocity using log law
    double uTau_guess = velMag * kappa / log(E * dist/max(ks, 1.0e-10));
    
    // Calculate y+
    double yPlus = dist * uTau_guess / max(nu, 1.0e-10);
    
    // Calculate ks+
    double ksPlus = ks * uTau_guess / max(nu, 1.0e-10);
    
    // Calculate wall stress based on appropriate wall model
    if(yPlus < yPlusViscous)
    {
        // Viscous sublayer - linear profile
        // tau_w = mu * dU/dy = mu * U/y
        wallStress = rho * nu * velMag / dist;
    }
    else
    {
        // Log law region with roughness effects
        double uTau; // Friction velocity
        
        if(ksPlus < ksPlusSmooth)
        {
            // Hydraulically smooth regime
            uTau = velMag * kappa / log(E * yPlus);
        }
        else if(ksPlus < ksPlusRough)
        {
            // Transitionally rough regime
            // Use Cebeci-Chang roughness model
            double delta_B = (1.0/kappa) * log(ksPlus) * sin(0.4258 * (log(ksPlus) - 0.811));
            uTau = velMag * kappa / (log(E * yPlus) - delta_B);
        }
        else
        {
            // Fully rough regime
            uTau = velMag * kappa / log(30.0 * dist/ks);
        }
        
        // Wall stress from friction velocity
        wallStress = rho * uTau * uTau;
    }
    
    // Apply wall stress with correct sign
    a->H(i,j,k) -= (wallStress * (velocity/velMag));
}




