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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cmath>

double rheology_f::yield_stress(lexer* p, fdm* a)
{
    double tau0;
    pressureval = MAX(0.0,pressureval);
    // tanphi >0 is enforced
    double relative_density = MAX(0.0,a->ro(i,j,k)-density_interstitial_fluid)/a->ro(i,j,k);
    gamma = strainterm(p,a);

    switch(p->W101)
    {
    default:
    case 0: // HB prescribed yield stress
        tau0=p->W96;
        break;

    case 1:  // Mohr-Coulomb
        tau0 = MAX(0.0,tanphi*pressureval + p->W102_c);
        break;
        
    case 2:  // Yield stress from viscoplastic fluid (bingham) Druker-Prager - tanphi*pressureval valid in 2D - W102_c is cohesion
        tau0 = (tanphi*pressureval + p->W102_c)*(1.0-exp(-p->W103*gamma));
        break;
        
    case 3:  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = (tanphi*pressureval*relative_density + p->W102_c)*(1.0-exp(-p->W103*gamma));    // rho_water = 1000.0, new input?
        break;

    case 4:  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*pressureval*exp(-p->W104*gamma)*relative_density + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_p is new input W 104 
        break;

    case 5:  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,pressureval*relative_density-p->W104*gamma) + p->W102_c)*(1.0-exp(-p->W103*gamma));    // m_u also use new input W 104
        break;
    }

    if(p->count==0)
    tau0 = p->W96;
    
    return tau0;
}

void rheology_f::yieldStressGradient(lexer* p, fdm* a, int ii, int jj, int kk)
{
    switch(p->W101)
    {
    default:
    case 0:
        tau01=p->W96;
        tau02=p->W96;
        break;

    case 1:  // Mohr-Coulomb
        tau01 = (tanphi*pressureval1 + p->W102_c);
        tau02 = (tanphi*pressureval2 + p->W102_c);
        break;
        
    case 2:  // Yield stress from viscoplastic fluid (bingham) Druker-Prager - tanphi*pressureval valid in 2D - W102_c is cohesion
        tau01 = (tanphi*pressureval1 + p->W102_c)*(1.0-exp(-p->W103*gamma));
        tau02 = (tanphi*pressureval2 + p->W102_c)*(1.0-exp(-p->W103*gamma));
        break;
        
    case 3:  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau01 = MAX(0.0,tanphi*pressureval1*MAX(0.0,a->ro(i,j,k)-density_interstitial_fluid)/a->ro(i,j,k) + p->W102_c);
        tau02 = MAX(0.0,tanphi*pressureval2*MAX(0.0,a->ro(i+1,j,k)-density_interstitial_fluid)/a->ro(i+1,j,k) + p->W102_c);
        break;
        
    case 4:  // HB-C shear rate generated excess pore pressure
        tau01 = MAX(0.0,tanphi*pressureval1*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-density_interstitial_fluid)/a->ro(i,j,k) + p->W102_c);
        tau02 = MAX(0.0,tanphi*pressureval2*exp(-p->W104*gamma)*MAX(0.0,a->ro(i+1,j,k)-density_interstitial_fluid)/a->ro(i+1,j,k) + p->W102_c);
        break;
        
    case 5:  // HB-C linear shear rate coupling, max given by pressure
        tau01 = MAX(0.0,tanphi*MAX(0.0,pressureval1*MAX(0.0,a->ro(i,j,k)-density_interstitial_fluid)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c);
        tau02 = MAX(0.0,tanphi*MAX(0.0,pressureval2*MAX(0.0,a->ro(i+1,j,k)-density_interstitial_fluid)/a->ro(i+1,j,k)-p->W104*gamma) + p->W102_c);
        break;
    }   

    if(p->count==0)
    {
        tau01=p->W96;
        tau02=p->W96;
    }
}

void rheology_f::pressurePhi(lexer* p, fdm* a, int ii, int jj, int kk, bool pressureGauge)
{
    phival = 0.5*(a->phi(i,j,k)+a->phi(i+ii,j+jj,k+kk));

    switch(p->W111)
    {
    default:
    case 1:
        pressureval=phival*0.5*(a->ro(i,j,k)+a->ro(i+ii,j+jj,k+kk))*gravity;
        break;
    
    case 2:
        pressureval=0.5*(a->press(i,j,k)+a->press(i+ii,j+jj,k+kk))-p->pressgage*pressureGauge;
        break;
    
    case 3:
        if(phival<p->W112*p->DXM)
            pressureval=0.5*(a->press(i,j,k)+a->press(i+ii,j+jj,k+kk))-p->pressgage*pressureGauge;
        else
            pressureval=phival*0.5*(a->ro(i,j,k)+a->ro(i+ii,j+jj,k+kk))*gravity;
        break;
    
    case 4:
        pressureval=0.5*(0.5*(a->press(i,j,k)+a->press(i+ii,j+jj,k+kk))-p->pressgage*pressureGauge + phival*0.5*(a->ro(i,j,k)+a->ro(i+ii,j+jj,k+kk))*gravity);
        break;
    
    case 5:
        if(p->count>=10)
            pressureval=0.5*(a->press(i,j,k)+a->press(i+ii,j+jj,k+kk))-p->pressgage*pressureGauge;
        else
            pressureval=phival*0.5*(a->ro(i,j,k)+a->ro(i+ii,j+jj,k+kk))*gravity;
        break;
    }
}

void rheology_f::pressurePhiGradient(lexer* p, fdm* a, int ii, int jj, int kk)
{
    phival = 0.5*(a->phi(i,j,k)+a->phi(i+ii,j+jj,k+kk));
       
    switch(p->W111)
    {
    default:
    case 1:
        pressureval1 = phival*a->ro(i,j,k)*gravity;
        pressureval2 = phival*a->ro(i+ii,j+jj,k+kk)*gravity;
        break;
    
    case 2:
        pressureval1 = a->press(i,j,k);
        pressureval2 = a->press(i+ii,j+jj,k+kk);
        break;
    
    case 3:
        if(phival<p->W112*p->DXM)
        {
            pressureval1 = a->press(i,j,k);
            pressureval2 = a->press(i+ii,j+jj,k+kk);
        }
        else
        {
            pressureval1 = phival*a->ro(i,j,k)*gravity;
            pressureval2 = phival*a->ro(i+ii,j+jj,k+kk)*gravity;
        }
        break;
    }
}
