/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Hans Bihs, Edgar Chavez
--------------------------------------------------------------------*/

#include"fnpf_runup.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include <math.h>

void fnpf_runup::fnpf_runup_calc(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{	
    double ru_eta, ru_L, ru_k, ru_m, ru_u, ru_g, ru_H, ru_Hm0, ru_Tp, ru_s0p, ru_mLA, ru_mLB, ru_mLC, ru_h, ru_Ur, ru_gammaD; //Free Surface, Wave length, Wave number, VSHT adjusment coefficient, particle velocity, gravity acceleration, Wave height, peak deep water steepness, m Level A, m Level B, m Level C, water depth, Ursell number, Diameter-run-up increasing factor
    
    k==0; //**Initial sigma grid**
    //Calculation of wvae Height and Period
    //ru_Hm0 = (c->eta(i,j)-p->F60) * 2; //Aproximation considering intial SWL and equal through and crest sizes
    
    //ru_Tp =
    
    
    // -> Run-Up Formulae
    ru_eta = c->eta(i,j); //free surface elevation, **eta is given in ALE coordinates not Eulerian?**
    ru_L = p->wL; // Wave length <-----------------Re-Calculate
    ru_k = 2*PI/(ru_L); //wave number
    
    R1 = ru_eta * sqrt(1 + pow(2 * ru_k * rc,2)); //McCarmy and Fuchs (1954)
    
    //------------------------------------------------------------------------------------------------------------------------
    ru_m = 1; //VSHT adjustment coefficient for long waves + LWT or regular waves + 2nd order ST 
    ru_u = acceleration(p,c,pgc) * p->dt; //Particle acceleration at the free surface in lagrangian coordinates
    ru_g = p->W22; //gravity acceleration
    
    R2 = ru_eta + ru_m * pow(ru_u,2)/(2 * ru_g); //Hallermeier (1976) - De Vos et al (2007) only when 2nd order stokes used. 
    
    //------------------------------------------------------------------------------------------------------------------------
    //ru_m = 6.83; //VSHT adjustment coefficient regular and irregular waves + LWT
    
    //R3 = ru_eta + ru_m * pow(ru_u,2)/(2 * ru_g); //Niedzwecki and Duggal (1992) **need H or SWL to calculate velocity at this point**
    
    //------------------------------------------------------------------------------------------------------------------------
    //ru_m = 6.52; //VSHT adjustment coefficient regular and irregular waves + LWT altering linear fit
    
    //R4 = ru_eta + ru_m * pow(ru_u,2)/(2 * ru_g); //Niedzwecki and Huston (1992) **need H or SWL to calculate velocity at this point**
    
    //------------------------------------------------------------------------------------------------------------------------
    ru_Hm0 = p->wHs; //Significant wave height <-----------------Re-Calculate with H per wave
    ru_Tp = p->wTp; //Peak wave period <-----------------Re-Calculate with T per wave
    ru_s0p = 2 * PI * ru_Hm0 / (ru_g * pow(ru_Tp,2)); //Peak deep water steepness
    
    if (ru_s0p <= 0.035)
    {
        ru_mLA = -66.667 * ru_s0p + 5.33; //Green water Run-up (Thick layer)
        ru_mLB = -93.333 * ru_s0p + 7.47; //Thin layer of water an air mixture
        ru_mLC = -200 * ru_s0p + 16; //Maximum spray
    }
    else{
        ru_mLA = 3;
        ru_mLB = 4.2;
        ru_mLC = 9;
    }
    
    ru_m = ru_mLB; //Consider Run-up height up to the thin layer
    
    R5 = ru_eta + ru_m * pow(ru_u,2)/(2 * ru_g); //Ramirez et al (2013)
    
    //------------------------------------------------------------------------------------------------------------------------
    ru_H = p->wd; //Wave height <-----------------Re-Calculate
    ru_h = p->wd; //Water depth <-----------------Consider new parameter
    ru_Ur = ru_H * pow(ru_L,2) / pow(ru_h,3); //Urell number
    ru_gammaD = 0.004 * log(251.8 * ru_H / (2 * rc) + 1); //Diameter-run-up increasing factor
    
    R6 = 7.39 * ru_gammaD * log(0.27 * ru_Ur +1); //Peng et al (2012)
    
    //Temporary calculations
    R1 = R1; //Fee surface
    R2 = R2;
    R3 = ru_eta + ru_mLA * pow(ru_u,2)/(2 * ru_g);
    R4 = ru_eta + ru_mLA * pow(ru_u,2)/(2 * ru_g);
    R5 = ru_eta + ru_mLA * pow(ru_u,2)/(2 * ru_g);
    R6 = R6;
    
    // Storing current time step information for next time step gradient calculation
    un = c->U[FIJK];
    vn = c->V[FIJK];
    etan = c->eta(i,j);
    
    cout<<" R1: "<<R1<<endl;
}


double fnpf_runup::acceleration(lexer *p, fdm_fnpf *c, ghostcell *pgc) //Particle total acceleration on the surface
{
    double ax1,ay1,ax2,ay2,ax3,ay3,dudsig_,dvdsig_,ax,ay;
    
    dudsig_= dudsig(p,c,pgc);
    dvdsig_= dvdsig(p,c,pgc);
 
    // Equation (9) o Pakozdi et al (2021)
    // Term 1
    ax1= (c->U[FIJK] - un)/(p->dt);
    ay1= (c->V[FIJK] - vn)/ (p->dt);
    // Term 2
    ax2 = c->U[FIJK] * (dudxi(p,c,pgc) + (dudsig_ * p->sigx[FIJK]));
    ay2 = c->V[FIJK] * (dvdxi(p,c,pgc) + (dvdsig_ * p->sigy[FIJK]));
    // Term 3
    ax3 = (c->W[FIJK] - (p->sig[FIJK] * dndt(p,c,pgc))) * dudsig_ * p->sigz[IJ];
    ay3 = (c->W[FIJK] - (p->sig[FIJK] * dndt(p,c,pgc))) * dvdsig_ * p->sigz[IJ];
    // Sum up acceleration
    ax = ax1 + ax2 + ax3;
    ay = ay1 + ay2 + ay3;
        
    return sqrt(pow(ax,2)+pow(ay,2));
}

double fnpf_runup::dndt(lexer *p, fdm_fnpf *c, ghostcell *pgc) //Free surface elevation derivated on time
{
    return (c->eta(i,j) - etan)/ p->dt;
}


double fnpf_runup::dudsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) //Velocity in X direction derivated on sigma grid
{
	return (c->U[FIJKp1] - c->U[FIJKm1])/(p->DZN[KP1] + p->DZN[KM1]);     
}


double fnpf_runup::dvdsig(lexer *p, fdm_fnpf *c, ghostcell *pgc) //Velocity in Y direction derivated on sigma grid
{
	return (c->V[FIJKp1] - c->V[FIJKm1])/(p->DZN[KP1] + p->DZN[KM1]);       
}


double fnpf_runup::dudxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) //Velocity in Y direction derivated on sigma grid
{
    return (c->U[FIp1JK] - c->U[FIm1JK])/(p->DXN[IP1] + p->DXN[IM1]); 
}


double fnpf_runup::dvdxi(lexer *p, fdm_fnpf *c, ghostcell *pgc) 	// getting dvdxi
{
    return (c->V[FIJp1K] - c->V[FIJm1K])/(p->DYN[JP1] + p->DYN[JM1]); 
}






