/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"fnpf_print_Hs.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

fnpf_print_Hs::fnpf_print_Hs(lexer *p, fdm_fnpf *c)
{

	
}

fnpf_print_Hs::~fnpf_print_Hs()
{

}

void fnpf_print_Hs::start(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    // DKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAF               
   // DKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAF
  /* 
    if (wfcall!=123456) // Initialise
    {
    wtime  = 0.0;
    T_sum  = 0.0;
    NumDT1 = 0.0;
    SLICELOOP4
    {
    c->ETAsum(i,j)        = 0.0;
	  c->ETAmean(i,j)       = 0.0;
	  c->ETA2sum(i,j)       = 0.0;
	  c->ETAvar(i,j)        = 0.0;
	  c->SigWaveHeight(i,j) = 0.0;
    }
    } // End initialise
   
    wfcall=123456;
    wtime  += p->dt; //DKAF
    if (wtime>stime) // Check we're past transients
    {
    T_sum  += p->dt; //DKAF
    NumDT1 += 1;    //DKAF
    
    //cout << "%%%%%%%%%%%%%%%%%%%%% " << NumDT1 <<endl;
    //cout << "%%%%%%%%%%%%%%%%%%%%% " << T_sum <<endl;
    
      SLICELOOP4
    {
	 // Here we do the wave-averaging NB: c->eta(i,j) is the FS
	 // variance equation with etamean initially unknown
      
    c->ETAsum(i,j)       += c->eta(i,j)*p->dt;
	  c->ETAmean(i,j)       = c->ETAsum(i,j)/T_sum;
	  c->ETA2sum(i,j)      += c->eta(i,j)*c->eta(i,j);
    //cout << "T_sum " << T_sum << " wtim " << wtime <<endl;
    //cin.get();  
    if(T_sum>=T_INTV_mean)
	  { 
	    c->ETAvar(i,j)        = 1.0/double(NumDT1-1)*c->ETA2sum(i,j)-c->ETAmean(i,j)*c->ETAmean(i,j)*double(NumDT1)/double(NumDT1-1);
	    c->SigWaveHeight(i,j) = 4.0*sqrt(c->ETAvar(i,j));
    }
	  
    } //end slice4loop	
 
	  if(T_sum>=T_INTV_mean)
	  {
	    T_sum  -= T_INTV_mean; // Doesn't need looping	
      NumDT1  = 0;
      SLICELOOP4
      {
	    c->ETAsum(i,j)  = 0.0;
      c->ETA2sum(i,j) = 0.0;
      }
    }	
    
   }// End check we're past transients
  // DKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAF
  // DKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAF
*/



/*
 SLICELOOP4
    {
    c->K(i,j) =  - c->Fx(i,j)*c->Ex(i,j) - c->Fy(i,j)*c->Ey(i,j)
    
                 + c->Fz(i,j)*(1.0 + pow(c->Ex(i,j),2.0) + pow(c->Ey(i,j),2.0));
  
  // DKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAF               
  // DKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAF
	// Here we do the wave-averaging NB: c->eta(i,j) is the FS
	// variance equation with etamean initially unknown
    
	  if(T_sum>=T_INTV_mean)
	  {
    // Compute Significant Wave Height at this time-step
    c->ETAsum(i,j)       += c->eta(i,j)*p->dt;
	  c->ETAmean(i,j)       = c->ETAsum(i,j)/T_sum;
	  c->ETA2sum(i,j)      += c->eta(i,j)*c->eta(i,j);
	  c->ETAvar(i,j)        = 1.0/double(NumDT1-1)*c->ETA2sum(i,j)-c->ETAmean(i,j)*c->ETAmean(i,j)*double(NumDT1)/double(NumDT1-1);
	  c->SigWaveHeight(i,j) = 4.0*sqrt(c->ETAvar(i,j));
	  }
	
	  if(T_sum>=T_INTV_mean)
	  {
	    T_sum  -= T_INTV_mean; // Doesn't need looping	
      NumDT1  = 0;
      SLICELOOP4
      {
	    c->ETAsum(i,j)  = 0.0;
      c->ETA2sum(i,j) = 0.0;
      }
    }	
	// DKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAF
  // DKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAFDKAF
    } //end slice4loop
                 */
}

void fnpf_print_Hs::fill_eta(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{

}

void fnpf_print_Hs::fill_deta(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    
}
    
void fnpf_print_Hs::fill_Uhorz(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f)
{
    
}


void fnpf_print_Hs::ini_location(lexer *p, fdm_fnpf *c)
{

}

int fnpf_print_Hs::conv(double a)
{



}
