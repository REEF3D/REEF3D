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
Author: Dave Kelly
--------------------------------------------------------------------*/

#include"fnpf_print_Hs.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

fnpf_print_Hs::fnpf_print_Hs(lexer *p, slice &Hs) : ETAsum(p), ETAmean(p), //DKAF
                                                    ETA2sum(p), ETAvar(p)
{
    NumDT1=0;      
    T_INTV_mean = 3600.5; // Averaging time for sig wave height
    dT_sum=0; 
    wfcall=0;     
    //T_INTV_mean = 3600.0; // Averaging time for sig. wave height
    wtime=0.0;
    stime = p->P111;        // Start avreging after transients
    
    // Initialise
    wtime  = 0.0;
    T_sum  = 0.0;
    NumDT1 = 0.0;
    
    SLICELOOP4
    {
    ETAsum(i,j)        = 0.0;
    ETAmean(i,j)       = 0.0;
    ETA2sum(i,j)       = 0.0;
    ETAvar(i,j)        = 0.0;
    Hs(i,j)            = 0.0;
    }
}

fnpf_print_Hs::~fnpf_print_Hs()
{
}

void fnpf_print_Hs::start(lexer *p, ghostcell *pgc, slice &eta, slice &Hs)
{
    // RK3

    wtime  += p->dt; //DKAF
    if (wtime>stime) // Check we're past transients
    {
    T_sum  += p->dt; //DKAF
    NumDT1 += 1;    //DKAF
    
    
    SLICELOOP4
    {
	 // Here we do the wave-averaging NB: c->eta(i,j) is the FS
	 // variance equation with etamean initially unknown
      
    ETAsum(i,j)      += eta(i,j)*p->dt;
    ETAmean(i,j)      = ETAsum(i,j)/T_sum;
    ETA2sum(i,j)     += eta(i,j)*eta(i,j);
    
    //cout << "T_sum " << T_sum << " wtim " << wtime <<endl;
    //cin.get();  
    
    //if(T_sum>=T_INTV_mean)
	  //{ 
	    ETAvar(i,j)        = (1.0/double(NumDT1-1))*ETA2sum(i,j)-ETAmean(i,j)*ETAmean(i,j)*(double(NumDT1)/double(NumDT1-1));
	    Hs(i,j)         = 4.0*sqrt(ETAvar(i,j));
    //}
	  
    } //end slice4loop	
    
   }
   
   pgc->gcsl_start4(p,Hs,1);
}


