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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"rheology_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h" 
#include"diff_void.h"
#include"ediff2.h"
#include"idiff2.h"
#include"idiff2_FS.h"


double rheology_f::Herschel_Bulkley(lexer *p, fdm *a, ghostcell *pgc)
{
	gamma = strainterm(p,a); 
    
    tau0=val=0.0;
    
        if(p->W110==1)
        {
        if(p->W101==0)  // HB
        tau0=p->W96;
        
        if(p->W101==1)  // HB-C dry sand
        tau0 = MAX(0.0,tanphi*(a->press(i,j,k)-p->pressgage) + p->W102_c);
        
        if(p->W101==2)  // HB-C dry sand, without MAX -> issues with negative viscosity and Hypre
        tau0 = (tanphi*(a->press(i,j,k)-p->pressgage) + p->W102_c);
        
        if(p->W101==3)  // HB-C hydrostatic  - MAX added for cells on the interface.
        tau0 = MAX(0.0,tanphi*(a->press(i,j,k)-p->pressgage)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // rho_water = 1000.0, new input?
        
        if(p->W101==4)  // HB-C shear rate generated excess pore pressure
        tau0 = MAX(0.0,tanphi*(a->press(i,j,k)-p->pressgage)*exp(-p->W104*gamma)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k) + p->W102_c);    // m_p is new input W 104 
        
        if(p->W101==5)  // HB-C linear shear rate coupling, max given by pressure
        tau0 = MAX(0.0,tanphi*MAX(0.0,(a->press(i,j,k)-p->pressgage)*MAX(0.0,a->ro(i,j,k)-1000.0)/a->ro(i,j,k)-p->W104*gamma) + p->W102_c);    // m_u also use new input W 104
        
        
        if(p->count==0)
        tau0=p->W96;
        }
    
    
    if(p->W110!=3 || p->W110!=5)
    {
    val =  (tau0/(gamma>1.0e-20?gamma:1.0e-20) + (p->W97)*pow(gamma,p->W98-1.0))/a->ro(i,j,k);
	
	val = MIN(val,p->W95);
    }
	
    return val;
    

}
