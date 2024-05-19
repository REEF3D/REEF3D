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

#include"weno_nug_func.h"
#include"lexer.h"
#include"fdm.h"

void weno_nug_func::precalc_cf(lexer* p)
{
    // array order
    // [knox][XN/XP][q_j][coeff]
     
// XN
    IBLOOP
    {
    // imin
    cfx[IP][0][0] = ((p->XN[IP1]-p->XN[IM2])/(p->XN[IP3]-p->XN[IM2])) 
                  * ((p->XN[IP1]-p->XN[IM1])/(p->XN[IP3]-p->XN[IM1]));
            
    cfx[IP][0][1] = ((p->XN[IP1]-p->XN[IM2])/(p->XN[IP3]-p->XN[IM2])) 
                  * ((p->XN[IP3]-p->XN[IP1])/(p->XN[IP3]-p->XN[IM1]))
                  * ((p->XN[IP3]-p->XN[IM1])/(p->XN[IP2]-p->XN[IM2])+1.0);
                  
    cfx[IP][0][2] = ((p->XN[IP2]-p->XN[IP1])/(p->XN[IP3]-p->XN[IM2])) 
                  * ((p->XN[IP3]-p->XN[IP1])/(p->XN[IP2]-p->XN[IM2]));
                  
    // imax
    cfx[IP][0][3] = ((p->XN[IP1]-p->XN[IP])/(p->XN[IP4]-p->XN[IM1])) 
                  * ((p->XN[IP1]-p->XN[IM1])/(p->XN[IP4]-p->XN[IP]));
                  
    cfx[IP][0][4] = ((p->XN[IP1]-p->XN[IM1])/(p->XN[IP4]-p->XN[IM1])) 
                  * ((p->XN[IP4]-p->XN[IP1])/(p->XN[IP4]-p->XN[IP]))
                  * ((p->XN[IP4]-p->XN[IP])/(p->XN[IP3]-p->XN[IM1])+1.0);
                  
    cfx[IP][0][5] = ((p->XN[IP3]-p->XN[IP1])/(p->XN[IP4]-p->XN[IM1])) 
                  * ((p->XN[IP4]-p->XN[IP1])/(p->XN[IP3]-p->XN[IM1]));
    }
                          
    JBLOOP
    {
    // jmin
    cfy[JP][0][0] = ((p->YN[JP1]-p->YN[JM2])/(p->YN[JP3]-p->YN[JM2])) 
                  * ((p->YN[JP1]-p->YN[JM1])/(p->YN[JP3]-p->YN[JM1]));
                    
    cfy[JP][0][1] = ((p->YN[JP1]-p->YN[JM2])/(p->YN[JP3]-p->YN[JM2])) 
                  * ((p->YN[JP3]-p->YN[JP1])/(p->YN[JP3]-p->YN[JM1]))
                  * ((p->YN[JP3]-p->YN[JM1])/(p->YN[JP2]-p->YN[JM2])+1.0);
                  
    cfy[JP][0][2] = ((p->YN[JP2]-p->YN[JP1])/(p->YN[JP3]-p->YN[JM2])) 
                  * ((p->YN[JP3]-p->YN[JP1])/(p->YN[JP2]-p->YN[JM2]));
                  
    // jmax
    cfy[JP][0][3] = ((p->YN[JP1]-p->YN[JP])/(p->YN[JP4]-p->YN[JM1])) 
                  * ((p->YN[JP1]-p->YN[JM1])/(p->YN[JP4]-p->YN[JP]));
                  
    cfy[JP][0][4] = ((p->YN[JP1]-p->YN[JM1])/(p->YN[JP4]-p->YN[JM1])) 
                  * ((p->YN[JP4]-p->YN[JP1])/(p->YN[JP4]-p->YN[JP]))
                  * ((p->YN[JP4]-p->YN[JP])/(p->YN[JP3]-p->YN[JM1])+1.0);
                  
    cfy[JP][0][5] = ((p->YN[JP3]-p->YN[JP1])/(p->YN[JP4]-p->YN[JM1])) 
                  * ((p->YN[JP4]-p->YN[JP1])/(p->YN[JP3]-p->YN[JM1]));
    }
    
    
    KBLOOP
    {
    // kmin
    cfz[KP][0][0] = ((p->ZN[KP1]-p->ZN[KM2])/(p->ZN[KP3]-p->ZN[KM2])) 
                  * ((p->ZN[KP1]-p->ZN[KM1])/(p->ZN[KP3]-p->ZN[KM1]));
                    
    cfz[KP][0][1] = ((p->ZN[KP1]-p->ZN[KM2])/(p->ZN[KP3]-p->ZN[KM2])) 
                  * ((p->ZN[KP3]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KM1]))
                  * ((p->ZN[KP3]-p->ZN[KM1])/(p->ZN[KP2]-p->ZN[KM2])+1.0);
                  
    cfz[KP][0][2] = ((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KM2])) 
                  * ((p->ZN[KP3]-p->ZN[KP1])/(p->ZN[KP2]-p->ZN[KM2]));
                                 
    // kmax
    cfz[KP][0][3] = ((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP4]-p->ZN[KM1])) 
                  * ((p->ZN[KP1]-p->ZN[KM1])/(p->ZN[KP4]-p->ZN[KP]));
                  
    cfz[KP][0][4] = ((p->ZN[KP1]-p->ZN[KM1])/(p->ZN[KP4]-p->ZN[KM1])) 
                  * ((p->ZN[KP4]-p->ZN[KP1])/(p->ZN[KP4]-p->ZN[KP]))
                  * ((p->ZN[KP4]-p->ZN[KP])/(p->ZN[KP3]-p->ZN[KM1])+1.0);
                  
    cfz[KP][0][5] = ((p->ZN[KP3]-p->ZN[KP1])/(p->ZN[KP4]-p->ZN[KM1])) 
                  * ((p->ZN[KP4]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KM1]));
    }
    
//-------------------------------------------------------------------------
// XP

    IBLOOP
    {
    // imin
    cfx[IP][1][0] = ((p->XP[IP1]-p->XP[IM2])/(p->XP[IP3]-p->XP[IM2])) 
                  * ((p->XP[IP1]-p->XP[IM1])/(p->XP[IP3]-p->XP[IM1]));
                  
    cfx[IP][1][1] = ((p->XP[IP1]-p->XP[IM2])/(p->XP[IP3]-p->XP[IM2])) 
                  * ((p->XP[IP3]-p->XP[IP1])/(p->XP[IP3]-p->XP[IM1]))
                  * ((p->XP[IP3]-p->XP[IM1])/(p->XP[IP2]-p->XP[IM2])+1.0);
                  
    cfx[IP][1][2] = ((p->XP[IP2]-p->XP[IP1])/(p->XP[IP3]-p->XP[IM2])) 
                  * ((p->XP[IP3]-p->XP[IP1])/(p->XP[IP2]-p->XP[IM2]));
                    
                  
    // imax
    cfx[IP][1][3] = ((p->XP[IP1]-p->XP[IP])/(p->XP[IP4]-p->XP[IM1])) 
                  * ((p->XP[IP1]-p->XP[IM1])/(p->XP[IP4]-p->XP[IP]));
                  
    cfx[IP][1][4] = ((p->XP[IP1]-p->XP[IM1])/(p->XP[IP4]-p->XP[IM1])) 
                  * ((p->XP[IP4]-p->XP[IP1])/(p->XP[IP4]-p->XP[IP]))
                  * ((p->XP[IP4]-p->XP[IP])/(p->XP[IP3]-p->XP[IM1])+1.0);
                  
    cfx[IP][1][5] = ((p->XP[IP3]-p->XP[IP1])/(p->XP[IP4]-p->XP[IM1])) 
                  * ((p->XP[IP4]-p->XP[IP1])/(p->XP[IP3]-p->XP[IM1]));
    }
                          
    JBLOOP
    {
    // jmin
    cfy[JP][1][0] = ((p->YP[JP1]-p->YP[JM2])/(p->YP[JP3]-p->YP[JM2])) 
                  * ((p->YP[JP1]-p->YP[JM1])/(p->YP[JP3]-p->YP[JM1]));
                    
    cfy[JP][1][1] = ((p->YP[JP1]-p->YP[JM2])/(p->YP[JP3]-p->YP[JM2])) 
                  * ((p->YP[JP3]-p->YP[JP1])/(p->YP[JP3]-p->YP[JM1]))
                  * ((p->YP[JP3]-p->YP[JM1])/(p->YP[JP2]-p->YP[JM2])+1.0);
                  
    cfy[JP][1][2] = ((p->YP[JP2]-p->YP[JP1])/(p->YP[JP3]-p->YP[JM2])) 
                  * ((p->YP[JP3]-p->YP[JP1])/(p->YP[JP2]-p->YP[JM2]));
                  
    // jmax
    cfy[JP][1][3] = ((p->YP[JP1]-p->YP[JP])/(p->YP[JP4]-p->YP[JM1])) 
                  * ((p->YP[JP1]-p->YP[JM1])/(p->YP[JP4]-p->YP[JP]));
                  
    cfy[JP][1][4] = ((p->YP[JP1]-p->YP[JM1])/(p->YP[JP4]-p->YP[JM1])) 
                  * ((p->YP[JP4]-p->YP[JP1])/(p->YP[JP4]-p->YP[JP]))
                  * ((p->YP[JP4]-p->YP[JP])/(p->YP[JP3]-p->YP[JM1])+1.0);
                  
    cfy[JP][1][5] = ((p->YP[JP3]-p->YP[JP1])/(p->YP[JP4]-p->YP[JM1])) 
                  * ((p->YP[JP4]-p->YP[JP1])/(p->YP[JP3]-p->YP[JM1]));
    }
    
    
    KBLOOP
    {
    // kmin
    cfz[KP][1][0] = ((p->ZP[KP1]-p->ZP[KM2])/(p->ZP[KP3]-p->ZP[KM2])) 
                  * ((p->ZP[KP1]-p->ZP[KM1])/(p->ZP[KP3]-p->ZP[KM1]));
                    
    cfz[KP][1][1] = ((p->ZP[KP1]-p->ZP[KM2])/(p->ZP[KP3]-p->ZP[KM2])) 
                  * ((p->ZP[KP3]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KM1]))
                  * ((p->ZP[KP3]-p->ZP[KM1])/(p->ZP[KP2]-p->ZP[KM2])+1.0);
                  
    cfz[KP][1][2] = ((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KM2])) 
                  * ((p->ZP[KP3]-p->ZP[KP1])/(p->ZP[KP2]-p->ZP[KM2]));
                  
    // kmax
    cfz[KP][1][3] = ((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP4]-p->ZP[KM1])) 
                  * ((p->ZP[KP1]-p->ZP[KM1])/(p->ZP[KP4]-p->ZP[KP]));
                  
    cfz[KP][1][4] = ((p->ZP[KP1]-p->ZP[KM1])/(p->ZP[KP4]-p->ZP[KM1])) 
                  * ((p->ZP[KP4]-p->ZP[KP1])/(p->ZP[KP4]-p->ZP[KP]))
                  * ((p->ZP[KP4]-p->ZP[KP])/(p->ZP[KP3]-p->ZP[KM1])+1.0);
                  
    cfz[KP][1][5] = ((p->ZP[KP3]-p->ZP[KP1])/(p->ZP[KP4]-p->ZP[KM1])) 
                  * ((p->ZP[KP4]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KM1]));
    }
}




