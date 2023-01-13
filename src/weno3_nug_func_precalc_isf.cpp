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

#include"weno3_nug_func.h"
#include"lexer.h"
#include"fdm.h"

void weno3_nug_func::precalc_isf(lexer* p)
{
    
// XN
    IBLOOP
    {
    // imin

    // is1    
    isfx[IP][0][0] = 4.0*pow((p->XN[IP1]-p->XN[IP])/(p->XN[IP2]-p->XN[IP]), 2.0);
                    
    // is2
    isfx[IP][0][1] = 4.0*pow((p->XN[IP1]-p->XN[IP])/(p->XN[IP1]-p->XN[IM1]), 2.0);
                      
                      
    // imax ---------
    
    // is1    
    isfx[IP][0][2] = 4.0*pow((p->XN[IP2]-p->XN[IP1])/(p->XN[IP3]-p->XN[IP1]), 2.0);
                    
    // is2
    isfx[IP][0][3] = 4.0*pow((p->XN[IP2]-p->XN[IP1])/(p->XN[IP2]-p->XN[IP]), 2.0);
    }
     
// YN
    JBLOOP
    {
    // imin

    // is1    
    isfy[JP][0][0] = 4.0*pow((p->YN[JP1]-p->YN[JP])/(p->YN[JP2]-p->YN[JP]), 2.0);
                    
    // is2
    isfy[JP][0][1] = 4.0*pow((p->YN[JP1]-p->YN[JP])/(p->YN[JP1]-p->YN[JM1]), 2.0);
                      
                      
    // imax ---------
    
    // is1    
    isfy[JP][0][2] = 4.0*pow((p->YN[JP2]-p->YN[JP1])/(p->YN[JP3]-p->YN[JP1]), 2.0);
                    
    // is2
    isfy[JP][0][3] = 4.0*pow((p->YN[JP2]-p->YN[JP1])/(p->YN[JP2]-p->YN[JP]), 2.0);
    }
    
// ZN
    KBLOOP
    {
    // imin

    // is1    
    isfz[KP][0][0] = 4.0*pow((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP2]-p->ZN[KP]), 2.0);
                    
    // is2
    isfz[KP][0][1] = 4.0*pow((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP1]-p->ZN[KM1]), 2.0);
                      
                      
    // imax ---------
    
    // is1    
    isfz[KP][0][2] = 4.0*pow((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KP1]), 2.0);
                    
    // is2
    isfz[KP][0][3] = 4.0*pow((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP2]-p->ZN[KP]), 2.0);
    }
    
// ---------------------------------------------------------------------

// XP
    IBLOOP
    {
    // imin

    // is1    
    isfx[IP][1][0] = 4.0*pow((p->XP[IP1]-p->XP[IP])/(p->XP[IP2]-p->XP[IP]), 2.0);
                    
    // is2
    isfx[IP][1][1] = 4.0*pow((p->XP[IP1]-p->XP[IP])/(p->XP[IP1]-p->XP[IM1]), 2.0);
                      
                      
    // imax ---------
    
    // is1    
    isfx[IP][1][2] = 4.0*pow((p->XP[IP2]-p->XP[IP1])/(p->XP[IP3]-p->XP[IP1]), 2.0);
                    
    // is2
    isfx[IP][1][3] = 4.0*pow((p->XP[IP2]-p->XP[IP1])/(p->XP[IP2]-p->XP[IP]), 2.0);
    }
     
// YP
    JBLOOP
    {
    // imin

    // is1    
    isfy[JP][1][0] = 4.0*pow((p->YP[JP1]-p->YP[JP])/(p->YP[JP2]-p->YP[JP]), 2.0);
                    
    // is2
    isfy[JP][1][1] = 4.0*pow((p->YP[JP1]-p->YP[JP])/(p->YP[JP1]-p->YP[JM1]), 2.0);
                      
                      
    // imax ---------
    
    // is1    
    isfy[JP][1][2] = 4.0*pow((p->YP[JP2]-p->YP[JP1])/(p->YP[JP3]-p->YP[JP1]), 2.0);
                    
    // is2
    isfy[JP][1][3] = 4.0*pow((p->YP[JP2]-p->YP[JP1])/(p->YP[JP2]-p->YP[JP]), 2.0);
    }
    
// ZP
    KBLOOP
    {
    // imin

    // is1    
    isfz[KP][1][0] = 4.0*pow((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP2]-p->ZP[KP]), 2.0);
                    
    // is2
    isfz[KP][1][1] = 4.0*pow((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP1]-p->ZP[KM1]), 2.0);
                      
                      
    // imax ---------
    
    // is1    
    isfz[KP][1][2] = 4.0*pow((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KP1]), 2.0);
                    
    // is2
    isfz[KP][1][3] = 4.0*pow((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP2]-p->ZP[KP]), 2.0);
    }
}
