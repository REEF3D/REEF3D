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


void weno_nug_func::precalc_qf(lexer* p)
{
    // [knox][XN/XP][q_j][coeff]
    
// XN
    IBLOOP
    {
    // imin
    qfx[IP][0][0][0] = ((p->XN[IP2]-p->XN[IP1])/(p->XN[IP3]-p->XN[IP])) 
                     * ((p->XN[IP3]-p->XN[IP1])/(p->XN[IP2]-p->XN[IP]));
                    
    qfx[IP][0][0][1] = ((p->XN[IP2]-p->XN[IP1])/(p->XN[IP3]-p->XN[IP])) 
                     *((p->XN[IP1]-p->XN[IP])/(p->XN[IP3]-p->XN[IP1]));
                          
    
    qfx[IP][0][1][0] = ((p->XN[IP1]-p->XN[IP])/(p->XN[IP2]-p->XN[IM1])) 
                     * ((p->XN[IP1]-p->XN[IM1])/(p->XN[IP2]-p->XN[IP]));
                    
    qfx[IP][0][1][1] = ((p->XN[IP1]-p->XN[IP])/(p->XN[IP2]-p->XN[IM1])) 
                     *((p->XN[IP2]-p->XN[IP1])/(p->XN[IP1]-p->XN[IM1]));
                     
                     
    qfx[IP][0][2][0] = ((p->XN[IP1]-p->XN[IP])/(p->XN[IP]-p->XN[IM2])) 
                     * ((p->XN[IP1]-p->XN[IM1])/(p->XN[IP1]-p->XN[IM2]));
                    
    qfx[IP][0][2][1] = (1.0 + (p->XN[IP1]-p->XN[IP])/(p->XN[IP1]-p->XN[IM1]) 
                            + (p->XN[IP1]-p->XN[IP])/(p->XN[IP1]-p->XN[IM2]));
                     
    //cout<<"qfx_min: "<<qfx[IP][0][0][0]<<" "<<qfx[IP][0][0][1]<<" "<<qfx[IP][0][1][0]<<" "<<qfx[IP][0][1][1]<<" "<<qfx[IP][0][2][0]<<" "<<qfx[IP][0][2][1]<<endl;
    // iplus
    qfx[IP][0][3][0] = (1.0 + (p->XN[IP2]-p->XN[IP1])/(p->XN[IP3]-p->XN[IP1]) 
                            + (p->XN[IP2]-p->XN[IP1])/(p->XN[IP4]-p->XN[IP1]));
                    
    qfx[IP][0][3][1] = ((p->XN[IP2]-p->XN[IP1])/(p->XN[IP4]-p->XN[IP1])) 
                     *((p->XN[IP3]-p->XN[IP1])/(p->XN[IP4]-p->XN[IP2]));
                     
    
    qfx[IP][0][4][0] = ((p->XN[IP2]-p->XN[IP1])/(p->XN[IP3]-p->XN[IP])) 
                     *((p->XN[IP3]-p->XN[IP1])/(p->XN[IP2]-p->XN[IP]));
                    
    qfx[IP][0][4][1] = ((p->XN[IP2]-p->XN[IP1])/(p->XN[IP3]-p->XN[IP])) 
                     *((p->XN[IP1]-p->XN[IP])/(p->XN[IP3]-p->XN[IP1]));
                     
                     
    qfx[IP][0][5][0] = ((p->XN[IP1]-p->XN[IP])/(p->XN[IP2]-p->XN[IM1])) 
                     *((p->XN[IP1]-p->XN[IM1])/(p->XN[IP2]-p->XN[IP]));
                    
    qfx[IP][0][5][1] = ((p->XN[IP1]-p->XN[IP])/(p->XN[IP2]-p->XN[IM1])) 
                     *((p->XN[IP2]-p->XN[IP1])/(p->XN[IP1]-p->XN[IM1]));
                     
                     
    //cout<<"qfx_max: "<<qfx[IP][0][3][0]<<" "<<qfx[IP][0][3][1]<<" "<<qfx[IP][0][4][0]<<" "<<qfx[IP][0][4][1]<<" "<<qfx[IP][0][5][0]<<" "<<qfx[IP][0][5][1]<<endl;
    
    }
    
    
    
    JBLOOP
    {
    // imin
    qfy[JP][0][0][0] = ((p->YN[JP2]-p->YN[JP1])/(p->YN[JP3]-p->YN[JP])) 
                     * ((p->YN[JP3]-p->YN[JP1])/(p->YN[JP2]-p->YN[JP]));
                    
    qfy[JP][0][0][1] = ((p->YN[JP2]-p->YN[JP1])/(p->YN[JP3]-p->YN[JP])) 
                     *((p->YN[JP1]-p->YN[JP])/(p->YN[JP3]-p->YN[JP1]));
                          
    
    qfy[JP][0][1][0] = ((p->YN[JP1]-p->YN[JP])/(p->YN[JP2]-p->YN[JM1])) 
                     * ((p->YN[JP1]-p->YN[JM1])/(p->YN[JP2]-p->YN[JP]));
                    
    qfy[JP][0][1][1] = ((p->YN[JP1]-p->YN[JP])/(p->YN[JP2]-p->YN[JM1])) 
                     *((p->YN[JP2]-p->YN[JP1])/(p->YN[JP1]-p->YN[JM1]));
    
                     
    qfy[JP][0][2][0] = ((p->YN[JP1]-p->YN[JP])/(p->YN[JP]-p->YN[JM2])) 
                     * ((p->YN[JP1]-p->YN[JM1])/(p->YN[JP1]-p->YN[JM2]));
                    
    qfy[JP][0][2][1] = (1.0 + (p->YN[JP1]-p->YN[JP])/(p->YN[JP1]-p->YN[JM1]) 
                            + (p->YN[JP1]-p->YN[JP])/(p->YN[JP1]-p->YN[JM2]));
                     
                     
    // iplus
    qfy[JP][0][3][0] = (1.0 + (p->YN[JP2]-p->YN[JP1])/(p->YN[JP3]-p->YN[JP1]) 
                            + (p->YN[JP2]-p->YN[JP1])/(p->YN[JP4]-p->YN[JP1]));
                    
    qfy[JP][0][3][1] = ((p->YN[JP2]-p->YN[JP1])/(p->YN[JP4]-p->YN[JP1])) 
                     *((p->YN[JP3]-p->YN[JP1])/(p->YN[JP4]-p->YN[JP2]));
                     
    
    qfy[JP][0][4][0] = ((p->YN[JP2]-p->YN[JP1])/(p->YN[JP3]-p->YN[JP])) 
                     *((p->YN[JP3]-p->YN[JP1])/(p->YN[JP2]-p->YN[JP]));
                    
    qfy[JP][0][4][1] = ((p->YN[JP2]-p->YN[JP1])/(p->YN[JP3]-p->YN[JP])) 
                     *((p->YN[JP1]-p->YN[JP])/(p->YN[JP3]-p->YN[JP1]));
                     
                     
    qfy[JP][0][5][0] = ((p->YN[JP1]-p->YN[JP])/(p->YN[JP2]-p->YN[JM1])) 
                     *((p->YN[JP1]-p->YN[JM1])/(p->YN[JP2]-p->YN[JP]));
                    
    qfy[JP][0][5][1] = ((p->YN[JP1]-p->YN[JP])/(p->YN[JP2]-p->YN[JM1])) 
                     *((p->YN[JP2]-p->YN[JP1])/(p->YN[JP1]-p->YN[JM1]));
    
    }
    
    
    
    KBLOOP
    {
    // imin
    qfz[KP][0][0][0] = ((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KP])) 
                     * ((p->ZN[KP3]-p->ZN[KP1])/(p->ZN[KP2]-p->ZN[KP]));
                    
    qfz[KP][0][0][1] = ((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KP])) 
                     *((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP3]-p->ZN[KP1]));
                          
    
    qfz[KP][0][1][0] = ((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP2]-p->ZN[KM1])) 
                     * ((p->ZN[KP1]-p->ZN[KM1])/(p->ZN[KP2]-p->ZN[KP]));
                    
    qfz[KP][0][1][1] = ((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP2]-p->ZN[KM1])) 
                     *((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP1]-p->ZN[KM1]));
    
    
    qfz[KP][0][2][0] = ((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP]-p->ZN[KM2])) 
                     * ((p->ZN[KP1]-p->ZN[KM1])/(p->ZN[KP1]-p->ZN[KM2]));
                    
    qfz[KP][0][2][1] = (1.0 + (p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP1]-p->ZN[KM1]) 
                            + (p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP1]-p->ZN[KM2]));
                     
                     
    // iplus
    qfz[KP][0][3][0] = (1.0 + (p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KP1]) 
                            + (p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP4]-p->ZN[KP1]));
                    
    qfz[KP][0][3][1] = ((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP4]-p->ZN[KP1])) 
                     *((p->ZN[KP3]-p->ZN[KP1])/(p->ZN[KP4]-p->ZN[KP2]));
                     
    
    qfz[KP][0][4][0] = ((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KP])) 
                     *((p->ZN[KP3]-p->ZN[KP1])/(p->ZN[KP2]-p->ZN[KP]));
                    
    qfz[KP][0][4][1] = ((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KP])) 
                     *((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP3]-p->ZN[KP1]));
                     
                     
    qfz[KP][0][5][0] = ((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP2]-p->ZN[KM1])) 
                     *((p->ZN[KP1]-p->ZN[KM1])/(p->ZN[KP2]-p->ZN[KP]));
                    
    qfz[KP][0][5][1] = ((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP2]-p->ZN[KM1])) 
                     *((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP1]-p->ZN[KM1]));
    }
    
//-------------------------------------------------------------------------
// XP
    IBLOOP
    {
    // imin
    qfx[IP][1][0][0] = ((p->XP[IP2]-p->XP[IP1])/(p->XP[IP3]-p->XP[IP])) 
                     * ((p->XP[IP3]-p->XP[IP1])/(p->XP[IP2]-p->XP[IP]));
                    
    qfx[IP][1][0][1] = ((p->XP[IP2]-p->XP[IP1])/(p->XP[IP3]-p->XP[IP])) 
                     *((p->XP[IP1]-p->XP[IP])/(p->XP[IP3]-p->XP[IP1]));
                          
    
    qfx[IP][1][1][0] = ((p->XP[IP1]-p->XP[IP])/(p->XP[IP2]-p->XP[IM1])) 
                     * ((p->XP[IP1]-p->XP[IM1])/(p->XP[IP2]-p->XP[IP]));
                    
    qfx[IP][1][1][1] = ((p->XP[IP1]-p->XP[IP])/(p->XP[IP2]-p->XP[IM1])) 
                     *((p->XP[IP2]-p->XP[IP1])/(p->XP[IP1]-p->XP[IM1]));
    

    qfx[IP][1][2][0] = ((p->XP[IP1]-p->XP[IP])/(p->XP[IP]-p->XP[IM2])) 
                     * ((p->XP[IP1]-p->XP[IM1])/(p->XP[IP1]-p->XP[IM2]));
                    
    qfx[IP][1][2][1] = (1.0 + (p->XP[IP1]-p->XP[IP])/(p->XP[IP1]-p->XP[IM1]) 
                            + (p->XP[IP1]-p->XP[IP])/(p->XP[IP1]-p->XP[IM2]));
                     
                     
    // iplus
    qfx[IP][1][3][0] = (1.0 + (p->XP[IP2]-p->XP[IP1])/(p->XP[IP3]-p->XP[IP1]) 
                            + (p->XP[IP2]-p->XP[IP1])/(p->XP[IP4]-p->XP[IP1]));
                    
    qfx[IP][1][3][1] = ((p->XP[IP2]-p->XP[IP1])/(p->XP[IP4]-p->XP[IP1])) 
                     *((p->XP[IP3]-p->XP[IP1])/(p->XP[IP4]-p->XP[IP2]));
                     
    
    qfx[IP][1][4][0] = ((p->XP[IP2]-p->XP[IP1])/(p->XP[IP3]-p->XP[IP])) 
                     *((p->XP[IP3]-p->XP[IP1])/(p->XP[IP2]-p->XP[IP]));
                    
    qfx[IP][1][4][1] = ((p->XP[IP2]-p->XP[IP1])/(p->XP[IP3]-p->XP[IP])) 
                     *((p->XP[IP1]-p->XP[IP])/(p->XP[IP3]-p->XP[IP1]));
                     
                     
    qfx[IP][1][5][0] = ((p->XP[IP1]-p->XP[IP])/(p->XP[IP2]-p->XP[IM1])) 
                     *((p->XP[IP1]-p->XP[IM1])/(p->XP[IP2]-p->XP[IP]));
                    
    qfx[IP][1][5][1] = ((p->XP[IP1]-p->XP[IP])/(p->XP[IP2]-p->XP[IM1])) 
                     *((p->XP[IP2]-p->XP[IP1])/(p->XP[IP1]-p->XP[IM1]));
    
    }
    
    
    
    JBLOOP
    {
    // imin
    qfy[JP][1][0][0] = ((p->YP[JP2]-p->YP[JP1])/(p->YP[JP3]-p->YP[JP])) 
                     * ((p->YP[JP3]-p->YP[JP1])/(p->YP[JP2]-p->YP[JP]));
                    
    qfy[JP][1][0][1] = ((p->YP[JP2]-p->YP[JP1])/(p->YP[JP3]-p->YP[JP])) 
                     *((p->YP[JP1]-p->YP[JP])/(p->YP[JP3]-p->YP[JP1]));
                          
    
    qfy[JP][1][1][0] = ((p->YP[JP1]-p->YP[JP])/(p->YP[JP2]-p->YP[JM1])) 
                     * ((p->YP[JP1]-p->YP[JM1])/(p->YP[JP2]-p->YP[JP]));
                    
    qfy[JP][1][1][1] = ((p->YP[JP1]-p->YP[JP])/(p->YP[JP2]-p->YP[JM1])) 
                     *((p->YP[JP2]-p->YP[JP1])/(p->YP[JP1]-p->YP[JM1]));
    
    
    qfy[JP][1][2][0] = ((p->YP[JP1]-p->YP[JP])/(p->YP[JP]-p->YP[JM2])) 
                     * ((p->YP[JP1]-p->YP[JM1])/(p->YP[JP1]-p->YP[JM2]));
                    
    qfy[JP][1][2][1] = (1.0 + (p->YP[JP1]-p->YP[JP])/(p->YP[JP1]-p->YP[JM1]) 
                            + (p->YP[JP1]-p->YP[JP])/(p->YP[JP1]-p->YP[JM2]));
                     
                     
    // iplus
    qfy[JP][1][3][0] = (1.0 + (p->YP[JP2]-p->YP[JP1])/(p->YP[JP3]-p->YP[JP1]) 
                            + (p->YP[JP2]-p->YP[JP1])/(p->YP[JP4]-p->YP[JP1]));
                    
    qfy[JP][1][3][1] = ((p->YP[JP2]-p->YP[JP1])/(p->YP[JP4]-p->YP[JP1])) 
                     *((p->YP[JP3]-p->YP[JP1])/(p->YP[JP4]-p->YP[JP2]));
                     
    
    qfy[JP][1][4][0] = ((p->YP[JP2]-p->YP[JP1])/(p->YP[JP3]-p->YP[JP])) 
                     *((p->YP[JP3]-p->YP[JP1])/(p->YP[JP2]-p->YP[JP]));
                    
    qfy[JP][1][4][1] = ((p->YP[JP2]-p->YP[JP1])/(p->YP[JP3]-p->YP[JP])) 
                     *((p->YP[JP1]-p->YP[JP])/(p->YP[JP3]-p->YP[JP1]));
                     
                     
    qfy[JP][1][5][0] = ((p->YP[JP1]-p->YP[JP])/(p->YP[JP2]-p->YP[JM1])) 
                     *((p->YP[JP1]-p->YP[JM1])/(p->YP[JP2]-p->YP[JP]));
                    
    qfy[JP][1][5][1] = ((p->YP[JP1]-p->YP[JP])/(p->YP[JP2]-p->YP[JM1])) 
                     *((p->YP[JP2]-p->YP[JP1])/(p->YP[JP1]-p->YP[JM1]));
    
    }
    
    
    
    KBLOOP
    {
    // imin
    qfz[KP][1][0][0] = ((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KP])) 
                     * ((p->ZP[KP3]-p->ZP[KP1])/(p->ZP[KP2]-p->ZP[KP]));
                    
    qfz[KP][1][0][1] = ((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KP])) 
                     *((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP3]-p->ZP[KP1]));
                     
    
    qfz[KP][1][1][0] = ((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP2]-p->ZP[KM1])) 
                     * ((p->ZP[KP1]-p->ZP[KM1])/(p->ZP[KP2]-p->ZP[KP]));
                    
    qfz[KP][1][1][1] = ((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP2]-p->ZP[KM1])) 
                     *((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP1]-p->ZP[KM1]));
    
    
    qfz[KP][1][2][0] = ((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP]-p->ZP[KM2])) 
                     * ((p->ZP[KP1]-p->ZP[KM1])/(p->ZP[KP1]-p->ZP[KM2]));
                    
    qfz[KP][1][2][1] = (1.0 + (p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP1]-p->ZP[KM1]) 
                            + (p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP1]-p->ZP[KM2]));
                     
                     
    // iplus
    qfz[KP][1][3][0] = (1.0 + (p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KP1]) 
                            + (p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP4]-p->ZP[KP1]));
                    
    qfz[KP][1][3][1] = ((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP4]-p->ZP[KP1])) 
                     *((p->ZP[KP3]-p->ZP[KP1])/(p->ZP[KP4]-p->ZP[KP2]));
                     
    
    qfz[KP][1][4][0] = ((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KP])) 
                     *((p->ZP[KP3]-p->ZP[KP1])/(p->ZP[KP2]-p->ZP[KP]));
                    
    qfz[KP][1][4][1] = ((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KP])) 
                     *((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP3]-p->ZP[KP1]));
                     
                     
    qfz[KP][1][5][0] = ((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP2]-p->ZP[KM1])) 
                     *((p->ZP[KP1]-p->ZP[KM1])/(p->ZP[KP2]-p->ZP[KP]));
                    
    qfz[KP][1][5][1] = ((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP2]-p->ZP[KM1])) 
                     *((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP1]-p->ZP[KM1]));
    }
    
}




