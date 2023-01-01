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

#include"weno_nug_func.h"
#include"lexer.h"
#include"fdm.h"

void weno_nug_func::precalc_isf(lexer* p)
{
    double fac;
    
    // [knox][XN/XP][is_j][coeff]
    
// XN
    IBLOOP
    {
    // imin

    // is1
    fac = 4.0*pow((p->XN[IP1]-p->XN[IP])/(p->XN[IP3]-p->XN[IP]), 2.0);
    
    isfx[IP][0][0][0] = fac*((10.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         + (p->XN[IP2]-p->XN[IP])*(p->XN[IP2]-p->XN[IP1]))
    
                      / pow(p->XN[IP3]-p->XN[IP1], 2.0));
                    
                    
    isfx[IP][0][0][1] = fac*((20.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         + 2.0*(p->XN[IP2]-p->XN[IP])*(p->XN[IP2]-p->XN[IP1])
                         + (p->XN[IP3]-p->XN[IP])*(2.0*p->XN[IP2]-p->XN[IP1]-p->XN[IP]))
    
                      / ((p->XN[IP3]-p->XN[IP1])*(p->XN[IP2]-p->XN[IP])));
                  
                  
    isfx[IP][0][0][2] = fac*((10.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         + (p->XN[IP3]+p->XN[IP2]-2.0*p->XN[IP])
                            *(p->XN[IP3]+p->XN[IP2]-p->XN[IP1]-p->XN[IP]))

                      / pow(p->XN[IP2]-p->XN[IP], 2.0));
                      
                      
    // is2
    fac = 4.0*pow((p->XN[IP1]-p->XN[IP])/(p->XN[IP2]-p->XN[IM1]), 2.0);
    
    isfx[IP][0][1][0] = fac*((10.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         + (p->XN[IP2]-p->XN[IP])*(p->XN[IP2]-p->XN[IP1]))
    
                      / pow(p->XN[IP1]-p->XN[IM1], 2.0));
                    
                    
    isfx[IP][0][1][1] = fac*((20.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         - (p->XN[IP2]-p->XN[IP1])*(p->XN[IP]-p->XN[IM1])
                         - (p->XN[IP2]-p->XN[IP])*(p->XN[IP1]-p->XN[IM1]))
    
                      / ((p->XN[IP2]-p->XN[IP])*(p->XN[IP1]-p->XN[IM1])));
                  
                  
    isfx[IP][0][1][2] = fac*((10.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         + (p->XN[IP]-p->XN[IM1])*(p->XN[IP1]-p->XN[IM1]))

                      / pow(p->XN[IP2]-p->XN[IP], 2.0));
                      
                      
    // is3
    fac = 4.0*pow((p->XN[IP1]-p->XN[IP])/(p->XN[IP1]-p->XN[IM2]), 2.0);
    
    isfx[IP][0][2][0] = fac*((10.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         + (p->XN[IP1]-p->XN[IM1])*(p->XN[IP]-p->XN[IM1]))
    
                      / pow(p->XN[IP]-p->XN[IM2], 2.0));
                    
                    
    isfx[IP][0][2][1] = fac*((20.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         + 2.0*(p->XN[IP1]-p->XN[IM1])*(p->XN[IP]-p->XN[IM1])
                         + (p->XN[IP1]-p->XN[IM2])*(p->XN[IP1]+p->XN[IP]-2.0*p->XN[IM1]))
    
                      / ((p->XN[IP1]-p->XN[IM1])*(p->XN[IP]-p->XN[IM2])));
                  
                  
    isfx[IP][0][2][2] = fac*((10.0*pow(p->XN[IP1]-p->XN[IP], 2.0) 
                         + (2.0*p->XN[IP1]-p->XN[IM2]-p->XN[IM1])
                            *(p->XN[IP1]+p->XN[IP]-p->XN[IM1]-p->XN[IM2]))

                      / pow(p->XN[IP1]-p->XN[IM1], 2.0));
                      
                      
                      
                      
    // imax ---------
    
    // is1
    fac = 4.0*pow((p->XN[IP2]-p->XN[IP1])/(p->XN[IP4]-p->XN[IP1]), 2.0);
    
    isfx[IP][0][3][0] = fac*((10.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         + (p->XN[IP3]-p->XN[IP1])*(p->XN[IP3]-p->XN[IP2]))
    
                      / pow(p->XN[IP4]-p->XN[IP2], 2.0));
                    
                    
    isfx[IP][0][3][1] = fac*((20.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         + 2.0*(p->XN[IP3]-p->XN[IP1])*(p->XN[IP3]-p->XN[IP2])
                         + (p->XN[IP4]-p->XN[IP1])*(2.0*p->XN[IP3]-p->XN[IP2]-p->XN[IP1]))
    
                      / ((p->XN[IP4]-p->XN[IP2])*(p->XN[IP3]-p->XN[IP1])));
                  
                  
    isfx[IP][0][3][2] = fac*((10.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         + (p->XN[IP4]+p->XN[IP3]-2.0*p->XN[IP1])
                            *(p->XN[IP4]+p->XN[IP3]-p->XN[IP2]-p->XN[IP1]))

                      / pow(p->XN[IP3]-p->XN[IP1], 2.0));
                      
            
    // is2
    fac = 4.0*pow((p->XN[IP2]-p->XN[IP1])/(p->XN[IP3]-p->XN[IP]), 2.0);
    
    isfx[IP][0][4][0] = fac*((10.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         + (p->XN[IP3]-p->XN[IP1])*(p->XN[IP3]-p->XN[IP2]))
    
                      / pow(p->XN[IP2]-p->XN[IP], 2.0));
                    
                    
    isfx[IP][0][4][1] = fac*((20.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         - (p->XN[IP3]-p->XN[IP2])*(p->XN[IP1]-p->XN[IP])
                         - (p->XN[IP3]-p->XN[IP1])*(p->XN[IP2]-p->XN[IP]))
    
                      / ((p->XN[IP3]-p->XN[IP1])*(p->XN[IP2]-p->XN[IP])));
                  
                  
    isfx[IP][0][4][2] = fac*((10.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         + (p->XN[IP1]-p->XN[IP])*(p->XN[IP2]-p->XN[IP]))

                      / pow(p->XN[IP3]-p->XN[IP1], 2.0));
                      
                      
    // is3
    fac = 4.0*pow((p->XN[IP2]-p->XN[IP1])/(p->XN[IP2]-p->XN[IM1]), 2.0);
    
    isfx[IP][0][5][0] = fac*((10.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         + (p->XN[IP2]-p->XN[IP])*(p->XN[IP1]-p->XN[IP]))
    
                      / pow(p->XN[IP1]-p->XN[IM1], 2.0));
                    
                    
    isfx[IP][0][5][1] = fac*((20.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         + 2.0*(p->XN[IP2]-p->XN[IP])*(p->XN[IP1]-p->XN[IP])
                         + (p->XN[IP2]-p->XN[IM1])*(p->XN[IP2]+p->XN[IP1]-2.0*p->XN[IP]))
    
                      / ((p->XN[IP2]-p->XN[IP])*(p->XN[IP1]-p->XN[IM1])));
                  
                  
    isfx[IP][0][5][2] = fac*((10.0*pow(p->XN[IP2]-p->XN[IP1], 2.0) 
                         + (2.0*p->XN[IP2]-p->XN[IM1]-p->XN[IP])
                            *(p->XN[IP2]+p->XN[IP1]-p->XN[IP]-p->XN[IM1]))

                      / pow(p->XN[IP2]-p->XN[IP], 2.0));

    }
     
// YN
    JBLOOP
    {
    // jmin
    
    // is1
    fac = 4.0*pow((p->YN[JP1]-p->YN[JP])/(p->YN[JP3]-p->YN[JP]), 2.0);
    
    isfy[JP][0][0][0] = fac*((10.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         + (p->YN[JP2]-p->YN[JP])*(p->YN[JP2]-p->YN[JP1]))
    
                      / pow(p->YN[JP3]-p->YN[JP1], 2.0));
                    
                    
    isfy[JP][0][0][1] = fac*((20.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         + 2.0*(p->YN[JP2]-p->YN[JP])*(p->YN[JP2]-p->YN[JP1])
                         + (p->YN[JP3]-p->YN[JP])*(2.0*p->YN[JP2]-p->YN[JP1]-p->YN[JP]))
    
                      / ((p->YN[JP3]-p->YN[JP1])*(p->YN[JP2]-p->YN[JP])));
                  
                  
    isfy[JP][0][0][2] = fac*((10.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         + (p->YN[JP3]+p->YN[JP2]-2.0*p->YN[JP])
                            *(p->YN[JP3]+p->YN[JP2]-p->YN[JP1]-p->YN[JP]))

                      / pow(p->YN[JP2]-p->YN[JP], 2.0));
                      
                      
    // is2
    fac = 4.0*pow((p->YN[JP1]-p->YN[JP])/(p->YN[JP2]-p->YN[JM1]), 2.0);
    
    isfy[JP][0][1][0] = fac*((10.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         + (p->YN[JP2]-p->YN[JP])*(p->YN[JP2]-p->YN[JP1]))
    
                      / pow(p->YN[JP1]-p->YN[JM1], 2.0));
                    
                    
    isfy[JP][0][1][1] = fac*((20.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         - (p->YN[JP2]-p->YN[JP1])*(p->YN[JP]-p->YN[JM1])
                         - (p->YN[JP2]-p->YN[JP])*(p->YN[JP1]-p->YN[JM1]))
    
                      / ((p->YN[JP2]-p->YN[JP])*(p->YN[JP1]-p->YN[JM1])));
                  
                  
    isfy[JP][0][1][2] = fac*((10.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         + (p->YN[JP]-p->YN[JM1])*(p->YN[JP1]-p->YN[JM1]))

                      / pow(p->YN[JP2]-p->YN[JP], 2.0));
                      
                      
     // is3
    fac = 4.0*pow((p->YN[JP1]-p->YN[JP])/(p->YN[JP1]-p->YN[JM2]), 2.0);
    
    isfy[JP][0][2][0] = fac*((10.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         + (p->YN[JP1]-p->YN[JM1])*(p->YN[JP]-p->YN[JM1]))
    
                      / pow(p->YN[JP]-p->YN[JM2], 2.0));
                    
                    
    isfy[JP][0][2][1] = fac*((20.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         + 2.0*(p->YN[JP1]-p->YN[JM1])*(p->YN[JP]-p->YN[JM1])
                         + (p->YN[JP1]-p->YN[JM2])*(p->YN[JP1]+p->YN[JP]-2.0*p->YN[JM1]))
    
                      / ((p->YN[JP1]-p->YN[JM1])*(p->YN[JP]-p->YN[JM2])));
                  
                  
    isfy[JP][0][2][2] = fac*((10.0*pow(p->YN[JP1]-p->YN[JP], 2.0) 
                         + (2.0*p->YN[JP1]-p->YN[JM2]-p->YN[JM1])
                            *(p->YN[JP1]+p->YN[JP]-p->YN[JM1]-p->YN[JM2]))

                      / pow(p->YN[JP1]-p->YN[JM1], 2.0));                 
                      
    // jmax ---------
    
    // is1
    fac = 4.0*pow((p->YN[JP2]-p->YN[JP1])/(p->YN[JP4]-p->YN[JP1]), 2.0);
    
    isfy[JP][0][3][0] = fac*((10.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         + (p->YN[JP3]-p->YN[JP1])*(p->YN[JP3]-p->YN[JP2]))
    
                      / pow(p->YN[JP4]-p->YN[JP2], 2.0));
                    
                    
    isfy[JP][0][3][1] = fac*((20.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         + 2.0*(p->YN[JP3]-p->YN[JP1])*(p->YN[JP3]-p->YN[JP2])
                         + (p->YN[JP4]-p->YN[JP1])*(2.0*p->YN[JP3]-p->YN[JP2]-p->YN[JP1]))
    
                      / ((p->YN[JP4]-p->YN[JP2])*(p->YN[JP3]-p->YN[JP1])));
                  
                  
    isfy[JP][0][3][2] = fac*((10.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         + (p->YN[JP4]+p->YN[JP3]-2.0*p->YN[JP1])
                            *(p->YN[JP4]+p->YN[JP3]-p->YN[JP2]-p->YN[JP1]))

                      / pow(p->YN[JP3]-p->YN[JP1], 2.0));
                      
    // is2
    fac = 4.0*pow((p->YN[JP2]-p->YN[JP1])/(p->YN[JP3]-p->YN[JP]), 2.0);
    
    isfy[JP][0][4][0] = fac*((10.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         + (p->YN[JP3]-p->YN[JP1])*(p->YN[JP3]-p->YN[JP2]))
    
                      / pow(p->YN[JP2]-p->YN[JP], 2.0));
                    
                    
    isfy[JP][0][4][1] = fac*((20.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         - (p->YN[JP3]-p->YN[JP2])*(p->YN[JP1]-p->YN[JP])
                         - (p->YN[JP3]-p->YN[JP1])*(p->YN[JP2]-p->YN[JP]))
    
                      / ((p->YN[JP3]-p->YN[JP1])*(p->YN[JP2]-p->YN[JP])));
                  
                  
    isfy[JP][0][4][2] = fac*((10.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         + (p->YN[JP1]-p->YN[JP])*(p->YN[JP2]-p->YN[JP]))

                      / pow(p->YN[JP3]-p->YN[JP1], 2.0));
                                     

    // is3
    fac = 4.0*pow((p->YN[JP2]-p->YN[JP1])/(p->YN[JP2]-p->YN[JM1]), 2.0);
    
    isfy[JP][0][5][0] = fac*((10.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         + (p->YN[JP2]-p->YN[JP])*(p->YN[JP1]-p->YN[JP]))
    
                      / pow(p->YN[JP1]-p->YN[JM1], 2.0));
                    
                    
    isfy[JP][0][5][1] = fac*((20.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         + 2.0*(p->YN[JP2]-p->YN[JP])*(p->YN[JP1]-p->YN[JP])
                         + (p->YN[JP2]-p->YN[JM1])*(p->YN[JP2]+p->YN[JP1]-2.0*p->YN[JP]))
    
                      / ((p->YN[JP2]-p->YN[JP])*(p->YN[JP1]-p->YN[JM1])));
                  
                  
    isfy[JP][0][5][2] = fac*((10.0*pow(p->YN[JP2]-p->YN[JP1], 2.0) 
                         + (2.0*p->YN[JP2]-p->YN[JM1]-p->YN[JP])
                            *(p->YN[JP2]+p->YN[JP1]-p->YN[JP]-p->YN[JM1]))

                      / pow(p->YN[JP2]-p->YN[JP], 2.0));
    }
    
// ZN
    KBLOOP
    {
    // kmin
    
    // is1
    fac = 4.0*pow((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP3]-p->ZN[KP]), 2.0);
    
    isfz[KP][0][0][0] = fac*((10.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         + (p->ZN[KP2]-p->ZN[KP])*(p->ZN[KP2]-p->ZN[KP1]))
    
                      / pow(p->ZN[KP3]-p->ZN[KP1], 2.0));
                    
                    
    isfz[KP][0][0][1] = fac*((20.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         + 2.0*(p->ZN[KP2]-p->ZN[KP])*(p->ZN[KP2]-p->ZN[KP1])
                         + (p->ZN[KP3]-p->ZN[KP])*(2.0*p->ZN[KP2]-p->ZN[KP1]-p->ZN[KP]))
    
                      / ((p->ZN[KP3]-p->ZN[KP1])*(p->ZN[KP2]-p->ZN[KP])));
                  
                  
    isfz[KP][0][0][2] = fac*((10.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         + (p->ZN[KP3]+p->ZN[KP2]-2.0*p->ZN[KP])
                            *(p->ZN[KP3]+p->ZN[KP2]-p->ZN[KP1]-p->ZN[KP]))

                      / pow(p->ZN[KP2]-p->ZN[KP], 2.0));
                      
    // is2
    fac = 4.0*pow((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP2]-p->ZN[KM1]), 2.0);
    
    isfz[KP][0][1][0] = fac*((10.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         + (p->ZN[KP2]-p->ZN[KP])*(p->ZN[KP2]-p->ZN[KP1]))
    
                      / pow(p->ZN[KP1]-p->ZN[KM1], 2.0));
                    
                    
    isfz[KP][0][1][1] = fac*((20.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         - (p->ZN[KP2]-p->ZN[KP1])*(p->ZN[KP]-p->ZN[KM1])
                         - (p->ZN[KP2]-p->ZN[KP])*(p->ZN[KP1]-p->ZN[KM1]))
    
                      / ((p->ZN[KP2]-p->ZN[KP])*(p->ZN[KP1]-p->ZN[KM1])));
                  
                  
    isfz[KP][0][1][2] = fac*((10.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         + (p->ZN[KP]-p->ZN[KM1])*(p->ZN[KP1]-p->ZN[KM1]))

                      / pow(p->ZN[KP2]-p->ZN[KP], 2.0));
                      
                      
     // is3
    fac = 4.0*pow((p->ZN[KP1]-p->ZN[KP])/(p->ZN[KP1]-p->ZN[KM2]), 2.0);
    
    isfz[KP][0][2][0] = fac*((10.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         + (p->ZN[KP1]-p->ZN[KM1])*(p->ZN[KP]-p->ZN[KM1]))
    
                      / pow(p->ZN[KP]-p->ZN[KM2], 2.0));
                    
                    
    isfz[KP][0][2][1] = fac*((20.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         + 2.0*(p->ZN[KP1]-p->ZN[KM1])*(p->ZN[KP]-p->ZN[KM1])
                         + (p->ZN[KP1]-p->ZN[KM2])*(p->ZN[KP1]+p->ZN[KP]-2.0*p->ZN[KM1]))
    
                      / ((p->ZN[KP1]-p->ZN[KM1])*(p->ZN[KP]-p->ZN[KM2])));
                  
                  
    isfz[KP][0][2][2] = fac*((10.0*pow(p->ZN[KP1]-p->ZN[KP], 2.0) 
                         + (2.0*p->ZN[KP1]-p->ZN[KM2]-p->ZN[KM1])
                            *(p->ZN[KP1]+p->ZN[KP]-p->ZN[KM1]-p->ZN[KM2]))

                      / pow(p->ZN[KP1]-p->ZN[KM1], 2.0));
                                       
                      
    // kmax ---------
    
    // is1
    fac = 4.0*pow((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP4]-p->ZN[KP1]), 2.0);
    
    isfz[KP][0][3][0] = fac*((10.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         + (p->ZN[KP3]-p->ZN[KP1])*(p->ZN[KP3]-p->ZN[KP2]))
    
                      / pow(p->ZN[KP4]-p->ZN[KP2], 2.0));
                    
                    
    isfz[KP][0][3][1] = fac*((20.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         + 2.0*(p->ZN[KP3]-p->ZN[KP1])*(p->ZN[KP3]-p->ZN[KP2])
                         + (p->ZN[KP4]-p->ZN[KP1])*(2.0*p->ZN[KP3]-p->ZN[KP2]-p->ZN[KP1]))
    
                      / ((p->ZN[KP4]-p->ZN[KP2])*(p->ZN[KP3]-p->ZN[KP1])));
                  
                  
    isfz[KP][0][3][2] = fac*((10.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         + (p->ZN[KP4]+p->ZN[KP3]-2.0*p->ZN[KP1])
                            *(p->ZN[KP4]+p->ZN[KP3]-p->ZN[KP2]-p->ZN[KP1]))

                      / pow(p->ZN[KP3]-p->ZN[KP1], 2.0));
                                       
                    
    // is2
    fac = 4.0*pow((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP3]-p->ZN[KP]), 2.0);
    
    isfz[KP][0][4][0] = fac*((10.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         + (p->ZN[KP3]-p->ZN[KP1])*(p->ZN[KP3]-p->ZN[KP2]))
    
                      / pow(p->ZN[KP2]-p->ZN[KP], 2.0));
                    
                    
    isfz[KP][0][4][1] = fac*((20.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         - (p->ZN[KP3]-p->ZN[KP2])*(p->ZN[KP1]-p->ZN[KP])
                         - (p->ZN[KP3]-p->ZN[KP1])*(p->ZN[KP2]-p->ZN[KP]))
    
                      / ((p->ZN[KP3]-p->ZN[KP1])*(p->ZN[KP2]-p->ZN[KP])));
                  
                  
    isfz[KP][0][4][2] = fac*((10.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         + (p->ZN[KP1]-p->ZN[KP])*(p->ZN[KP2]-p->ZN[KP]))

                      / pow(p->ZN[KP3]-p->ZN[KP1], 2.0));
                      
    // is3
    fac = 4.0*pow((p->ZN[KP2]-p->ZN[KP1])/(p->ZN[KP2]-p->ZN[KM1]), 2.0);
    
    isfz[KP][0][5][0] = fac*((10.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         + (p->ZN[KP2]-p->ZN[KP])*(p->ZN[KP1]-p->ZN[KP]))
    
                      / pow(p->ZN[KP1]-p->ZN[KM1], 2.0));
                    
                    
    isfz[KP][0][5][1] = fac*((20.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         + 2.0*(p->ZN[KP2]-p->ZN[KP])*(p->ZN[KP1]-p->ZN[KP])
                         + (p->ZN[KP2]-p->ZN[KM1])*(p->ZN[KP2]+p->ZN[KP1]-2.0*p->ZN[KP]))
    
                      / ((p->ZN[KP2]-p->ZN[KP])*(p->ZN[KP1]-p->ZN[KM1])));
                  
                  
    isfz[KP][0][5][2] = fac*((10.0*pow(p->ZN[KP2]-p->ZN[KP1], 2.0) 
                         + (2.0*p->ZN[KP2]-p->ZN[KM1]-p->ZN[KP])
                            *(p->ZN[KP2]+p->ZN[KP1]-p->ZN[KP]-p->ZN[KM1]))

                      / pow(p->ZN[KP2]-p->ZN[KP], 2.0));
                      
    /*cout<<" is1: "<<isfz[KP][0][3][0]<<" "<<isfz[KP][0][3][1]<<" "<<isfz[KP][0][3][2]<<endl;
    cout<<" is2: "<<isfz[KP][0][4][0]<<" "<<isfz[KP][0][4][1]<<" "<<isfz[KP][0][4][2]<<endl;
    cout<<" is3: "<<isfz[KP][0][5][0]<<" "<<isfz[KP][0][5][1]<<" "<<isfz[KP][0][5][2]<<endl;*/
    
    }
    
// ---------------------------------------------------------------------

// XP
    IBLOOP
    {
    // imin
    
    // is1
    fac = 4.0*pow((p->XP[IP1]-p->XP[IP])/(p->XP[IP3]-p->XP[IP]), 2.0);
    
    isfx[IP][1][0][0] = fac*((10.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         + (p->XP[IP2]-p->XP[IP])*(p->XP[IP2]-p->XP[IP1]))
    
                      / pow(p->XP[IP3]-p->XP[IP1], 2.0));
                    
                    
    isfx[IP][1][0][1] = fac*((20.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         + 2.0*(p->XP[IP2]-p->XP[IP])*(p->XP[IP2]-p->XP[IP1])
                         + (p->XP[IP3]-p->XP[IP])*(2.0*p->XP[IP2]-p->XP[IP1]-p->XP[IP]))
    
                      / ((p->XP[IP3]-p->XP[IP1])*(p->XP[IP2]-p->XP[IP])));
                  
                  
    isfx[IP][1][0][2] = fac*((10.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         + (p->XP[IP3]+p->XP[IP2]-2.0*p->XP[IP])
                            *(p->XP[IP3]+p->XP[IP2]-p->XP[IP1]-p->XP[IP]))

                      / pow(p->XP[IP2]-p->XP[IP], 2.0));
                      
                      
    // is2
    fac = 4.0*pow((p->XP[IP1]-p->XP[IP])/(p->XP[IP2]-p->XP[IM1]), 2.0);
    
    isfx[IP][1][1][0] = fac*((10.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         + (p->XP[IP2]-p->XP[IP])*(p->XP[IP2]-p->XP[IP1]))
    
                      / pow(p->XP[IP1]-p->XP[IM1], 2.0));
                    
                    
    isfx[IP][1][1][1] = fac*((20.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         - (p->XP[IP2]-p->XP[IP1])*(p->XP[IP]-p->XP[IM1])
                         - (p->XP[IP2]-p->XP[IP])*(p->XP[IP1]-p->XP[IM1]))
    
                      / ((p->XP[IP2]-p->XP[IP])*(p->XP[IP1]-p->XP[IM1])));
                  
                  
    isfx[IP][1][1][2] = fac*((10.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         + (p->XP[IP]-p->XP[IM1])*(p->XP[IP1]-p->XP[IM1]))

                      / pow(p->XP[IP2]-p->XP[IP], 2.0));
       
    // is3
    fac = 4.0*pow((p->XP[IP1]-p->XP[IP])/(p->XP[IP1]-p->XP[IM2]), 2.0);
    
    isfx[IP][1][2][0] = fac*((10.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         + (p->XP[IP1]-p->XP[IM1])*(p->XP[IP]-p->XP[IM1]))
    
                      / pow(p->XP[IP]-p->XP[IM2], 2.0));
                    
                    
    isfx[IP][1][2][1] = fac*((20.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         + 2.0*(p->XP[IP1]-p->XP[IM1])*(p->XP[IP]-p->XP[IM1])
                         + (p->XP[IP1]-p->XP[IM2])*(p->XP[IP1]+p->XP[IP]-2.0*p->XP[IM1]))
    
                      / ((p->XP[IP1]-p->XP[IM1])*(p->XP[IP]-p->XP[IM2])));
                  
                  
    isfx[IP][1][2][2] = fac*((10.0*pow(p->XP[IP1]-p->XP[IP], 2.0) 
                         + (2.0*p->XP[IP1]-p->XP[IM2]-p->XP[IM1])
                            *(p->XP[IP1]+p->XP[IP]-p->XP[IM1]-p->XP[IM2]))

                      / pow(p->XP[IP1]-p->XP[IM1], 2.0));
                                        
          
    // imax ---------
    // is1
    fac = 4.0*pow((p->XP[IP2]-p->XP[IP1])/(p->XP[IP4]-p->XP[IP1]), 2.0);
    
    isfx[IP][1][3][0] = fac*((10.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         + (p->XP[IP3]-p->XP[IP1])*(p->XP[IP3]-p->XP[IP2]))
    
                      / pow(p->XP[IP4]-p->XP[IP2], 2.0));
                    
                    
    isfx[IP][1][3][1] = fac*((20.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         + 2.0*(p->XP[IP3]-p->XP[IP1])*(p->XP[IP3]-p->XP[IP2])
                         + (p->XP[IP4]-p->XP[IP1])*(2.0*p->XP[IP3]-p->XP[IP2]-p->XP[IP1]))
    
                      / ((p->XP[IP4]-p->XP[IP2])*(p->XP[IP3]-p->XP[IP1])));
                  
                  
    isfx[IP][1][3][2] = fac*((10.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         + (p->XP[IP4]+p->XP[IP3]-2.0*p->XP[IP1])
                            *(p->XP[IP4]+p->XP[IP3]-p->XP[IP2]-p->XP[IP1]))

                      / pow(p->XP[IP3]-p->XP[IP1], 2.0));
     //cout<<" is1: "<<isfx[IP][0][0][0]<<" "<<isfx[IP][0][0][1]<<" "<<isfx[IP][0][0][2]<<endl;
                 
    // is2
    fac = 4.0*pow((p->XP[IP2]-p->XP[IP1])/(p->XP[IP3]-p->XP[IP]), 2.0);
    
    isfx[IP][1][4][0] = fac*((10.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         + (p->XP[IP3]-p->XP[IP1])*(p->XP[IP3]-p->XP[IP2]))
    
                      / pow(p->XP[IP2]-p->XP[IP], 2.0));
                    
                    
    isfx[IP][1][4][1] = fac*((20.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         - (p->XP[IP3]-p->XP[IP2])*(p->XP[IP1]-p->XP[IP])
                         - (p->XP[IP3]-p->XP[IP1])*(p->XP[IP2]-p->XP[IP]))
    
                      / ((p->XP[IP3]-p->XP[IP1])*(p->XP[IP2]-p->XP[IP])));
                  
                  
    isfx[IP][1][4][2] = fac*((10.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         + (p->XP[IP1]-p->XP[IP])*(p->XP[IP2]-p->XP[IP]))

                      / pow(p->XP[IP3]-p->XP[IP1], 2.0));
      
    //cout<<" is2: "<<isfx[IP][0][1][0]<<" "<<isfx[IP][0][1][1]<<" "<<isfx[IP][0][1][2]<<endl;     
              
    // is3
    fac = 4.0*pow((p->XP[IP2]-p->XP[IP1])/(p->XP[IP2]-p->XP[IM1]), 2.0);
    
    isfx[IP][1][5][0] = fac*((10.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         + (p->XP[IP2]-p->XP[IP])*(p->XP[IP1]-p->XP[IP]))
    
                      / pow(p->XP[IP1]-p->XP[IM1], 2.0));
                    
                    
    isfx[IP][1][5][1] = fac*((20.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         + 2.0*(p->XP[IP2]-p->XP[IP])*(p->XP[IP1]-p->XP[IP])
                         + (p->XP[IP2]-p->XP[IM1])*(p->XP[IP2]+p->XP[IP1]-2.0*p->XP[IP]))
    
                      / ((p->XP[IP2]-p->XP[IP])*(p->XP[IP1]-p->XP[IM1])));
                  
                  
    isfx[IP][1][5][2] = fac*((10.0*pow(p->XP[IP2]-p->XP[IP1], 2.0) 
                         + (2.0*p->XP[IP2]-p->XP[IM1]-p->XP[IP])
                            *(p->XP[IP2]+p->XP[IP1]-p->XP[IP]-p->XP[IM1]))

                      / pow(p->XP[IP2]-p->XP[IP], 2.0));
                      
    //cout<<" is3: "<<isfx[IP][0][2][0]<<" "<<isfx[IP][0][2][1]<<" "<<isfx[IP][0][2][2]<<endl;   
    }
    
    
// YP
    JBLOOP
    {
    // jmin
    // is1
    fac = 4.0*pow((p->YP[JP1]-p->YP[JP])/(p->YP[JP3]-p->YP[JP]), 2.0);
    
    isfy[JP][1][0][0] = fac*((10.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         + (p->YP[JP2]-p->YP[JP])*(p->YP[JP2]-p->YP[JP1]))
    
                      / pow(p->YP[JP3]-p->YP[JP1], 2.0));
                    
                    
    isfy[JP][1][0][1] = fac*((20.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         + 2.0*(p->YP[JP2]-p->YP[JP])*(p->YP[JP2]-p->YP[JP1])
                         + (p->YP[JP3]-p->YP[JP])*(2.0*p->YP[JP2]-p->YP[JP1]-p->YP[JP]))
    
                      / ((p->YP[JP3]-p->YP[JP1])*(p->YP[JP2]-p->YP[JP])));
                  
                  
    isfy[JP][1][0][2] = fac*((10.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         + (p->YP[JP3]+p->YP[JP2]-2.0*p->YP[JP])
                            *(p->YP[JP3]+p->YP[JP2]-p->YP[JP1]-p->YP[JP]))

                      / pow(p->YP[JP2]-p->YP[JP], 2.0));
                      
    // is2
    fac = 4.0*pow((p->YP[JP1]-p->YP[JP])/(p->YP[JP2]-p->YP[JM1]), 2.0);
    
    isfy[JP][1][1][0] = fac*((10.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         + (p->YP[JP2]-p->YP[JP])*(p->YP[JP2]-p->YP[JP1]))
    
                      / pow(p->YP[JP1]-p->YP[JM1], 2.0));
                    
                    
    isfy[JP][1][1][1] = fac*((20.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         - (p->YP[JP2]-p->YP[JP1])*(p->YP[JP]-p->YP[JM1])
                         - (p->YP[JP2]-p->YP[JP])*(p->YP[JP1]-p->YP[JM1]))
    
                      / ((p->YP[JP2]-p->YP[JP])*(p->YP[JP1]-p->YP[JM1])));
                  
                  
    isfy[JP][1][1][2] = fac*((10.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         + (p->YP[JP]-p->YP[JM1])*(p->YP[JP1]-p->YP[JM1]))

                      / pow(p->YP[JP2]-p->YP[JP], 2.0));
                      
    // is3
    fac = 4.0*pow((p->YP[JP1]-p->YP[JP])/(p->YP[JP1]-p->YP[JM2]), 2.0);
    
    isfy[JP][1][2][0] = fac*((10.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         + (p->YP[JP1]-p->YP[JM1])*(p->YP[JP]-p->YP[JM1]))
    
                      / pow(p->YP[JP]-p->YP[JM2], 2.0));
                    
                    
    isfy[JP][1][2][1] = fac*((20.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         + 2.0*(p->YP[JP1]-p->YP[JM1])*(p->YP[JP]-p->YP[JM1])
                         + (p->YP[JP1]-p->YP[JM2])*(p->YP[JP1]+p->YP[JP]-2.0*p->YP[JM1]))
    
                      / ((p->YP[JP1]-p->YP[JM1])*(p->YP[JP]-p->YP[JM2])));
                  
                  
    isfy[JP][1][2][2] = fac*((10.0*pow(p->YP[JP1]-p->YP[JP], 2.0) 
                         + (2.0*p->YP[JP1]-p->YP[JM2]-p->YP[JM1])
                            *(p->YP[JP1]+p->YP[JP]-p->YP[JM1]-p->YP[JM2]))

                      / pow(p->YP[JP1]-p->YP[JM1], 2.0));                  
    
                  
    // jmax ---------
    
    // is1
    fac = 4.0*pow((p->YP[JP2]-p->YP[JP1])/(p->YP[JP4]-p->YP[JP1]), 2.0);
    
    isfy[JP][1][3][0] = fac*((10.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         + (p->YP[JP3]-p->YP[JP1])*(p->YP[JP3]-p->YP[JP2]))
    
                      / pow(p->YP[JP4]-p->YP[JP2], 2.0));
                    
                    
    isfy[JP][1][3][1] = fac*((20.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         + 2.0*(p->YP[JP3]-p->YP[JP1])*(p->YP[JP3]-p->YP[JP2])
                         + (p->YP[JP4]-p->YP[JP1])*(2.0*p->YP[JP3]-p->YP[JP2]-p->YP[JP1]))
    
                      / ((p->YP[JP4]-p->YP[JP2])*(p->YP[JP3]-p->YP[JP1])));
                  
                  
    isfy[JP][1][3][2] = fac*((10.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         + (p->YP[JP4]+p->YP[JP3]-2.0*p->YP[JP1])
                            *(p->YP[JP4]+p->YP[JP3]-p->YP[JP2]-p->YP[JP1]))

                      / pow(p->YP[JP3]-p->YP[JP1], 2.0));
                      
    // is2
    fac = 4.0*pow((p->YP[JP2]-p->YP[JP1])/(p->YP[JP3]-p->YP[JP]), 2.0);
    
    isfy[JP][1][4][0] = fac*((10.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         + (p->YP[JP3]-p->YP[JP1])*(p->YP[JP3]-p->YP[JP2]))
    
                      / pow(p->YP[JP2]-p->YP[JP], 2.0));
                    
                    
    isfy[JP][1][4][1] = fac*((20.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         - (p->YP[JP3]-p->YP[JP2])*(p->YP[JP1]-p->YP[JP])
                         - (p->YP[JP3]-p->YP[JP1])*(p->YP[JP2]-p->YP[JP]))
    
                      / ((p->YP[JP3]-p->YP[JP1])*(p->YP[JP2]-p->YP[JP])));
                  
                  
    isfy[JP][1][4][2] = fac*((10.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         + (p->YP[JP1]-p->YP[JP])*(p->YP[JP2]-p->YP[JP]))

                      / pow(p->YP[JP3]-p->YP[JP1], 2.0));
                      

    // is3
    fac = 4.0*pow((p->YP[JP2]-p->YP[JP1])/(p->YP[JP2]-p->YP[JM1]), 2.0);
    
    isfy[JP][1][5][0] = fac*((10.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         + (p->YP[JP2]-p->YP[JP])*(p->YP[JP1]-p->YP[JP]))
    
                      / pow(p->YP[JP1]-p->YP[JM1], 2.0));
                    
                    
    isfy[JP][1][5][1] = fac*((20.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         + 2.0*(p->YP[JP2]-p->YP[JP])*(p->YP[JP1]-p->YP[JP])
                         + (p->YP[JP2]-p->YP[JM1])*(p->YP[JP2]+p->YP[JP1]-2.0*p->YP[JP]))
    
                      / ((p->YP[JP2]-p->YP[JP])*(p->YP[JP1]-p->YP[JM1])));
                  
                  
    isfy[JP][1][5][2] = fac*((10.0*pow(p->YP[JP2]-p->YP[JP1], 2.0) 
                         + (2.0*p->YP[JP2]-p->YP[JM1]-p->YP[JP])
                            *(p->YP[JP2]+p->YP[JP1]-p->YP[JP]-p->YP[JM1]))

                      / pow(p->YP[JP2]-p->YP[JP], 2.0));
    }
    
// ZP
    KBLOOP
    {
    // kmin
    
    // is1
    fac = 4.0*pow((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP3]-p->ZP[KP]), 2.0);
    
    isfz[KP][1][0][0] = fac*((10.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         + (p->ZP[KP2]-p->ZP[KP])*(p->ZP[KP2]-p->ZP[KP1]))
    
                      / pow(p->ZP[KP3]-p->ZP[KP1], 2.0));
                    
                    
    isfz[KP][1][0][1] = fac*((20.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         + 2.0*(p->ZP[KP2]-p->ZP[KP])*(p->ZP[KP2]-p->ZP[KP1])
                         + (p->ZP[KP3]-p->ZP[KP])*(2.0*p->ZP[KP2]-p->ZP[KP1]-p->ZP[KP]))
    
                      / ((p->ZP[KP3]-p->ZP[KP1])*(p->ZP[KP2]-p->ZP[KP])));
                  
                  
    isfz[KP][1][0][2] = fac*((10.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         + (p->ZP[KP3]+p->ZP[KP2]-2.0*p->ZP[KP])
                            *(p->ZP[KP3]+p->ZP[KP2]-p->ZP[KP1]-p->ZP[KP]))

                      / pow(p->ZP[KP2]-p->ZP[KP], 2.0));
                      
    // is2
    fac = 4.0*pow((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP2]-p->ZP[KM1]), 2.0);
    
    isfz[KP][1][1][0] = fac*((10.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         + (p->ZP[KP2]-p->ZP[KP])*(p->ZP[KP2]-p->ZP[KP1]))
    
                      / pow(p->ZP[KP1]-p->ZP[KM1], 2.0));
                    
                    
    isfz[KP][1][1][1] = fac*((20.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         - (p->ZP[KP2]-p->ZP[KP1])*(p->ZP[KP]-p->ZP[KM1])
                         - (p->ZP[KP2]-p->ZP[KP])*(p->ZP[KP1]-p->ZP[KM1]))
    
                      / ((p->ZP[KP2]-p->ZP[KP])*(p->ZP[KP1]-p->ZP[KM1])));
                  
                  
    isfz[KP][1][1][2] = fac*((10.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         + (p->ZP[KP]-p->ZP[KM1])*(p->ZP[KP1]-p->ZP[KM1]))

                      / pow(p->ZP[KP2]-p->ZP[KP], 2.0));
                      
    // is3
    fac = 4.0*pow((p->ZP[KP1]-p->ZP[KP])/(p->ZP[KP1]-p->ZP[KM2]), 2.0);
    
    isfz[KP][1][2][0] = fac*((10.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         + (p->ZP[KP1]-p->ZP[KM1])*(p->ZP[KP]-p->ZP[KM1]))
    
                      / pow(p->ZP[KP]-p->ZP[KM2], 2.0));
                    
                    
    isfz[KP][1][2][1] = fac*((20.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         + 2.0*(p->ZP[KP1]-p->ZP[KM1])*(p->ZP[KP]-p->ZP[KM1])
                         + (p->ZP[KP1]-p->ZP[KM2])*(p->ZP[KP1]+p->ZP[KP]-2.0*p->ZP[KM1]))
    
                      / ((p->ZP[KP1]-p->ZP[KM1])*(p->ZP[KP]-p->ZP[KM2])));
                  
                  
    isfz[KP][1][2][2] = fac*((10.0*pow(p->ZP[KP1]-p->ZP[KP], 2.0) 
                         + (2.0*p->ZP[KP1]-p->ZP[KM2]-p->ZP[KM1])
                            *(p->ZP[KP1]+p->ZP[KP]-p->ZP[KM1]-p->ZP[KM2]))

                      / pow(p->ZP[KP1]-p->ZP[KM1], 2.0));                  
                      
    // kmax ---------
    
    // is1
    fac = 4.0*pow((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP4]-p->ZP[KP1]), 2.0);
    
    isfz[KP][1][3][0] = fac*((10.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         + (p->ZP[KP3]-p->ZP[KP1])*(p->ZP[KP3]-p->ZP[KP2]))
    
                      / pow(p->ZP[KP4]-p->ZP[KP2], 2.0));
                    
                    
    isfz[KP][1][3][1] = fac*((20.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         + 2.0*(p->ZP[KP3]-p->ZP[KP1])*(p->ZP[KP3]-p->ZP[KP2])
                         + (p->ZP[KP4]-p->ZP[KP1])*(2.0*p->ZP[KP3]-p->ZP[KP2]-p->ZP[KP1]))
    
                      / ((p->ZP[KP4]-p->ZP[KP2])*(p->ZP[KP3]-p->ZP[KP1])));
                  
                  
    isfz[KP][1][3][2] = fac*((10.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         + (p->ZP[KP4]+p->ZP[KP3]-2.0*p->ZP[KP1])
                            *(p->ZP[KP4]+p->ZP[KP3]-p->ZP[KP2]-p->ZP[KP1]))

                      / pow(p->ZP[KP3]-p->ZP[KP1], 2.0));
    
                 
    // is2
    fac = 4.0*pow((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP3]-p->ZP[KP]), 2.0);
    
    isfz[KP][1][4][0] = fac*((10.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         + (p->ZP[KP3]-p->ZP[KP1])*(p->ZP[KP3]-p->ZP[KP2]))
    
                      / pow(p->ZP[KP2]-p->ZP[KP], 2.0));
                    
                    
    isfz[KP][1][4][1] = fac*((20.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         - (p->ZP[KP3]-p->ZP[KP2])*(p->ZP[KP1]-p->ZP[KP])
                         - (p->ZP[KP3]-p->ZP[KP1])*(p->ZP[KP2]-p->ZP[KP]))
    
                      / ((p->ZP[KP3]-p->ZP[KP1])*(p->ZP[KP2]-p->ZP[KP])));
                  
                  
    isfz[KP][1][4][2] = fac*((10.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         + (p->ZP[KP1]-p->ZP[KP])*(p->ZP[KP2]-p->ZP[KP]))

                      / pow(p->ZP[KP3]-p->ZP[KP1], 2.0));
                      
     // is3
    fac = 4.0*pow((p->ZP[KP2]-p->ZP[KP1])/(p->ZP[KP2]-p->ZP[KM1]), 2.0);
    
    isfz[KP][1][5][0] = fac*((10.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         + (p->ZP[KP2]-p->ZP[KP])*(p->ZP[KP1]-p->ZP[KP]))
    
                      / pow(p->ZP[KP1]-p->ZP[KM1], 2.0));
                    
                    
    isfz[KP][1][5][1] = fac*((20.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         + 2.0*(p->ZP[KP2]-p->ZP[KP])*(p->ZP[KP1]-p->ZP[KP])
                         + (p->ZP[KP2]-p->ZP[KM1])*(p->ZP[KP2]+p->ZP[KP1]-2.0*p->ZP[KP]))
    
                      / ((p->ZP[KP2]-p->ZP[KP])*(p->ZP[KP1]-p->ZP[KM1])));
                  
                  
    isfz[KP][1][5][2] = fac*((10.0*pow(p->ZP[KP2]-p->ZP[KP1], 2.0) 
                         + (2.0*p->ZP[KP2]-p->ZP[KM1]-p->ZP[KP])
                            *(p->ZP[KP2]+p->ZP[KP1]-p->ZP[KP]-p->ZP[KM1]))

                      / pow(p->ZP[KP2]-p->ZP[KP], 2.0));                 
    }
}
