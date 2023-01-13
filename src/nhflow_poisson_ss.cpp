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

#include"nhflow_pjm_ss.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"hypre_struct_fnpf.h"
#include"hypre_sstruct_fnpf.h"

void nhflow_pjm_ss::poisson2D(lexer* p, fdm_nhf *d, field &f)
{
    double sigxyz2,ab,denom;
    double fbxm,fbxp,fbym,fbyp;
    p->poissoniter=0;
    p->poissontime=0.0;
    
/*
 * {{0,0}, {-1,0}, {1,0},  {0,-1}, {0,1}, {-1,-1},{-1,1},{1,-1},{1,1}};
p  0
s  1
n  2
b  3
t  4 -
sb 5
st 6 -
nb 7
nt 8 -
*/

	n=0;
    KJILOOP
	{
        if(p->flag4[IJK]>0 && p->wet[IJ]==1)
        {
        sigxyz2 = pow(0.5*(p->sigx[FIJK]+p->sigx[FIJKp1]),2.0) + pow(p->sigz[IJ],2.0);
        
        M[n*9]  =         (CPORNH*PORVALNH)/(p->W1*p->DXP[IP]*p->DXN[IP])
                        + (CPORNH*PORVALNH)/(p->W1*p->DXP[IM1]*p->DXN[IP])
                        
                        + (sigxyz2*CPORNH*PORVALNH)/(p->W1*p->DZP[KP]*p->DZN[KP])
                        + (sigxyz2*CPORNH*PORVALNH)/(p->W1*p->DZP[KM1]*p->DZN[KP]);

        
        M[n*9+1] = -(CPORNH*PORVALNH)/(p->W1*p->DXP[IP]*p->DXN[IP]);
        
        M[n*9+2] = -(CPORNH*PORVALNH)/(p->W1*p->DXP[IM1]*p->DXN[IP]);
        
        M[n*9+3] = -(sigxyz2*CPORNH*PORVALNH)/(p->W1*p->DZP[KP]*p->DZN[KP])     
                        + CPORNH*PORVALNH*0.5*(p->sigxx[FIJK]+p->sigxx[FIJKp1])/(p->W1*(p->DZN[KP]+p->DZN[KM1]));
                        
        M[n*9+4] = -(sigxyz2*CPORNH*PORVALNH)/(p->W1*p->DZP[KM1]*p->DZN[KP]) 
                        - CPORNH*PORVALNH*0.5*(p->sigxx[FIJK]+p->sigxx[FIJKp1])/(p->W1*(p->DZN[KP]+p->DZN[KM1]));
        
      
        M[n*9+5]  = -CPORNH*PORVALNH*(p->sigx[FIJK]+p->sigx[FIJKp1])/(p->W1*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
        M[n*9+6]  =  CPORNH*PORVALNH*(p->sigx[FIJK]+p->sigx[FIJKp1])/(p->W1*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
        M[n*9+7]  =  CPORNH*PORVALNH*(p->sigx[FIJK]+p->sigx[FIJKp1])/(p->W1*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
        M[n*9+8]  = -CPORNH*PORVALNH*(p->sigx[FIJK]+p->sigx[FIJKp1])/(p->W1*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1])); 
        }
        
        if(p->flag4[IJK]<0 || p->wet[IJ]==0)
        {
        M[n*9]   = 1.0;
        M[n*9+1] = 0.0;
        M[n*9+2] = 0.0;
        M[n*9+3] = 0.0;
        M[n*9+4] = 0.0;
        M[n*9+5] = 0.0;
        M[n*9+6] = 0.0;
        M[n*9+7] = 0.0;
        M[n*9+8] = 0.0;
        }
	++n;
	}
    
    n=0;
	KJILOOP
	{
        if(p->flag4[IJK]>0 && p->wet[IJ]==1)
        {
            // south
            if((p->flag4[Im1JK]<0 || p->wet[Im1J]==0))
            {
            rhs[n] -= M[n*9+1]*f(i-1,j,k);
            M[n*9+1] = 0.0;       
            }
            
            // north
            if((p->flag4[Ip1JK]<0 || p->wet[Ip1J]==0))
            {
            rhs[n] -= M[n*9+2]*f(i+1,j,k);
            M[n*9+2] = 0.0;
            }

            // KFSFBC
            if(p->flag4[IJKp1]<0)
            {
                if(p->D37==1)
                {
                rhs[n] -= M[n*9+4]*f(i,j,k+1);
                M[n*9+4] = 0.0;
                }
            }
   
        // diagonal entries
            // st
                // fsfbc
            if(p->flag4[Im1JKp1]<0 && p->flag4[IJKp1]<0) // fsfbc
            {
                if(p->D37==1)
                {
                rhs[n] -= M[n*9+6]*f(i-1,j,k+1);
                M[n*9+6] = 0.0;
                }
            }
                // wall
            if((p->flag4[Im1JKp1]<0 && p->flag4[IJKp1]>0)) //
            {
                if(p->D37==1)
                {
                rhs[n] -= M[n*9+6]*f(i-1,j,k+1);
                M[n*9+6] = 0.0;
                }
            }
            
            // nt
                // fsfbc
            if(p->flag4[Ip1JKp1]<0 && p->flag4[IJKp1]<0) 
            {
                
                if(p->D37==1)
                {
                rhs[n] -= M[n*9+8]*f(i+1,j,k+1);
                M[n*9+8] = 0.0; 
                }
                
            }
            
                // wall
            if(p->flag4[Ip1JKp1]<0 && p->flag4[IJKp1]>0)
            {
                if(p->D37==1)
                {
                rhs[n] -= M[n*9+8]*f(i+1,j,k+1);
                M[n*9+8] = 0.0; 
                }

            }
            
            // sb 
                // wall
            if(((p->flag4[Im1JKm1]<0 && p->flag4[IJKm1]>0)|| p->wet[Im1J]==0))
            {
            rhs[n] -= M[n*9+5]*f(i-1,j,k-1);
            M[n*9+5] = 0.0;       
            }
        
            // nb 
                // wall
            if(((p->flag4[Ip1JKm1]<0 && p->flag4[IJKm1]>0)|| p->wet[Ip1J]==0))
            {
            rhs[n] -= M[n*9+7]*f(i+1,j,k-1);
            M[n*9+7] = 0.0;       
            }
                
            // sb KBEDBC
            if(p->flag4[Im1JKm1]<0 && p->flag4[IJKm1]<0)
            {
            rhs[n] -= M[n*9+5]*f(i-1,j,k-1);
            M[n*9+5] = 0.0;
            }
            
            // nb KBEDBC
            if(p->flag4[Ip1JKm1]<0 && p->flag4[IJKm1]<0)
            {
            rhs[n] -= M[n*9+7]*f(i+1,j,k-1);
            M[n*9+7] = 0.0;
            }
 
            // KBEDBC
            if(p->flag4[IJKm1]<0)
            {
            rhs[n] -= M[n*9+3]*f(i,j,k-1);
            M[n*9+3] = 0.0;

            }
        }
        
	++n;
	}
}

void nhflow_pjm_ss::poisson3D(lexer* p, fdm_nhf *d, field &f)
{
    
}

