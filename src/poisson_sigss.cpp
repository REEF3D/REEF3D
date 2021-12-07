/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
--------------------------------------------------------------------*/

#include"pjm_sigss.h"
#include"lexer.h"
#include"fdm.h"
#include"heat.h"
#include"concentration.h"
#include"density_f.h"
#include"density_comp.h"
#include"density_conc.h"
#include"density_heat.h"
#include"density_vof.h"
#include"hypre_struct_fnpf.h"
#include"hypre_sstruct_fnpf.h"

void pjm_sigss::poisson2D(lexer* p, fdm *a, field &f)
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
t  4
sb 5
st 6
nb 7
nt 8
*/

	n=0;
    KJILOOP
	{
        if(p->flag4[IJK]>0 && a->wet(i,j)==1)
        {
        sigxyz2 = pow(0.5*(p->sigx[FIJK]+p->sigx[FIJKp1]),2.0) + pow(0.5*(p->sigy[FIJK]+p->sigy[FIJKp1]),2.0) + pow(p->sigz[IJ],2.0);
        
        M[n*9]  =         (CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DXP[IP]*p->DXN[IP])
                        + (CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DXP[IM1]*p->DXN[IP])
                        
                        + (sigxyz2*CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DZP[KP]*p->DZN[KP])
                        + (sigxyz2*CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*p->DZP[KM1]*p->DZN[KP]);

        
        M[n*9+1] = -(CPOR1*PORVAL1)/(pd->roface(p,a,1,0,0)*p->DXP[IP]*p->DXN[IP]);
        M[n*9+2] = -(CPOR1m*PORVAL1m)/(pd->roface(p,a,-1,0,0)*p->DXP[IM1]*p->DXN[IP]);
        
        M[n*9+3] = -(sigxyz2*CPOR3*PORVAL3)/(pd->roface(p,a,0,0,1)*p->DZP[KP]*p->DZN[KP])     
                        + CPOR4*PORVAL4*0.5*(p->sigxx[FIJK]+p->sigxx[FIJKp1])/(a->ro(i,j,k)*(p->DZN[KP]+p->DZN[KM1]));
        M[n*9+4] = -(sigxyz2*CPOR3m*PORVAL3m)/(pd->roface(p,a,0,0,-1)*p->DZP[KM1]*p->DZN[KP]) 
                        - CPOR4*PORVAL4*0.5*(p->sigxx[FIJK]+p->sigxx[FIJKp1])/(a->ro(i,j,k)*(p->DZN[KP]+p->DZN[KM1]));
        
      
        M[n*9+5]  = -CPOR4*PORVAL4*(p->sigx[FIJK]+p->sigx[FIJKp1])/(a->ro(i,j,k)*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
        M[n*9+6]  =  CPOR4*PORVAL4*(p->sigx[FIJK]+p->sigx[FIJKp1])/(a->ro(i,j,k)*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
        M[n*9+7]  =  CPOR4*PORVAL4*(p->sigx[FIJK]+p->sigx[FIJKp1])/(a->ro(i,j,k)*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
        M[n*9+8]  = -CPOR4*PORVAL4*(p->sigx[FIJK]+p->sigx[FIJKp1])/(a->ro(i,j,k)*(p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1])); 
        }
        
        if(p->flag4[IJK]<0 || a->wet(i,j)==0)
        {
        M[n*9]  =  1.0;
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
        if(p->flag4[IJK]>0 && a->wet(i,j)==1)
        {
            // south
            if((p->flag4[Im1JK]<0 || a->wet(i-1,j)==0))
            {
            rhs[n] -= M[n*9+1]*f[Im1JK];
            M[n*9+1] = 0.0;       
            }
            
            // north
            if((p->flag4[Ip1JK]<0 || a->wet(i+1,j)==0))
            {
            rhs[n] -= M[n*9+2]*f[Ip1JK];
            M[n*9+2] = 0.0;
            }

            // top
            if(p->flag4[IJKp1]<0)
            {
            rhs[n] -= M[n*9+4]*f[IJKp1];
            M[n*9+4] = 0.0;
            }
   
        // diagonal entries
            // st
                // fsfbc
            if(p->flag4[Im1JKp1]<0 && p->flag4[IJKp1]<0) // fsfbc
            {
            rhs[n] -= M[n*9+6]*f[Im1JKp1];
            M[n*9+6] = 0.0;
            }
                // wall
            if((p->flag4[Im1JKp1]<0 && p->flag4[IJKp1]>0)) //
            {
            rhs[n] -= M[n*9+6]*f[Im1JKp1];
            M[n*9+6] = 0.0;  
            }
            
            // nt
                // fsfbc
            if(p->flag4[Ip1JKp1]<0 && p->flag4[IJKp1]<0) 
            {
            rhs[n] -= M[n*9+8]*f(i+1,j,k+1);
            M[n*9+8] = 0.0; 
            }
            
                // wall
            if(p->flag4[Ip1JKp1]<0 && p->flag4[IJKp1]>0)
            {
            rhs[n] -= M[n*9+8]*f(i+1,j,k+1);
            M[n*9+8] = 0.0; 
            }
            
            // sb 
                // wall
            if(((p->flag4[Im1JKm1]<0 && p->flag4[IJKm1]>0)|| a->wet(i-1,j)==0))
            {
            rhs[n] -= M[n*9+5]*f(i-1,j,k-1);
            M[n*9+5] = 0.0;       
            }
        
            // nb 
                // wall
            if(((p->flag4[Ip1JKm1]<0 && p->flag4[IJKm1]>0)|| a->wet(i+1,j)==0))
            {
            rhs[n] -= M[n*9+7]*f(i+1,j,k-1);
            M[n*9+7] = 0.0;       
            }
                
            // sb KBEDBC
            if(p->flag4[Im1JKm1]<0 && p->flag4[IJKm1]<0)
            { 
            /*ab = -2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]))*p->x_dir;
            
            denom = p->sigz[Im1J] + c->Bx(i-1,j)*p->sigx[FIm1JK];

                    if(c->wet(i+1,j)==1 && c->wet(i-1,j)==1)
                    {
                    M[n*9+2] +=  ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*9+1] += -ab*2.0*p->DZN[KP]*c->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                
                M[n*9+6] += ab;
                M[n*9+5] = 0.0;*/
                
            rhs[n] -= M[n*9+5]*f(i-1,j,k-1);
            M[n*9+5] = 0.0;
            }
            
            // nb KBEDBC
            if(p->flag4[Ip1JKm1]<0 && p->flag4[IJKm1]<0)
            {
            /*ab = 2.0*p->sigx[FIJK]/((p->DXN[IP]+p->DXN[IM1])*(p->DZN[KP]+p->DZN[KM1]));
            
            denom = p->sigz[IJ] + a->Bx(i,j)*p->sigx[FIJK];

                    if(a->wet(i+1,j)==1 && a->wet(i-1,j)==1)
                    {
                    M[n*9+2] +=  ab*2.0*p->DZN[KP]*a->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*9+1] += -ab*2.0*p->DZN[KP]*a->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                
                M[n*9+8] += ab;
            M[n*9+7] = 0.0;*/
            
            rhs[n] -= M[n*9+7]*f(i+1,j,k-1);
            M[n*9+7] = 0.0;
            }
 
            // KBEDBC
            if(p->flag4[IJKm1]<0)
            {
            /*sigxyz2 = pow(p->sigx[FIJK],2.0) + pow(p->sigy[FIJK],2.0) + pow(p->sigz[IJ],2.0);
            
            ab = -(sigxyz2/(p->DZP[KM1]*p->DZN[KM1]) - p->sigxx[FIJK]/(p->DZN[KP]+p->DZN[KM1]));
            
            denom = p->sigz[IJ] + a->Bx(i,j)*p->sigx[FIJK];

                    if(a->wet(i+1,j)==1 && a->wet(i-1,j)==1)
                    {
                    M[n*9+2] +=  ab*2.0*p->DZN[KP]*a->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    M[n*9+1] += -ab*2.0*p->DZN[KP]*a->Bx(i,j)/(denom*(p->DXP[IP] + p->DXP[IM1]));
                    }
                    
                M[n*9+4] += ab;
                M[n*9+3] = 0.0;*/
                
            
            rhs[n] -= M[n*9+3]*f(i,j,k-1);
            M[n*9+3] = 0.0;

            }
        }
        
	++n;
	}
    
    
    

}


void pjm_sigss::poisson3D(lexer* p, fdm *a, field &f)
{
    
}
