/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2018 Hans Bihs

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

#include"poisson_4th.h"
#include"lexer.h"
#include"fdm.h"

poisson_4th::poisson_4th(lexer * p) : density(p)
{
}

poisson_4th::~poisson_4th()
{
}

void poisson_4th::estart(lexer* p, fdm *a, field &f)
{	
    double X1,X2,X3,X4,X0;
    double Y1,Y2,Y3,Y4,Y0;
    double Z1,Z2,Z3,Z4,Z0;
   
	// see p. 1130-1132
    n=0;
    LOOP
	{
    X0 = -0.5*p->XN[IP2] + 13.0*p->XN[IP1] - 13.0*p->XN[IM1] + 0.5*p->XN[IM2];
    X1 = (-p->XN[IP3] + 27.0*p->XN[IP2] -27.0*p->XN[IP1] + p->XN[IP])*X0;
    X2 = (-p->XN[IP2] + 27.0*p->XN[IP1] -27.0*p->XN[IP] + p->XN[IM1])*X0;
    X3 = (-p->XN[IP1] + 27.0*p->XN[IP] -27.0*p->XN[IM1] + p->XN[IM2])*X0;
    X4 = (-p->XN[IP] + 27.0*p->XN[IM1] -27.0*p->XN[IM2] + p->XN[IM3])*X0;
    
    Y0 = -0.5*p->YN[JP2] + 13.0*p->YN[JP1] - 13.0*p->YN[JM1] + 0.5*p->YN[JM2];
    Y1 = (-p->YN[JP3] + 27.0*p->YN[JP2] -27.0*p->YN[JP1] + p->YN[JP])*Y0;
    Y2 = (-p->YN[JP2] + 27.0*p->YN[JP1] -27.0*p->YN[JP] + p->YN[JM1])*Y0;
    Y3 = (-p->YN[JP1] + 27.0*p->YN[JP] -27.0*p->YN[JM1] + p->YN[JM2])*Y0;
    Y4 = (-p->YN[JP] + 27.0*p->YN[JM1] -27.0*p->YN[JM2] + p->YN[JM3])*Y0;
    
    Z0 = -0.5*p->ZN[KP2] + 13.0*p->ZN[KP1] - 13.0*p->ZN[KM1] + 0.5*p->ZN[KM2];
    Z1 = (-p->ZN[KP3] + 27.0*p->ZN[KP2] -27.0*p->ZN[KP1] + p->ZN[KP])*Z0;
    Z2 = (-p->ZN[KP2] + 27.0*p->ZN[KP1] -27.0*p->ZN[KP] + p->ZN[KM1])*Z0;
    Z3 = (-p->ZN[KP1] + 27.0*p->ZN[KP] -27.0*p->ZN[KM1] + p->ZN[KM2])*Z0;
    Z4 = (-p->ZN[KP] + 27.0*p->ZN[KM1] -27.0*p->ZN[KM2] + p->ZN[KM3])*Z0;
    

	a->M.p[n] = (1.0/X1 + 729.0/X2 + 729.0/X3 + 1.0/X4)*p->x_dir 
             + (1.0/Y1 + 729.0/Y2 + 729.0/Y3 + 1.0/Y4)*p->y_dir 
             + (1.0/Z1 + 729.0/Z2 + 729.0/Z3 + 1.0/Z4)*p->z_dir;
    
    
   	a->M.n[n] = -(27.0/X1 + 729.0/X2 + 27.0/X3)*p->x_dir;
	a->M.s[n] = -(27.0/X2 + 729.0/X3 + 27.0/X4)*p->x_dir;

	a->M.w[n] = -(27.0/Y1 + 729.0/Y2 + 27.0/Y3)*p->y_dir;
	a->M.e[n] = -(27.0/Y2 + 729.0/Y3 + 27.0/Y4)*p->y_dir;

	a->M.t[n] = -(27.0/Z1 + 729.0/Z2 + 27.0/Z3)*p->z_dir;
	a->M.b[n] = -(27.0/Z2 + 729.0/Z3 + 27.0/Z4)*p->z_dir;
    
    
    a->M.nn[n] = (27.0/X1 + 27.0/X2)*p->x_dir; 
    a->M.ss[n] = (27.0/X3 + 27.0/X4)*p->x_dir;
    
    a->M.ww[n] = (27.0/Y1 + 27.0/Y2)*p->y_dir;
    a->M.ee[n] = (27.0/Y3 + 27.0/Y4)*p->y_dir;
     
    a->M.tt[n] = (27.0/Z1 + 27.0/Z2)*p->z_dir; 
    a->M.bb[n] = (27.0/Z3 + 27.0/Z4)*p->z_dir;
                  
    a->rhsvec.V[n] +=   (f[FIp3JK]*(1.0/X1) + f[FIm3JK]*(1.0/X4))*p->x_dir
                    +  (f[FIJp3K]*(1.0/Y1) + f[FIJm3K]*(1.0/Y4))*p->y_dir
                    +  (f[FIJKp3]*(1.0/Z1) + f[FIJKm3]*(1.0/Z4))*p->z_dir;
        
	++n;
	}
    
    n=0;
	LOOP
	{
        if(p->flag4[FIJK]>0)
        {
    
            if(p->flag4[FIm1JK]<0)
            {
            a->rhsvec.V[n] -= a->M.s[n]*f[FIJK];
            a->M.s[n] = 0.0;
            }
            
            if(p->flag4[FIp1JK]<0)
            {
            a->rhsvec.V[n] -= a->M.n[n]*f[FIJK];
            a->M.n[n] = 0.0;
            }
            
            if(p->flag4[FIJm1K]<0)
            {
            a->rhsvec.V[n] -= a->M.e[n]*f[FIJK];
            a->M.e[n] = 0.0;
            }
            
            if(p->flag4[FIJp1K]<0)
            {
            a->rhsvec.V[n] -= a->M.w[n]*f[FIJK];
            a->M.w[n] = 0.0;
            }
            
            if(p->flag4[FIJKm1]<0)
            {
            a->rhsvec.V[n] -= a->M.b[n]*f[FIJKm1];
            a->M.t[n] = 0.0;
            }
            
            if(p->flag4[FIJKp1]<0)
            {
            a->rhsvec.V[n] -= a->M.t[n]*f[FIJKp1];
            a->M.t[n] = 0.0;
            }
            
            
            //--
            if(p->flag4[FIm2JK]<0)
            {
            a->rhsvec.V[n] -= a->M.ss[n]*f[FIJK];
            a->M.ss[n] = 0.0;
            }
            
            if(p->flag4[FIp2JK]<0)
            {
            a->rhsvec.V[n] -= a->M.nn[n]*f[FIJK];
            a->M.nn[n] = 0.0;
            }
            
            if(p->flag4[FIJm2K]<0)
            {
            a->rhsvec.V[n] -= a->M.ee[n]*f[FIJK];
            a->M.ee[n] = 0.0;
            }
            
            if(p->flag4[FIJp2K]<0)
            {
            a->rhsvec.V[n] -= a->M.ww[n]*f[FIJK];
            a->M.ww[n] = 0.0;
            }
            
            if(p->flag4[FIJKm2]<0)
            {
            a->rhsvec.V[n] -= a->M.bb[n]*f[FIJKm2];
            a->M.bb[n] = 0.0;
            }
            
            if(p->flag4[FIJKp2]<0)
            {
            a->rhsvec.V[n] -= a->M.tt[n]*f[FIJKp2];
            a->M.tt[n] = 0.0;
            }
            
        
        }
	++n;
	}
}

void poisson_4th::istart(lexer* p, fdm* a, field &apu, field &apv, field &apw, field &pcorr)
{
}

