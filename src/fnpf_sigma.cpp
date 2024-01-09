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

#include"fnpf_sigma.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_fsf.h"

#define WLVL (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0e-5) // keep as is for wetting-drying

#define WLVLDRY (0.01*c->wd_criterion)

fnpf_sigma::fnpf_sigma(lexer *p, fdm_fnpf *c, ghostcell *pgc) 
{
}

fnpf_sigma::~fnpf_sigma()
{
}

void fnpf_sigma::sigma_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc, fnpf_fsf *pf, slice &eta)
{	
    
    
    FLOOP
    p->sig[FIJK] =  p->ZN[KP];
    
    // bc
    SLICELOOP4
    {
        k=0;
    
            p->sig[FIJKm1] = p->ZN[KM1];
            p->sig[FIJKm2] = p->ZN[KM2];
            p->sig[FIJKm3] = p->ZN[KM3];
        
        k=p->knoz;

            p->sig[FIJKp1] = p->ZN[KP1];
            p->sig[FIJKp2] = p->ZN[KP2];
            p->sig[FIJKp3] = p->ZN[KP3];
    }
    
    
    SLICELOOP4
	c->bed(i,j) = p->bed[IJ];
    
    SLICELOOP4
    c->WL(i,j) = MAX(0.0,c->eta(i,j) + p->wd - c->bed(i,j));
    
    SLICEBASELOOP
    {
    PSLICECHECK4
	c->depth(i,j) = p->wd - c->bed(i,j);
    
    SSLICECHECK4
	c->depth(i,j) = p->wd - p->bed[IJ];
    }
    
    pgc->gcsl_start4(p,c->depth,50);
    
    SLICELOOP4
    p->sigz[IJ] = 1.0/WLVL;
}

void fnpf_sigma::sigma_update(lexer *p, fdm_fnpf *c, ghostcell *pgc, fnpf_fsf *pf, slice &eta)
{
    // sigx
    FBASELOOP
    p->sigx[FIJK] = (1.0 - p->sig[FIJK])*(c->Bx(i,j)/WLVL) - p->sig[FIJK]*(c->Ex(i,j)/WLVL);
    
    // sigy
    FBASELOOP
    p->sigy[FIJK] = (1.0 - p->sig[FIJK])*(c->By(i,j)/WLVL) - p->sig[FIJK]*(c->Ey(i,j)/WLVL);
    
    // sigz
    SLICEBASELOOP
    {
        PSLICECHECK4
        p->sigz[IJ] = 1.0/WLVL;
        
        SSLICECHECK4
        p->sigz[IJ] = 1.0/(p->wd-p->bed[IJ]);
    }
    
    // sigxx
    FBASELOOP
    {
    p->sigxx[FIJK] = ((1.0 - p->sig[FIJK])/WLVL)*(c->Bxx(i,j) - pow(c->Bx(i,j),2.0)/WLVL) // xx
    
                  - (p->sig[FIJK]/WLVL)*(c->Exx(i,j) - pow(c->Ex(i,j),2.0)/WLVL)
                  
                  - (p->sigx[FIJK]/WLVL)*(c->Bx(i,j) + c->Ex(i,j))
                  
                  - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(c->Bx(i,j)*c->Ex(i,j))
                  
                  
                  + ((1.0 - p->sig[FIJK])/WLVL)*(c->Byy(i,j) - pow(c->By(i,j),2.0)/WLVL) // yy
    
                  - (p->sig[FIJK]/WLVL)*(c->Eyy(i,j) - pow(c->Ey(i,j),2.0)/WLVL)
                  
                  - (p->sigy[FIJK]/WLVL)*(c->By(i,j) + c->Ey(i,j))
                  
                  - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(c->By(i,j)*c->Ey(i,j));
    }
    
    /*
    FBASELOOP
    {
    p->sigxx[FIJK] = - 2.0*((1.0 - p->sig[FIJK])*pow(c->Bx(i,j),2.0)/pow(WLVL,2.0))
    
                     + 2.0*(2.0*p->sig[FIJK]-1.0)*(c->Ex(i,j)*c->Bx(i,j))/pow(WLVL,2.0)
    
                    + (1.0 - p->sig[FIJK])*c->Bxx(i,j)/WLVL + 2.0*p->sig[FIJK]*pow(c->Ex(i,j),2.0)/pow(WLVL,2.0)
                    
                    - p->sig[FIJK]*c->Exx(i,j) /WLVL;              
    }
    */

    // sig BC
    SLICELOOP4
    {
        k=0;
            p->sigx[FIJKm1] = p->sigx[FIJK];
            p->sigx[FIJKm2] = p->sigx[FIJK];
            p->sigx[FIJKm3] = p->sigx[FIJK];
            
            
        
        k=p->knoz;
            p->sigx[FIJKp1] = p->sigx[FIJK];
            p->sigx[FIJKp2] = p->sigx[FIJK];
            p->sigx[FIJKp3] = p->sigx[FIJK];
    }
    
    SLICELOOP4
    {
        k=0;
            p->sigy[FIJKm1] = p->sigy[FIJK];
            p->sigy[FIJKm2] = p->sigy[FIJK];
            p->sigy[FIJKm3] = p->sigy[FIJK];
        
        k=p->knoz;
            p->sigy[FIJKp1] = p->sigy[FIJK];
            p->sigy[FIJKp2] = p->sigy[FIJK];
            p->sigy[FIJKp3] = p->sigy[FIJK];
    }
    
    SLICELOOP4
    {
        k=0;
            p->sigxx[FIJKm1] = p->sigxx[FIJK];
            p->sigxx[FIJKm2] = p->sigxx[FIJK];
            p->sigxx[FIJKm3] = p->sigxx[FIJK];

        
        k=p->knoz;
            p->sigxx[FIJKp1] = p->sigxx[FIJK];
            p->sigxx[FIJKp2] = p->sigxx[FIJK];
            p->sigxx[FIJKp3] = p->sigxx[FIJK];
    }
    /*
    k=p->knoz;
    SLICELOOP4
    {
        if(p->flag7[FIm1JK]<0 || i==0)
        p->sigz[Im1J] = p->sigz[IJ];
        
        if(p->flag7[FIp1JK]<0 || i==p->knox-1)
        p->sigz[Ip1J] = p->sigz[IJ];
        
        if(p->flag7[FIJm1K]<0 || j==0)
        p->sigz[IJm1] = p->sigz[IJ];
        
        if(p->flag7[FIJp1K]<0 || j==p->knoy-1)
        p->sigz[IJp1] = p->sigz[IJ];
    }*/
    
    
    FLOOP
    {
    FPCHECK
    p->ZSN[FIJK] = p->ZN[KP]*c->WL(i,j) + c->bed(i,j); 
    
    FSCHECK
    p->ZSN[FIJK] = p->ZN[KP]*(p->wd);
    }
}




