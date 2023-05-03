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

#include"grid_sigma.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"fnpf_ddx_cds2.h"
#include"fnpf_ddx_cds4.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"
#include"grid_sigma_data.h"

#define WLVL (fabs(d->WL(i,j))>1.0e-20?d->WL(i,j):1.0-20)

grid_sigma::grid_sigma(lexer *p) 
{
}

grid_sigma::~grid_sigma()
{
}

void grid_sigma::sigma_coord_ini(lexer *p)
{
    double L, ZN0temp;
    
    L = p->ZN[p->knoz+marge] - p->ZN[0+marge];
    
    ZN0temp = p->ZN[0+marge];
    
    for(k=-marge;k<p->knoz+marge;++k)
    {
    p->ZN[KP] = (p->ZN[KP]-ZN0temp)/L;
    }
}
    
void grid_sigma::sigma_ini(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &eta)
{	
    // generate Ex,Bx slices
    pd = new grid_sigma_data(p);
    
    // generate discretization 
    if(p->A312==1||p->A312==2)
    {
    pddx = new fnpf_ddx_cds2(p);
    pdx  = new fnpf_cds2(p);
    }
    
    if(p->A312==3)
    {
    pddx = new fnpf_ddx_cds4(p);
    pdx  = new fnpf_cds4(p);
    }
    
    SLICELOOP4
    p->wet[IJ]=1;
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    
    d->wd_criterion=0.00005;
    
    p->Darray(p->sig, p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigx,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigy,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigz,p->imax*p->jmax);
    p->Darray(p->sigt,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigx1,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigy2,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigx4,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigy4,p->imax*p->jmax*(p->kmax+1));
    
    p->Darray(p->sigxx,p->imax*p->jmax*(p->kmax+1));
    

    FLOOP
    p->sig[FIJK] =  p->ZN[KP];
    
    // bc
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sig[FIJKm1] = p->ZN[FIJK];
            p->sig[FIJKm2] = p->ZN[FIJK];
            p->sig[FIJKm3] = p->ZN[FIJK];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sig[FIJKp1] = p->ZN[FIJK];
            p->sig[FIJKp2] = p->ZN[FIJK];
            p->sig[FIJKp3] = p->ZN[FIJK];
        } 
    }
    
    pgc->start7S(p,p->sig,1);

    
    SLICELOOP4
	d->bed(i,j) = p->bed[IJ];
    
    
    SLICELOOP4
    d->WL(i,j) = MAX(0.0, d->eta(i,j) + p->wd - d->bed(i,j));
    
    
    SLICELOOP4
	d->depth(i,j) = p->wd - d->bed(i,j);
    
    SLICELOOP4
    {
    d->Bx(i,j) = 0.0;
    d->By(i,j) = 0.0;
    }
    
    pgc->gcsl_start4(p,d->depth,50);
    pgc->gcsl_start4(p,d->Bx,50);
    pgc->gcsl_start4(p,d->By,50);
    
    SLICELOOP4
    p->sigz[IJ] = 1.0/WLVL;
    
    SLICELOOP4
    p->sigt[FIJK] = 0.0;

}

double grid_sigma::sigmax(lexer *p, int ipol)
{    
    if(ipol==1)
    sig = 0.25*(p->sigx[FIJK] + p->sigx[FIp1JK] + p->sigx[FIJKp1] + p->sigx[FIp1JKp1]);

    if(ipol==2)
    sig = 0.25*(p->sigx[FIJK] + p->sigx[FIJp1K] + p->sigx[FIJKp1] + p->sigx[FIJp1Kp1]);
    
    if(ipol==3)
    sig = p->sigx[FIJKp1];

    if(ipol==4)
    sig = 0.5*(p->sigx[FIJK] + p->sigx[FIJKp1]);

    return sig;
}

double grid_sigma::sigmay(lexer *p, int ipol)
{  
    if(ipol==1)
    sig = 0.25*(p->sigy[FIJK] + p->sigy[FIp1JK] + p->sigy[FIJKp1] + p->sigy[FIp1JKp1]);

    if(ipol==2)
    sig = 0.25*(p->sigy[FIJK] + p->sigy[FIJp1K] + p->sigy[FIJKp1] + p->sigy[FIJp1Kp1]);
    
    if(ipol==3)
    sig = p->sigy[FIJKp1];

    if(ipol==4)
    sig = 0.5*(p->sigy[FIJK] + p->sigy[FIJKp1]);

    return sig;
}

double grid_sigma::sigmaz(lexer *p, int ipol)
{    
    if(ipol==1)
    sig = 0.5*(p->sigz[IJ] + p->sigz[Ip1J]);

    if(ipol==2)
    sig = 0.5*(p->sigz[IJ] + p->sigz[IJp1]);

    if(ipol==3)
    sig = p->sigz[IJ];

    if(ipol==4)
    sig = p->sigz[IJ];

    return sig;
}

double grid_sigma::sigmat(lexer *p, int ipol)
{    
    if(ipol==1)
    sig = 0.25*(p->sigt[FIJK] + p->sigt[FIp1JK] + p->sigt[FIJKp1] + p->sigt[FIp1JKp1]);

    if(ipol==2)
    sig = 0.25*(p->sigt[FIJK] + p->sigt[FIJp1K] + p->sigt[FIJKp1] + p->sigt[FIJp1Kp1]);
    
    if(ipol==3)
    sig = p->sigt[FIJKp1];

    if(ipol==4)
    sig = 0.5*(p->sigt[FIJK] + p->sigt[FIJKp1]);

    return sig;
}


