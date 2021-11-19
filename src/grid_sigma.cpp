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

#include"grid_sigma.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fnpf_ddx_cds2.h"
#include"fnpf_ddx_cds4.h"
#include"fnpf_cds2.h"
#include"fnpf_cds4.h"
#include"grid_sigma_data.h"

#define WLVL (fabs(a->WL(i,j))>1.0e-20?a->WL(i,j):1.0-20)

#define WLVLDRY (0.01*a->wd_criterion)

grid_sigma::grid_sigma(lexer *p) 
{
    
}

grid_sigma::~grid_sigma()
{
}

void grid_sigma::sigma_coord_ini(lexer *p)
{
    double L;
    
    L = p->ZN[p->knoz+marge] - p->ZN[0+marge];

    
    for(k=-marge;k<p->knoz+marge;++k)
    p->ZN[KP] = (p->ZN[KP]-p->ZN[0+marge])/L;
}
    
void grid_sigma::sigma_ini(lexer *p, fdm *a, ghostcell *pgc, slice &eta)
{	
    
    // generate Ex,Bx slices
    pd = new grid_sigma_data(p);
    
    // generate discretization 
    if(p->A312==2)
    {
    pddx = new fnpf_ddx_cds2(p);
    pdx = new fnpf_cds2(p);
    }
    
    if(p->A312==3)
    {
    pddx = new fnpf_ddx_cds4(p);
    pdx = new fnpf_cds4(p);
    }
    
    SLICELOOP4
    a->wet(i,j)=1;
    
    pgc->gcsl_start4int(p,a->wet,50);
    
    
    a->wd_criterion=0.00005;
    
    p->Darray(p->sig, p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigx,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigy,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigz,p->imax*p->jmax);
    p->Darray(p->sigt,p->imax*p->jmax*(p->kmax+1));
    
    p->Darray(p->sigxx,p->imax*p->jmax*(p->kmax+1));
    

    // p->sig[FIJK]
    FLOOP
    p->sig[FIJK] =  p->ZN[KP];
    
    // bc
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sig[FIJKm1] = p->ZN[KM1];
            p->sig[FIJKm2] = p->ZN[KM2];
            p->sig[FIJKm3] = p->ZN[KM3];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sig[FIJKp1] = p->ZN[KP1];
            p->sig[FIJKp2] = p->ZN[KP2];
            p->sig[FIJKp3] = p->ZN[KP3];
        } 
    }

    
    SLICELOOP4
	a->bed(i,j) = p->bed[IJ];
    
    
    SLICELOOP4
    a->WL(i,j) = MAX(0.0, a->eta(i,j) + p->wd - a->bed(i,j));
    
    
    SLICELOOP4
	a->depth(i,j) = p->wd - a->bed(i,j);
    
    pgc->gcsl_start4(p,a->depth,50);
    
    SLICELOOP4
    p->sigz[IJ] = 1.0/WLVL;
    
    SLICELOOP4
    p->sigt[FIJK] = 0.0;

}

void grid_sigma::sigma_update(lexer *p, fdm *a, ghostcell *pgc, slice &eta)
{
    SLICELOOP4
    {
    a->WL_n(i,j) = a->WL(i,j);
    a->WL(i,j) = MAX(0.0, a->eta(i,j) + p->wd - a->bed(i,j));
    }
    
    // calculate: Ex,Ey,Exx,Eyy
    // 3D
    if(p->i_dir==1 && p->j_dir==1)
    SLICELOOP4
    {
    pd->Ex(i,j) = pdx->sx(p,a->eta,1.0);
    pd->Ey(i,j) = pdx->sy(p,a->eta,1.0);
    
    pd->Exx(i,j) = pddx->sxx(p,a->eta);
    pd->Eyy(i,j) = pddx->syy(p,a->eta);
    }
    
    // 2D
    if(p->i_dir==1 && p->j_dir==0)
    SLICELOOP4
    {
    pd->Ex(i,j) = pdx->sx(p,a->eta,1.0);    
    pd->Exx(i,j) = pddx->sxx(p,a->eta);
    }
    
    pgc->gcsl_start4(p,pd->Ex,1);
    pgc->gcsl_start4(p,pd->Ey,1);
    
    // calculate: Bx,By,Bxx,Byy
    // 3D
    if(p->i_dir==1 && p->j_dir==1)
    SLICELOOP4
    {
    pd->Bx(i,j) = pdx->sx(p,a->depth,1.0);
    pd->By(i,j) = pdx->sy(p,a->depth,1.0);
    
    pd->Bxx(i,j) = pddx->sxx(p,a->depth);
    pd->Byy(i,j) = pddx->syy(p,a->depth);
    }

    // 2D
    if(p->i_dir==1 && p->j_dir==0)
    SLICELOOP4
    {
    pd->Bx(i,j) = pdx->sx(p,a->depth,1.0);    
    pd->Bxx(i,j) = pddx->sxx(p,a->depth);
    }
    
    pgc->gcsl_start4(p,pd->Bx,1);
    pgc->gcsl_start4(p,pd->By,1);
    
    
    // sigx
    FLOOP
    {
    if(a->wet(i,j)==1)
    p->sigx[FIJK] = (1.0 - p->sig[FIJK])*(pd->Bx(i,j)/WLVL) - p->sig[FIJK]*(pd->Ex(i,j)/WLVL);
    
    if(a->wet(i,j)==0)
    p->sigx[FIJK] = (1.0 - p->sig[FIJK])*(pd->Bx(i,j)/WLVLDRY) - p->sig[FIJK]*(pd->Bx(i,j)/WLVLDRY);
    }
    
    // sigy
    FLOOP
    {
    if(a->wet(i,j)==1)
    p->sigy[FIJK] = (1.0 - p->sig[FIJK])*(pd->By(i,j)/WLVL) - p->sig[FIJK]*(pd->Ey(i,j)/WLVL);
    
    if(a->wet(i,j)==0)
    p->sigy[FIJK] = (1.0 - p->sig[FIJK])*(pd->By(i,j)/WLVLDRY) - p->sig[FIJK]*(pd->By(i,j)/WLVLDRY);
    }
    
    // sigz
    SLICELOOP4
    {
    if(a->wet(i,j)==1)
    p->sigz[IJ] = 1.0/WLVL;
    
    if(a->wet(i,j)==0)
    p->sigz[IJ] = 1.0/WLVLDRY;
    }
    
    // sigt
    FLOOP
    {
    if(a->wet(i,j)==1)
    p->sigt[FIJK] = -(p->sig[FIJK]/WLVL)*(a->WL(i,j)-a->WL_n(i,j))/p->dt;
    
    if(a->wet(i,j)==0)
    p->sigt[FIJK] = 0.0;
    }
    
    // sigxx
    FLOOP
    {
        if(a->wet(i,j)==1)
        {
        p->sigxx[FIJK] = ((1.0 - p->sig[FIJK])/WLVL)*(pd->Bxx(i,j) - pow(pd->Bx(i,j),2.0)/WLVL) // xx
        
                      - (p->sig[FIJK]/WLVL)*(pd->Exx(i,j) - pow(pd->Ex(i,j),2.0)/WLVL)
                      
                      - (p->sigx[FIJK]/WLVL)*(pd->Bx(i,j) + pd->Ex(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(pd->Bx(i,j)*pd->Ex(i,j))
                      
                      
                      + ((1.0 - p->sig[FIJK])/WLVL)*(pd->Byy(i,j) - pow(pd->By(i,j),2.0)/WLVL) // yy
        
                      - (p->sig[FIJK]/WLVL)*(pd->Eyy(i,j) - pow(pd->Ey(i,j),2.0)/WLVL)
                      
                      - (p->sigy[FIJK]/WLVL)*(pd->By(i,j) + pd->Ey(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(pd->By(i,j)*pd->Ey(i,j));
        }
                      
                      
        if(a->wet(i,j)==0)
        {
        p->sigxx[FIJK] = ((1.0 - p->sig[FIJK])/WLVLDRY)*(pd->Bxx(i,j) - pow(pd->Bx(i,j),2.0)/WLVLDRY) // xx
        
                      - (p->sig[FIJK]/WLVLDRY)*(pd->Bxx(i,j) - pow(pd->Bx(i,j),2.0)/WLVLDRY)
                      
                      - (p->sigx[FIJK]/WLVLDRY)*(pd->Bx(i,j) + pd->Bx(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVLDRY,2.0))*(pd->Bx(i,j)*pd->Bx(i,j))
                      
                      
                      + ((1.0 - p->sig[FIJK])/WLVLDRY)*(pd->Byy(i,j) - pow(pd->By(i,j),2.0)/WLVLDRY) // yy
        
                      - (p->sig[FIJK]/WLVLDRY)*(pd->Byy(i,j) - pow(pd->By(i,j),2.0)/WLVLDRY)
                      
                      - (p->sigy[FIJK]/WLVLDRY)*(pd->By(i,j) + pd->By(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVLDRY,2.0))*(pd->By(i,j)*pd->By(i,j));
        }
    }
    
    // sig BC
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigx[FIJKm1] = p->sigx[KP];
            p->sigx[FIJKm2] = p->sigx[KP];
            p->sigx[FIJKm3] = p->sigx[KP];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigx[FIJKp1] = p->sigx[KP];
            p->sigx[FIJKp2] = p->sigx[KP];
            p->sigx[FIJKp3] = p->sigx[KP];
        } 
    }
    
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigy[FIJKm1] = p->sigy[KP];
            p->sigy[FIJKm2] = p->sigy[KP];
            p->sigy[FIJKm3] = p->sigy[KP];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigy[FIJKp1] = p->sigy[KP];
            p->sigy[FIJKp2] = p->sigy[KP];
            p->sigy[FIJKp3] = p->sigy[KP];
        } 
    }
    
    SLICELOOP4
    {
        k=0;
        if(p->nb5==-2)
        {
            p->sigxx[FIJKm1] = p->sigxx[KP];
            p->sigxx[FIJKm2] = p->sigxx[KP];
            p->sigxx[FIJKm3] = p->sigxx[KP];
        }
        
        k=p->knoz;
        if(p->nb6==-2)
        {
            p->sigxx[FIJKp1] = p->sigxx[KP];
            p->sigxx[FIJKp2] = p->sigxx[KP];
            p->sigxx[FIJKp3] = p->sigxx[KP];
        } 
    }
    
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
    }
    
    
    FLOOP
    p->ZSN[FIJK] = p->ZN[KP]*a->WL(i,j) + a->bed(i,j);
    
    
    LOOP
    p->ZSP[IJK]  = p->ZP[KP]*a->WL(i,j) + a->bed(i,j);
        
}

double grid_sigma::sigmax(lexer *p, field &f, int ipol)
{    
    if(ipol==1)
    sig = 0.25*(p->sigx[FIJKm1] + p->sigx[FIp1JKm1] + p->sigx[FIJK] + p->sigx[FIp1JK]);

    if(ipol==2)
    sig = 0.25*(p->sigx[FIJKm1] + p->sigx[FIJp1Km1] + p->sigx[FIJK] + p->sigx[FIJp1K]);
    
    if(ipol==3)
    sig = p->sigx[FIJKm1];

    if(ipol==4)
    sig = 0.5*(p->sigx[FIJKm1] + p->sigx[FIJK]);

    return sig;
}

double grid_sigma::sigmay(lexer *p, field &f, int ipol)
{  
    if(ipol==1)
    sig = 0.25*(p->sigy[FIJKm1] + p->sigy[FIp1JKm1] + p->sigy[FIJK] + p->sigy[FIp1JK]);

    if(ipol==2)
    sig = 0.25*(p->sigy[FIJKm1] + p->sigy[FIJp1Km1] + p->sigy[FIJK] + p->sigy[FIJp1K]);
    
    if(ipol==3)
    sig = p->sigy[FIJKm1];

    if(ipol==4)
    sig = 0.5*(p->sigy[FIJKm1] + p->sigy[FIJK]);

    return sig;
}

double grid_sigma::sigmaz(lexer *p, field &f, int ipol)
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

double grid_sigma::sigmat(lexer *p, field &f, int ipol)
{    
    if(ipol==1)
    sig = 0.25*(p->sigt[FIJKm1] + p->sigt[FIp1JKm1] + p->sigt[FIJK] + p->sigt[FIp1JK]);

    if(ipol==2)
    sig = 0.25*(p->sigt[FIJKm1] + p->sigt[FIJp1Km1] + p->sigt[FIJK] + p->sigt[FIJp1K]);
    
    if(ipol==3)
    sig = p->sigt[FIJKm1];

    if(ipol==4)
    sig = 0.5*(p->sigt[FIJKm1] + p->sigt[FIJK]);

    return sig;
}



