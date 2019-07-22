/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#define WLVL (fabs(a->WL(i,j))>1.0e-20?a->WL(i,j):1.0-20)

#define WLVLDRY (0.01*a->wd_criterion)

grid_sigma::grid_sigma(lexer *p) : Ex(p),Ey(p),Bx(p),By(p),Exx(p),Eyy(p),Bxx(p),Byy(p)
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
    a->wd_criterion=0.00005;
    
    /*
    if(p->A444==1)
    wd_criterion=p->A244_val;
    
    if(p->A445==1)
    wd_criterion=p->A245_val*p->dx;*/
    
    p->Darray(p->sig,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigx,p->imax*p->jmax*(p->kmax+1));
    p->Darray(p->sigy,p->imax*p->jmax*(p->kmax+1));
    
    p->Darray(p->sigxx,p->imax*p->jmax*(p->kmax+1));
    
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
    a->WL(i,j) = MAX(0.0,a->eta(i,j) + p->wd - a->bed(i,j));
    
    SLICELOOP4
	a->depth(i,j) = p->wd - a->bed(i,j);
    
    pgc->gcsl_start4(p,a->depth,50);
    
    SLICELOOP4
    p->sigz[IJ] = 1.0/WLVL;
}

void grid_sigma::sigma_update(lexer *p, fdm *a, ghostcell *pgc, slice &eta)
{
    // sigx
    FLOOP
    {
    if(a->wet(i,j)==1)
    p->sigx[FIJK] = (1.0 - p->sig[FIJK])*(Bx(i,j)/WLVL) - p->sig[FIJK]*(Ex(i,j)/WLVL);
    
    if(a->wet(i,j)==0)
    p->sigx[FIJK] = (1.0 - p->sig[FIJK])*(Bx(i,j)/WLVLDRY) - p->sig[FIJK]*(Bx(i,j)/WLVLDRY);
    }
    
    // sigy
    FLOOP
    {
    if(a->wet(i,j)==1)
    p->sigy[FIJK] = (1.0 - p->sig[FIJK])*(By(i,j)/WLVL) - p->sig[FIJK]*(Ey(i,j)/WLVL);
    
    if(a->wet(i,j)==0)
    p->sigy[FIJK] = (1.0 - p->sig[FIJK])*(By(i,j)/WLVLDRY) - p->sig[FIJK]*(By(i,j)/WLVLDRY);
    }
    
    // sigz
    SLICELOOP4
    {
    if(a->wet(i,j)==1)
    p->sigz[IJ] = 1.0/WLVL;
    
    if(a->wet(i,j)==0)
    p->sigz[IJ] = 1.0/WLVLDRY;
    }
    
    // sigxx
    FLOOP
    {
        if(a->wet(i,j)==1)
        {
        p->sigxx[FIJK] = ((1.0 - p->sig[FIJK])/WLVL)*(Bxx(i,j) - pow(Bx(i,j),2.0)/WLVL) // xx
        
                      - (p->sig[FIJK]/WLVL)*(Exx(i,j) - pow(Ex(i,j),2.0)/WLVL)
                      
                      - (p->sigx[FIJK]/WLVL)*(Bx(i,j) + Ex(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(Bx(i,j)*Ex(i,j))
                      
                      
                      + ((1.0 - p->sig[FIJK])/WLVL)*(Byy(i,j) - pow(By(i,j),2.0)/WLVL) // yy
        
                      - (p->sig[FIJK]/WLVL)*(Eyy(i,j) - pow(Ey(i,j),2.0)/WLVL)
                      
                      - (p->sigy[FIJK]/WLVL)*(By(i,j) + Ey(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVL,2.0))*(By(i,j)*Ey(i,j));
        }
                      
                      
        if(a->wet(i,j)==0)
        {
        p->sigxx[FIJK] = ((1.0 - p->sig[FIJK])/WLVLDRY)*(Bxx(i,j) - pow(Bx(i,j),2.0)/WLVLDRY) // xx
        
                      - (p->sig[FIJK]/WLVLDRY)*(Bxx(i,j) - pow(Bx(i,j),2.0)/WLVLDRY)
                      
                      - (p->sigx[FIJK]/WLVLDRY)*(Bx(i,j) + Bx(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVLDRY,2.0))*(Bx(i,j)*Bx(i,j))
                      
                      
                      + ((1.0 - p->sig[FIJK])/WLVLDRY)*(Byy(i,j) - pow(By(i,j),2.0)/WLVLDRY) // yy
        
                      - (p->sig[FIJK]/WLVLDRY)*(Byy(i,j) - pow(By(i,j),2.0)/WLVLDRY)
                      
                      - (p->sigy[FIJK]/WLVLDRY)*(By(i,j) + By(i,j))
                      
                      - ((1.0 - 2.0*p->sig[FIJK])/pow(WLVLDRY,2.0))*(By(i,j)*By(i,j));
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
    p->ZSN[FIJK] = p->ZN[KP]*(eta(i,j) + p->wd);
        
}

double grid_sigma::sigmax(lexer *p, field &f, int ipol)
{    
    if(ipol==1)
    sig = 0.25*(p->sigx[FIJKm1] + p->sigx[FIp1JKm1] + p->sigx[FIJK] + p->sigx[FIp1JK])
        
        *((f(i,j,k+1)-f(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]));

    if(ipol==2)
    sig = 0.25*(p->sigx[FIJKm1] + p->sigx[FIJp1Km1] + p->sigx[FIJK] + p->sigx[FIJp1K])
    
        *((f(i,j,k+1)-f(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]));
    
    if(ipol==3)
    sig = p->sigx[FIJKm1]
    
        *((f(i,j,k+1)-f(i,j,k-1))/(p->DZN[KP]+p->DZN[KM1]));

    if(ipol==4)
    sig = 0.5*(p->sigx[FIJKm1] + p->sigx[FIJK])
    
        *((f(i,j,k+1)-f(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]));

    return sig;
}


double grid_sigma::sigmay(lexer *p, field &f, int ipol)
{  
    if(ipol==1)
    sig = 0.25*(p->sigy[FIJKm1] + p->sigy[FIp1JKm1] + p->sigy[FIJK] + p->sigy[FIp1JK])
    
        *((f(i,j,k+1)-f(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]));

    if(ipol==2)
    sig = 0.25*(p->sigy[FIJKm1] + p->sigy[FIJp1Km1] + p->sigy[FIJK] + p->sigy[FIJp1K])
    
        *((f(i,j,k+1)-f(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]));

    if(ipol==3)
    sig = p->sigy[FIJKm1]
    
        *((f(i,j,k+1)-f(i,j,k-1))/(p->DZN[KP]+p->DZN[KM1]));

    if(ipol==4)
    sig = 0.5*(p->sigy[FIJKm1] + p->sigy[FIJK])
    
        *((f(i,j,k+1)-f(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]));

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



