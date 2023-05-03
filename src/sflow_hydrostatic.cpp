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

#include"sflow_hydrostatic.h"
#include"lexer.h"
#include"fdm2D.h" 
#include"ghostcell.h"
#include"poisson.h"
#include"solver2D.h"
#include"momentum.h"
#include"ioflow.h"
#include"patchBC_interface.h"

#define HX (fabs(b->hx(i,j))>1.0e-20?b->hx(i,j):1.0e20)
#define HXP (fabs(0.5*(b->WL(i,j)+b->WL(i+1,j)))>1.0e-20?0.5*(b->WL(i,j)+b->WL(i+1,j)):1.0e20)
#define HY (fabs(b->hy(i,j))>1.0e-20?b->hy(i,j):1.0e20)
 
sflow_hydrostatic::sflow_hydrostatic(lexer* p, fdm2D *b, patchBC_interface *ppBC)
{
    pBC = ppBC;
}

sflow_hydrostatic::~sflow_hydrostatic()
{
}

void sflow_hydrostatic::start(lexer *p, fdm2D *b, ghostcell *pgc, solver2D *psolv, ioflow *pflow, slice &P, slice &Q, slice &Pn, slice &Qn, slice &ws, slice &eta, double alpha)
{
    
}

void sflow_hydrostatic::ucorr(lexer* p, fdm2D* b, slice& uvel, slice &eta,double alpha)
{	
}

void sflow_hydrostatic::vcorr(lexer* p, fdm2D* b, slice& vvel, slice &eta,double alpha)
{	
}

void sflow_hydrostatic::wcorr(lexer* p, fdm2D* b,double alpha, slice &uvel, slice &vvel, slice &ws)
{	
}

void sflow_hydrostatic::wcalc(lexer* p, fdm2D* b,double alpha, slice &uvel, slice &vvel, slice &ws)
{	
}

void sflow_hydrostatic::upgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{
        if(p->A221==1)
        SLICELOOP1
        WETDRY1
        {
        b->F(i,j) -= fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) 
                                     - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j) )/(p->DXM);
        }
        
        
        if(p->A221==2)
        SLICELOOP1
        WETDRY1
        b->F(i,j) -= fabs(p->W22)*(1.0/HX)*
    
                    (0.5*(pow(eta(i+1,j),2.0) - pow(eta(i,j),2.0))/p->DXP[IP]
                    
                    + ((p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j))*b->depth(i+1,j) - (p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j))*b->depth(i,j))/p->DXP[IP]
                    
                    - 0.5*((p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j)) + (p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j)))*(b->depth(i+1,j)-b->depth(i,j))/p->DXP[IP]);
        
        if(p->B77==2)
        for(n=0;n<p->gcslout_count;n++)
        {
        i=p->gcslout[n][0]-1;
        j=p->gcslout[n][1];
        
        b->F(i,j) += fabs(p->W22)*(p->A223*eta(i+1,j) + (1.0-p->A223)*eta_n(i+1,j) 
                                     - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j) )/(p->DXM);
                                     
        b->F(i,j) -= fabs(p->W22)*(p->A223*(b->bed(i,j)-p->wd) + (1.0-p->A223)*(b->bed(i,j)-p->wd)
                                     - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j) )/(p->DXM);
                                     
        }
        
        pBC->patchBC_pressure2D_ugrad(p,b,eta,eta_n);
}

void sflow_hydrostatic::vpgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{
        if(p->A221==1)
        SLICELOOP2
        WETDRY2
        b->G(i,j) -= fabs(p->W22)*(p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1) 
                                 - p->A223*eta(i,j) - (1.0-p->A223)*eta_n(i,j) )/(p->DXM);
                                 
        if(p->A221==2)
        SLICELOOP1
        WETDRY2
        b->G(i,j) -= fabs(p->W22)*(1.0/HY)*
    
                    (0.5*(pow(eta(i,j+1),2.0) - pow(eta(i,j),2.0))/p->DYP[JP]
                    
                    + ((p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1))*b->depth(i,j+1) - (p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j))*b->depth(i,j))/p->DYP[JP]
                    
                    - 0.5*((p->A223*eta(i,j) + (1.0-p->A223)*eta_n(i,j)) + (p->A223*eta(i,j+1) + (1.0-p->A223)*eta_n(i,j+1)))*(b->depth(i,j+1)-b->depth(i,j))/p->DYP[JP]);
                                 
        pBC->patchBC_pressure2D_vgrad(p,b,eta,eta_n);
}

void sflow_hydrostatic::wpgrad(lexer*p, fdm2D* b, slice &eta, slice &eta_n)
{	    
}





