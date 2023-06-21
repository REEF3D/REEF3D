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

#include"sflow_eta.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"sflow_eta_weno.h"
#include"sflow_hxy_weno.h"
#include"sflow_hxy_hires.h"
#include"sflow_hxy_cds.h"
#include"sflow_hxy_fou.h"
#include"patchBC_interface.h"

sflow_eta::sflow_eta(lexer *p, fdm2D *b , ghostcell *pgc, patchBC_interface *ppBC) : Lab(p)
{   
    pBC = ppBC;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
    
    
    p->phimean=p->F60;
    
    SLICELOOP4
    b->eta(i,j) = 0.0;

    pgc->gcsl_start4(p,b->eta,gcval_eta);
	
	pconvec = new sflow_eta_weno(p);
	
	if(p->A241==1)
	phxy = new sflow_hxy_fou(p,pBC);
	
	if(p->A241==2)
	phxy = new sflow_hxy_cds(p,pBC);
	
	if(p->A241==4)
	phxy = new sflow_hxy_weno(p,pBC);
    
    if(p->A241>=6)
	phxy = new sflow_hxy_hires(p,pBC,p->A241);
    
    wd_criterion=0.00005;
    
    if(p->A244==1)
    wd_criterion=p->A244_val;
    
    if(p->A245==1)
    wd_criterion=p->A245_val*p->DXM;
    
}

sflow_eta::~sflow_eta()
{
}

void sflow_eta::start(lexer* p, fdm2D* b, ghostcell* pgc, ioflow* pflow, slice &P, slice &Q, double alpha)
{
    starttime=pgc->timer();
    
    if(p->A240==2)
    {
    wetdry(p,b,pgc,b->P,b->Q,b->ws);
    depth_update(p,b,pgc,b->P,b->Q,b->ws,b->eta);
    
    SLICELOOP4
     b->eta(i,j)  =      b->eta(i,j) 
    
                -      p->dt*(b->P(i,j)*b->hx(i,j) - b->P(i-1,j)*b->hx(i-1,j)
                       +      b->Q(i,j)*b->hy(i,j) - b->Q(i,j-1)*b->hy(i,j-1))/p->DXM;
                
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    depth_update(p,b,pgc,b->P,b->Q,b->ws,b->eta);
    breaking(p,b,pgc,b->eta,b->eta,1.0);
    pflow->eta_relax(p,pgc,b->eta);
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    }
}

void sflow_eta::depth_update(lexer *p, fdm2D *b , ghostcell *pgc, slice &P, slice &Q, slice &ws, slice &etark)
{
    double factor=1.0;
    
	// cell center
	SLICELOOP4
	b->depth(i,j) = p->wd - b->bed(i,j);
    
    
    SLICELOOP4
    if(etark(i,j)< -p->wd  + b->bed(i,j)-factor*wd_criterion+1.0e-20)
    etark(i,j) = -p->wd  + b->bed(i,j)-1.0*wd_criterion - 1.0e-15;

    if(p->A243>=2)
    {
        SLICELOOP4
        if(etark(i,j)< -p->wd  + b->bed(i,j)-factor*wd_criterion)
        etark(i,j) = -p->wd  + b->bed(i,j)-factor*wd_criterion;
        
        
        SLICELOOP4
        if(etark(i,j)>=-p->wd + b->bed(i,j)+wd_criterion)
        {
            if(etark(i+1,j)<-p->wd + b->bed(i+1,j)-factor*wd_criterion)
            etark(i+1,j) = etark(i,j);
            
            if(etark(i-1,j)<-p->wd + b->bed(i-1,j)-factor*wd_criterion)
            etark(i-1,j) = etark(i,j);
            
            if(etark(i,j+1)<-p->wd + b->bed(i,j+1)-factor*wd_criterion)
            etark(i,j+1) = etark(i,j);

            if(etark(i,j-1)<-p->wd + b->bed(i,j-1)-factor*wd_criterion)
            etark(i,j-1) = etark(i,j);
        }
    }
	
	pgc->gcsl_start4(p,b->depth,50);

	SLICELOOP4
	b->hp(i,j) = MAX(etark(i,j) + p->wd - b->bed(i,j),0.0);
    
    
	pgc->gcsl_start4(p,b->hp,gcval_eta);
	
	phxy->start(p,b->hx,b->hy,b->depth,p->wet,etark,P,Q);
    
    
    SLICELOOP1
    b->hx(i,j) = MAX(b->hx(i,j), 0.0);
    
    SLICELOOP2
    b->hy(i,j) = MAX(b->hy(i,j), 0.0);
    

	pgc->gcsl_start1(p,b->hx,gcval_eta);    
	pgc->gcsl_start2(p,b->hy,gcval_eta);
    pgc->gcsl_start4(p,b->depth,1);
    

    if(p->A243>=1)
    wetdry(p,b,pgc,P,Q,ws);
}
	
void sflow_eta::ini(lexer *p, fdm2D *b , ghostcell *pgc, ioflow *pflow)
{
    
    p->phimean=p->F60;
    
    SLICELOOP4
    b->eta(i,j) = 0.0;
    
    pflow->eta_relax(p,pgc,b->eta);
    
    pgc->gcsl_start4(p,b->eta,gcval_eta);

}

