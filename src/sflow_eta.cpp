/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

	if(p->A241==1)
	phxy = new sflow_hxy_fou(p,pBC);
	
	if(p->A241==2)
	phxy = new sflow_hxy_cds(p,pBC);
	
	if(p->A241==4)
	phxy = new sflow_hxy_weno(p,pBC);
    
    if(p->A241>=6)
	phxy = new sflow_hxy_hires(p,pBC,p->A241);
    
    wd_criterion=p->A244;
    
    p->Iarray(temp,p->imax*p->jmax);
    
}

sflow_eta::~sflow_eta()
{
}

void sflow_eta::start(lexer* p, fdm2D* b, ghostcell* pgc, ioflow* pflow, slice &P, slice &Q, double alpha)
{
}

void sflow_eta::disc(lexer *p, fdm2D *b , ghostcell *pgc, slice &P, slice &Q, slice &ws, slice &eta)
{
    double factor=1.0;
    double eps=0.0;
    
	phxy->start(p,b->hx,b->hy,b->depth,p->wet,eta,P,Q);
    
    SLICELOOP1
    b->hx(i,j) = MAX(b->hx(i,j), 0.0);
    
    SLICELOOP2
    b->hy(i,j) = MAX(b->hy(i,j), 0.0);
    

	pgc->gcsl_start1(p,b->hx,gcval_eta);    
	pgc->gcsl_start2(p,b->hy,gcval_eta);  
}

void sflow_eta::depth_update(lexer *p, fdm2D *b , ghostcell *pgc, slice &P, slice &Q, slice &ws, slice &eta)
{
    double factor=1.0;
    
	// cell center
	SLICELOOP4
	b->depth(i,j) = p->wd - b->bed(i,j);
    
    pgc->gcsl_start4(p,b->depth,1);
    
    
    // set fsf outflow
    double wsfout=p->phimean;
    double f=1.0;
    
    if(p->F62>1.0e-20)
    {
        if(p->F64==0)
        wsfout=p->F62;
        
        if(p->F64>0)
        {
        if(p->count<p->F64)
        f = 0.5*cos(PI + PI*double(p->count)/double(p->F64)) + 0.5;
        
        if(p->count>=p->F64)
        f = 1.0;
        
        wsfout = f*p->F62 + (1.0-f)*p->F60;
        }
    }
    
    
    if(p->F50==2 || p->F50==3)
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
    
    if(p->wet[IJ]==1)
    eta(i,j) = wsfout-p->wd;
    }
    
    // wetdry
    if(p->A243==1)
    wetdry(p,b,pgc,eta,P,Q,ws);

    if(p->A243==2)
    wetdry_nb(p,b,pgc,eta,P,Q,ws);
}
	
void sflow_eta::ini(lexer *p, fdm2D *b , ghostcell *pgc, ioflow *pflow)
{
    p->phimean=p->F60;
    
    SLICELOOP4
    b->eta(i,j) = 0.0;
    
    pflow->eta_relax(p,pgc,b->eta);
    
    pgc->gcsl_start4(p,b->eta,gcval_eta);
}

