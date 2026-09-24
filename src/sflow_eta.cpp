/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
#include"patchBC_interface.h"

sflow_eta::sflow_eta(lexer *p, fdm2D *b , ghostcell *pgc, patchBC_interface *ppBC) : eps(1.0e-6), K(p)
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
    
    wd_criterion=p->A244;
    
    p->Iarray(temp,p->imax*p->jmax);
}

sflow_eta::~sflow_eta()
{
}

void sflow_eta::update(lexer *p, fdm2D *b, ghostcell *pgc, ioflow *pflow, slice &WLout, slice &WL0, slice &WLs, double a)
{
    starttime=pgc->timer();
    
    // continuity: dWL/dt = -div(FE)
    SLICELOOP4
    K(i,j) = 0.0;
    
    SLICELOOP4
    WETDRY
    K(i,j) = -(b->FEx(i,j) - b->FEx(i-1,j))/p->DXN[IP] 
             -(b->FEy(i,j) - b->FEy(i,j-1))/p->DYN[JP]*p->y_dir;
    
    SLICELOOP4
    WLout(i,j) = a*WL0(i,j) + (1.0-a)*(WLs(i,j) + p->dt*K(i,j));
    
    // eta + boundary conditions
    SLICELOOP4
    b->eta(i,j) = WLout(i,j) - b->depth(i,j);
    
    pflow->waterlevel2D(p,b,pgc,b->eta);
    pflow->eta_relax(p,pgc,b->eta);
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    
    SLICELOOP4
    WLout(i,j) = MAX(b->eta(i,j) + b->depth(i,j), 0.0);
    
    pgc->gcsl_start4(p,WLout,gcval_eta);
    
    // wetting & drying
    wetdry(p,b,pgc,WLout);
    
    p->lsmtime+=pgc->timer()-starttime;
}

void sflow_eta::depth_update(lexer *p, fdm2D *b , ghostcell *pgc, slice &WL)
{
	// still water depth
	SLICELOOP4
	b->depth(i,j) = p->wd - b->bed(i,j);
    
    pgc->gcsl_start4(p,b->depth,50);
    
    // the water column is kept: eta follows the bed
    SLICELOOP4
    b->eta(i,j) = WL(i,j) - b->depth(i,j);
    
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
    b->eta(i,j) = wsfout-p->wd;
    }
    
    pgc->gcsl_start4(p,b->eta,gcval_eta);
    
    SLICELOOP4
    WL(i,j) = MAX(b->eta(i,j) + b->depth(i,j), 0.0);
    
    pgc->gcsl_start4(p,WL,gcval_eta);
    
    // wetdry
    wetdry(p,b,pgc,WL);
}

void sflow_eta::ini(lexer *p, fdm2D *b , ghostcell *pgc, ioflow *pflow)
{
    p->phimean=p->F60;
    
    SLICELOOP4
    b->eta(i,j) = 0.0;
    
    pflow->eta_relax(p,pgc,b->eta);
    
    pgc->gcsl_start4(p,b->eta,gcval_eta);
}
