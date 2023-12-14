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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

ioflow_f::ioflow_f(lexer *p, ghostcell *pgc, patchBC_interface *ppBC) 
{
	pBC = ppBC;
    
	walldin_size=1;
	walldout_size=1;
	
	p->Darray(walldin, walldin_size);
    p->Darray(walldout, walldout_size);

	p->Darray(betaB71,p->B71);

	p->Darray(tan_betaB71,p->B71);

	p->Darray(dist_B71,p->B71);
	
	for(n=0;n<p->B71;++n)
	betaB71[n] = (p->B71_b[n]+90.0)*(PI/180.0);
	
	for(n=0;n<p->B71;++n)
	tan_betaB71[n] = tan(betaB71[n]);
	
	
	kinval = 0.0;	
	
	if(p->T10==1 || p->T10==11 || p->T10==21)
	{
    epsval=(pow(0.09,0.75)*pow(kinval,1.5))/(0.5*0.4*p->DXM);
	eddyval = p->cmu*kinval*kinval/epsval;
	}

    if(p->T10==2 || p->T10==12 || p->T10==22)
	{
    epsval=(pow(0.09,0.75)*pow(kinval,0.5))/(0.5*0.4*p->DXM);
	eddyval = kinval/epsval;
	}

    if(p->T10==3 || p->T10==13)
	{
    epsval=(pow(0.09,0.75)*pow(kinval,0.5))/(0.5*0.4*p->DXM);	
	eddyval = kinval/epsval;
	}
    
    if(p->B60==2||p->B60==4)
    {
    hydrograph_in_read(p,pgc);
    p->Ui=hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)/(Ai>1.0e-20?Ai:1.0e20);    
    }
	
	if(p->B60==3||p->B60==4)
    {
    hydrograph_out_read(p,pgc);
    p->Uo=hydrograph_ipol(p,pgc,hydro_out,hydro_out_count)/(Ai>1.0e-20?Ai:1.0e20);    
    }
    
    epsi1 = 1.6;    
    epsi2 = 3.6;
    
    p->fsfoutval=p->fsfout;
    p->fsfinval=p->fsfin;
    
    iter0=-1;
	
}

ioflow_f::~ioflow_f()
{
}
