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

#include"bedconc.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedconc::bedconc(lexer *p)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    shields=p->S30;
    visc=p->W2;
    kappa=0.4;
    ks=2.5*d50;
    adist=0.5*d50;
    deltab=3.0*d50;
    Rstar=(rhosed-rhowat)/rhowat;
}

bedconc::~bedconc()
{
}

void bedconc::start(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    SLICELOOP4
    s->cbn(i,j) = s->cbe(i,j);
    
    // cb* van Rijn
    SLICELOOP4
    {
	
    Ti=MAX((s->tau_eff(i,j)-s->tau_crit(i,j))/s->tau_crit(i,j),0.0);


    Ds= d50*pow((Rstar*g)/(visc*visc),1.0/3.0);
    
    if(s->active(i,j)==1)
    s->cbe(i,j) = (0.015*d50*pow(Ti,1.5))/(pow(Ds,0.3)*adist);
    
    if(s->active(i,j)==0)
    s->cbe(i,j) = 0.0;
    }
    
    // cb at first cell center
    SLICELOOP4
    {
        k=s->bedk(i,j);
        

        zdist = (p->ZP[KP]-s->bedzh(i,j));

        s->cb(i,j) = s->cbe(i,j)*pow(((s->waterlevel(i,j)-zdist)/zdist)*(adist/(s->waterlevel(i,j)-adist)),zdist);
        
        //cout<<"CB: "<<s->cbe(i,j)<<" "<<s->cb(i,j)<<" "<<zdist<<" "<<adist<<" "<<s->waterlevel(i,j)<<endl;
    }
    
    if(p->S34==2)
    {
    SLICELOOP4
    s->qbe(i,j) += s->conc(i,j);
    
    pgc->gcsl_start4(p,s->qbe,1);
    }
}




