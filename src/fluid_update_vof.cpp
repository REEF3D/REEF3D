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
Author: Hans Bihs, Fabian Knoblauch
--------------------------------------------------------------------*/


#include"density_vof.h"
#include"fluid_update_vof.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

fluid_update_vof::fluid_update_vof(lexer *p, fdm* a, ghostcell* pgc) : dx(p->DXM),
												visc_air(p->W4),visc_water(p->W2),ro_air(p->W3),ro_water(p->W1),visc_body(p->X44)
{
    gcval_ro=1;
	gcval_visc=1;
}

fluid_update_vof::~fluid_update_vof()
{
}

void fluid_update_vof::start(lexer *p, fdm* a, ghostcell* pgc,field& uvel, field& vvel, field& wvel)
{
	double H=0.0, Hro=0.0;
    double H_fb=0.0;
    double psiro;
	p->volume1=0.0;
	p->volume2=0.0;
    
    if(p->count>iter)
    iocheck=0;
	iter=p->count;
    double phival;


	LOOP
	{
        if(p->F92==1)
        {
            phival = a->phi(i,j,k);

            if(phival>p->psi)
                H=1.0;

            if(phival<-p->psi)
                H=0.0;

            if(fabs(phival)<=p->psi)
                H=0.5*(1.0 + phival/(p->psi) + (1.0/PI)*sin((PI*phival)/(p->psi)));
                
            psiro = p->psi;
            if(phival>psiro)
                Hro=1.0;

            if(phival<-psiro)
                Hro=0.0;

            if(fabs(phival)<=psiro)
                Hro=0.5*(1.0 + phival/(psiro) + (1.0/PI)*sin((PI*phival)/(psiro)));
        }
            
        if(p->F92>1)
        {
            H=a->vof(i,j,k);
            if(H>1.0)
                H=1.0;
            if(H<0.0)
                H=0.0;
        }
            
        a->ro(i,j,k)=     ro_water*H +   ro_air*(1.0-H);
        a->visc(i,j,k)= visc_water*H + visc_air*(1.0-H);
            
        p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(H-(1.0-PORVAL4));
        p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*(1.0-H-(1.0-PORVAL4));
	}
    
	pgc->start4(p,a->ro,gcval_ro);
	pgc->start4(p,a->visc,gcval_visc);

	p->volume1 = pgc->globalsum(p->volume1);
	p->volume2 = pgc->globalsum(p->volume2);
    
    if(p->mpirank==0 && iocheck==0 && (p->count%p->P12==0))
    {
	cout<<"Volume 1: "<<p->volume1<<endl;
	cout<<"Volume 2: "<<p->volume2<<endl;
    }
    ++iocheck;

}


int fluid_update_vof::iocheck;
int fluid_update_vof::iter;
