/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"multiphase_fluid_update_rheology.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"rheology_f.h"

multiphase_fluid_update_rheology::multiphase_fluid_update_rheology(lexer *p, fdm* a, ghostcell* pgc) : dx(p->dx),
												visc3(p->W7),visc2(p->W4),visc1(p->W2),ro3(p->W6),ro2(p->W3),ro1(p->W1)
{
    gcval_ro=1;
	gcval_visc=1;
	
	eps12 = p->F321;
	eps13 = p->F322;
	eps23 = p->F323;
    
    prheo = new rheology_f(p,a);
}

multiphase_fluid_update_rheology::~multiphase_fluid_update_rheology()
{
}

void multiphase_fluid_update_rheology::start(lexer *p, fdm* a, ghostcell* pgc, field &ls1, field &ls2)
{
	double H1=0.0;
	double H2=0.0;
	double H3=0.0;
    
	p->volume1=0.0;
	p->volume2=0.0;
	p->volume3=0.0;
    
    if(p->count>iter)
    iocheck=0;
	iter=p->count;
    
    if(p->j_dir==0)        
    epsi = p->F45*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    epsi = p->F45*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);

	LOOP
	{
        
		// water
		if(ls1(i,j,k)>epsi)
		{
		H1=1.0;
        visc1 = prheo->viscosity(p,a,pgc);
		H2=0.0;
		H3=0.0;
		}
		
		if(fabs(ls1(i,j,k))<=epsi)
		{
		H1=0.5*(1.0 + ls1(i,j,k)/epsi + (1.0/PI)*sin((PI*ls1(i,j,k))/epsi));
        visc1 = prheo->viscosity(p,a,pgc);
		
			
			// water-oil
			if(ls2(i,j,k)>epsi)
			{
			H2=0.0;
			H3=1.0-H1;
			}
			
			// water-oil
			if(ls2(i,j,k)<=epsi && ls2(i,j,k)>=0.0)
			{
			H2=0.0;
			H3=1.0-H1;
			}
			
			// water-air
			if(ls2(i,j,k)<-epsi)
			{
			H2=1.0-H1;
			H3=0.0;
			}
			
			// water-air
			if(ls2(i,j,k)>=-epsi && ls2(i,j,k)<0.0)
			{
			H2=1.0-H1;
			H3=0.0;
			}
		}
		
		if(ls1(i,j,k)<-epsi)
		{
		H1=0.0;
        visc1=0.0;
		
			// oil
			if(ls2(i,j,k)>epsi)
			{
			H2=0.0;
			H3=1.0;
			}
			
			// air
			if(ls2(i,j,k)<-epsi)
			{
			H2=1.0;
			H3=0.0;
			}
			
			// oil-air
			if(fabs(ls2(i,j,k))<=epsi)
			{
			H3=0.5*(1.0 + ls2(i,j,k)/epsi + (1.0/PI)*sin((PI*ls2(i,j,k))/epsi));
			H2=1.0-H3;
			}
		}

		a->ro(i,j,k) =     ro1*H1 + ro2*H2 + ro3*H3;
		a->visc(i,j,k) =   visc1*H1 + visc2*H2 + visc3*H3;
        
        p->volume1 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*H1;
        p->volume2 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*H2;
        p->volume3 += p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*H3;
	}
	

	pgc->start4(p,a->ro,gcval_ro);
	pgc->start4(p,a->visc,gcval_visc);

	p->volume1 = pgc->globalsum(p->volume1);
	p->volume2 = pgc->globalsum(p->volume2);
	p->volume3 = pgc->globalsum(p->volume3);

    if(p->mpirank==0 && p->count%p->P12==0)
    {
	cout<<"Volume 1: "<<p->volume1<<endl;
	cout<<"Volume 2: "<<p->volume2<<endl;
	cout<<"Volume 3: "<<p->volume3<<endl;
    }
    ++iocheck;
}


int multiphase_fluid_update_rheology::iocheck;
int multiphase_fluid_update_rheology::iter;


