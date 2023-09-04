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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void iowave::full_initialize(lexer *p, fdm*a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"full NWT initialize"<<endl;
    
    
	FLUIDLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
		
		a->phi(i,j,k) = (wave_h(p,pgc,xg,yg,0.0)-p->pos_z());
	}
	
	UFLUIDLOOP
    {
		xg = xgen1(p);
        yg = ygen1(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
		
		phival = 0.5*(a->phi(i,j,k)+a->phi(i-1,j,k));
        
        if(phival<0.0)
        z = wave_h(p,pgc,xg,yg,0.0)-p->pos_z();
        
        if(phival>=-psi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(phival>=-epsi && phival<-psi)
		H=0.5*(1.0 + phival/fabs(epsi) + (1.0/PI)*sin((PI*phival)/fabs(epsi)));
	
		a->u(i,j,k) = wave_u(p,pgc,xg,yg,z)*H;
	}	
	
	
	VFLUIDLOOP
    {
        xg = xgen2(p);
        yg = ygen2(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
		
		phival = 0.5*(a->phi(i,j,k)+a->phi(i,j-1,k));
        
        if(phival<0.0)
        z = wave_h(p,pgc,xg,yg,0.0)-p->pos_z();
        
        if(phival>=-psi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(phival>=-epsi && phival<-psi)
		H=0.5*(1.0 + phival/fabs(epsi) + (1.0/PI)*sin((PI*phival)/fabs(epsi)));
		
		a->v(i,j,k) = wave_v(p,pgc,xg,yg,z)*H;
	}
	
	WFLUIDLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos3_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos3_z()));
		
		phival = 0.5*(a->phi(i,j,k)+a->phi(i,j,k-1));
        
        if(phival<0.0)
        z = wave_h(p,pgc,xg,yg,0.0)-p->pos_z();
        
        if(phival>=-psi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(phival>=-epsi && phival<-psi)
		H=0.5*(1.0 + phival/fabs(epsi) + (1.0/PI)*sin((PI*phival)/fabs(epsi)));
		
		a->w(i,j,k) = wave_w(p,pgc,xg,yg,z)*H;
	}
	
	FLUIDLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
		
		a->press(i,j,k) = (wave_h(p,pgc,xg,0.0,0.0) - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
	}
	
	
	FLUIDLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p); 
        
		if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
		
		phival = a->phi(i,j,k);

        if(phival>=-psi)
		H=1.0;

		if(phival<-epsi)
		H=0.0;

		if(phival>=-epsi && phival<-psi)
		H=0.5*(1.0 + phival/fabs(epsi) + (1.0/PI)*sin((PI*phival)/fabs(epsi)));
        
        if(p->A10==3)
        H=1.0;
		
		a->Fi(i,j,k) = wave_fi(p,pgc,xg,yg,z)*H;
	}
    
    // eta
	SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);

		a->eta(i,j) = wave_eta(p,pgc,xg,yg);
    }
    
    // Fifsf
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
        z = a->eta(i,j);

		a->Fifsf(i,j) = wave_fi(p,pgc,xg,yg,z);
    }
}

void iowave::full_initialize_ptf(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"full NWT initialize"<<endl;
    
    // eta
	SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);

		a->eta(i,j) = wave_eta(p,pgc,xg,yg);

    }
    
    // Fifsf
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
        z = a->eta(i,j);

		a->Fifsf(i,j) = wave_fi(p,pgc,xg,yg,z);
    }

    
    // Fi
    LOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
        
        
        a->Fi(i,j,k) = wave_fi(p,pgc,xg,yg,z);
      
    }
    

    pgc->start4(p,a->Fi,50);

}