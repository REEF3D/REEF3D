/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Author: Alexander Hanke
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"partres.h"

#include"lexer.h"
#include"ghostcell.h"
#include"fdm.h"
#include"vrans_f.h"
#include"reinitopo.h"
#include"ioflow.h"
#include"bedshear.h"
#include"sediment_fdm.h"

/// @brief Initializes everything in the sediment for the CFD solver
/// Determines cell which should be filled with particles
/// Allocates memory for the particles
/// Seeds the particles
/// Prepares particles and particle related variables for simulation
/// Initializes VRANS
void sediment_part::ini_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    // vrans
    pvrans->sed_update(p,a,pgc);
    ALOOP
	{
        if(a->topo(i,j,k)<0.0)
	        a->porosity(i,j,k)= p->S24;
        else
	        a->porosity(i,j,k)=1.0;
    }

    pst->setup(p,*a,PP.d50);

    volume0 = pst->volume(p,*a,PP);
    volume0 = pgc->globalsum(volume0);
    volume = volume0;


    double h;
    ILOOP
    JLOOP
	{
		KLOOP
		PBASECHECK
		{
        if(a->topo(i,j,k-1)<0.0 && a->topo(i,j,k)>0.0)
        h = -(a->topo(i,j,k-1)*p->DZP[KP])/(a->topo(i,j,k)-a->topo(i,j,k-1)) + p->pos_z()-p->DZP[KP];
		}
		s.bedzh(i,j)=h;
        s.bedzh0(i,j)=h;
	}
    /*
    k=0;
    SLICEBASELOOP
    s.bedzh0(i,j) = s.bedzh(i,j) = p->ZN[KP] + 0.5*p->DZP[KP]-a->topo(i,j,k);*/
    
    pgc->gcsl_start4(p,s.bedzh0,50); 

    if(p->I40!=1)
    {
        // seed
        seed_ini(p,a,pgc);
        PP.reserve(maxparticle);
        pst->setParticleMax(maxparticle);
        seed(p,a);
        // make_stationary(p,a,&PP);
    }

    gparticle_active = pgc->globalisum(PP.size);

    fill_PQ_cfd(p,a,pgc);
    // pstmove(p,*a,*pgc,PP,s,*pturb);
    if(p->Q13==1)
        pst->update(p,*a,*pgc,PP);
    
    // print
    if((p->I40!=1)||(p->I40==1&&inicount>0))
    print_particles(p);
    
    if(p->mpirank==0)
    {
        if(p->I40!=1)
            cout<<"Sediment particles: "<<gparticle_active<<"\n";
        else if (inicount>0)
            cout<<"Loaded particles "<<gparticle_active<<" from state file.\n";
        cout<<"Initial bed volume: "<<std::setprecision(prec)<<volume0<<" m^3"<<endl;
    }
    
    SLICELOOP4
    s.bedk(i,j)=0;
    
    SLICELOOP4
    {
        KLOOP
            PBASECHECK
            if(a->topo(i,j,k)<0.0 && a->topo(i,j,k+1)>=0.0)
                s.bedk(i,j)=k+1;
        s.reduce(i,j)=0.3;
    }
    pbedshear.taubed(p,a,pgc,&s);
    pgc->gcsl_start4(p,s.tau_eff,1);
    pbedshear.taucritbed(p,a,pgc,&s);
    pgc->gcsl_start4(p,s.tau_crit,1);
    
    ++inicount;
    debug(p,a,pgc);

    // if(inicount==1)
    // {
    //     particles_obj Dummy(1000);
    //     seedDummy(p,a,Dummy);
    //     printDummyVTP(p,Dummy);
    // }
    
    pgc->gcdf_update(p,a);
}